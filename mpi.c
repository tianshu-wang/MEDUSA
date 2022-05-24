#include "decs.h"

#if (USE_MPI==TRUE)
void ND_FUNC(build_types,NDIM)();

static MPI_Datatype *neighbor_send_type0, *neighbor_send_type1;
static MPI_Datatype *neighbor_recv_type0, *neighbor_recv_type1;
static MPI_Datatype *neighbor_send_shock_type, *neighbor_recv_shock_type;
static MPI_Datatype mpi_tracer_type;
//static proc_info *proc;
static int64_t *neighbor_send_offsets;
static int64_t *neighbor_recv_offsets;
static MPI_Comm axis_comm;
static int n_axis_comm;
static MPI_Request *mpi_tracer_requests;
static MPI_Status  *mpi_tracer_status;
static int n_tracer_neighbors,*tracer_neighbors;

#endif

void init_mpi_setup()
{
  int dd,ideal=1,actual=1,worst=1;
  double efficiency;

  #if (USE_MPI==FALSE)

  myrank   = 0;
  numprocs = 1;
  iorank   = 0;
  ND_ARRAY_SET(global_grid_dims,n1,n2,n3);
  ND_ARRAY_SET(istart,0,0,0);
  ND_ARRAY_SET(istop,n1,n2,n3);
  ND_ARRAY_SET(my_grid_dims,n1,n2,n3);

  #else // WE'RE USING MPI, SET IT UP!

  ND_ARRAY_SET(global_grid_dims,n1,n2,n3);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

  proc = split_domain(&ideal,&actual,&worst);
  efficiency = 1.0 - (actual-ideal)/(1.0*ideal);
  if (mpi_io_proc()) fprintf(stderr,"done!  cells/rank: %d (ideal)  %d (actual)  => efficiency = %g\n",ideal,actual,efficiency);

  // if (efficiency<0.50) {
  //   if (mpi_io_proc()) {
  //     fprintf(stderr,"[init_mpi_setup]:  Inefficient decomposition!\n");
  //     fprintf(stderr,"worst = %d\n",worst);
  //   }
  //   exit(1);
  // }

  iorank = myrank;

  DLOOP {
    istart[dd]       = proc[myrank].istart[dd];
    istop[dd]        = proc[myrank].istop[dd];
    my_grid_dims[dd] = istop[dd] - istart[dd];
  }

  for (dd=NDIM; dd<SPACEDIM; dd++) {
    istart[dd]       = 0;
    istop[dd]        = 0;
    my_grid_dims[dd] = 0;
  }

  mpi_requests = malloc_rank1(4*proc[myrank].nneigh, sizeof *mpi_requests);
  mpi_status   = malloc_rank1(4*proc[myrank].nneigh, sizeof *mpi_status);


  // save the ranks of all my potential tracer neighbors
  n_tracer_neighbors = 0;
  tracer_neighbors = NULL;
  // add direct face neighbors first as an optimization for later searches
  for(int i=0;i<proc[myrank].nneigh;i++) {
      int irank = proc[myrank].neigh_id[i];
      tracer_neighbors = realloc(tracer_neighbors, (n_tracer_neighbors+1)*sizeof(int));
      tracer_neighbors[n_tracer_neighbors] = irank;
      n_tracer_neighbors++;
  }
  for(int i=0;i<proc[myrank].nneigh;i++) {
      int irank = proc[myrank].neigh_id[i];
      #if(NDIM > 1)
      int idir = abs(proc[myrank].neigh_dir[i]);
      for(int j=0;j<proc[irank].nneigh;j++) {
          int jrank = proc[irank].neigh_id[j];
          int jdir = abs(proc[irank].neigh_dir[j]);
          if(idir==jdir) continue; // skip neighbors that are 0 or 2 shifts in a particular direction
          tracer_neighbors = realloc(tracer_neighbors, (n_tracer_neighbors+1)*sizeof(int));
          tracer_neighbors[n_tracer_neighbors] = jrank;
          n_tracer_neighbors++;
      }
      #endif
  }



  mpi_tracer_requests = malloc_rank1(n_tracer_neighbors, sizeof(MPI_Request));
  mpi_tracer_status   = malloc_rank1(n_tracer_neighbors, sizeof(MPI_Status));

  // make MPI derived type to send tracers
  tracer_mpi_packet_t packet;
  int packet_count = 2;
  int blocklengths[] = {2*NDIM+1, 1};//, 1};
  MPI_Aint begin,end;
  MPI_Get_address(&packet.d[0], &begin);
  MPI_Get_address(&packet.id, &end);
  MPI_Aint indices[] = {0, end-begin};
  //fprintf(stderr,"***\n\n\n  indices[] = { %ld , %ld }\n\n\n***", indices[0], indices[1]);
  MPI_Datatype old_types[] = {MPI_DOUBLE, MPI_INT};//, MPI_CHAR};
  MPI_Type_create_struct(packet_count, blocklengths, indices, old_types, &mpi_tracer_type);
  MPI_Type_commit(&mpi_tracer_type);


  #endif /* USE_MPI */

  if (mpi_io_proc()) fprintf(stderr,"numprocs = %d\n",numprocs);
}


int is_outside_rank(int index[], int rank)
{
#if (USE_MPI==TRUE)
    if(index[0] < proc[rank].istart[0] || index[0] >= proc[rank].istop[0]) return 1;
    #if(NDIM>1)
    if(index[1] < JS(index[0],proc[rank].istart[1]) || index[1] >= JS(index[0],proc[rank].istop[1])) return 2;
    #endif
    #if(NDIM==3)
    if(index[2] < KS(index[0],index[1],proc[rank].istart[2]) || index[2] >= KS(index[0],index[1],proc[rank].istop[2])) return 3;
    #endif
#endif
    return 0;
}

int x_to_rank(double x[], int rank_guess) {

    int index[NDIM];

    x_to_ijk(x,index);

    // first check if the guess is right
    if(!is_outside_rank(index, rank_guess)) return rank_guess;

#if (USE_MPI==TRUE)
    for(int i=0;i<n_tracer_neighbors;i++) {
        if(!is_outside_rank(index, tracer_neighbors[i])) return tracer_neighbors[i];
    }
#endif

#if (NDIM==1)
    fprintf(stderr,"problem coming %d: %d\n", myrank, index[0]);
#elif (NDIM==2)
    fprintf(stderr,"problem coming %d: %d %d\n", myrank, index[0], index[1]*DJS(index[0]));
#else
    fprintf(stderr,"problem coming %d: %d %d %d\n", myrank, index[0], index[1]*DJS(index[0]), index[2]*DKS(index[0],index[1]));
#endif


    return -1;
}


#if (NDIM==1)
void build_types1()
{
  #if (USE_MPI==TRUE)
  int i,count,blocks;
  MPI_Aint base0,base1,offsets0,offsets1;
  proc_info* myproc = &proc[myrank];
  proc_info* nproc = NULL;

  neighbor_send_type0 = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_type0);
  neighbor_recv_type0 = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_type0);
  neighbor_send_type1 = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_type1);
  neighbor_recv_type1 = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_type1);

  for (i=0; i < myproc->nneigh; i++) {
    nproc = &proc[myproc->neigh_id[i]];

    count    = 1;
    blocks   = NG*nvars;
    offsets0 = 0;
    offsets1 = 0;

    MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,istart[0],jj,kk,0), &base0);
    MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,istart[0],jj,kk,0), &base1);

    if (myproc->istart[0] < nproc->istart[0]) {
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,istart[0]+my_grid_dims[0]-NG,jj,kk,0), &offsets0);
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,istart[0]+my_grid_dims[0]-NG,jj,kk,0), &offsets1);
      offsets0 -= base0;
      offsets1 -= base1;
      MPI_Type_create_hindexed(count, &blocks, &offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
      MPI_Type_create_hindexed(count, &blocks, &offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,istart[0]+my_grid_dims[0],jj,kk,0), &offsets0);
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,istart[0]+my_grid_dims[0],jj,kk,0), &offsets1);
      offsets0 -= base0;
      offsets1 -= base1;
      MPI_Type_create_hindexed(count, &blocks, &offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
      MPI_Type_create_hindexed(count, &blocks, &offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
    } else {
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,istart[0],jj,kk,0), &offsets0);
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,istart[0],jj,kk,0), &offsets1);
      offsets0 -= base0;
      offsets1 -= base1;
      MPI_Type_create_hindexed(count, &blocks, &offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
      MPI_Type_create_hindexed(count, &blocks, &offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,istart[0]-NG,jj,kk,0), &offsets0);
      MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,istart[0]-NG,jj,kk,0), &offsets1);
      offsets0 -= base0;
      offsets1 -= base1;
      MPI_Type_create_hindexed(count, &blocks, &offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
      MPI_Type_create_hindexed(count, &blocks, &offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
    }

    MPI_Type_commit(&neighbor_send_type0[i]);
    MPI_Type_commit(&neighbor_send_type1[i]);
    MPI_Type_commit(&neighbor_recv_type0[i]);
    MPI_Type_commit(&neighbor_recv_type1[i]);
  }
  #endif

  return;
}
#endif


#if (NDIM==2)
void build_types2()
{
  #if (USE_MPI==TRUE)
  int i,j,idir,ip,c,ix,iy;
  int i0,j0,k0;
  MPI_Aint base0,base1,*offsets0,*offsets1,base_shock,*offsets_shock;
  int overlap_start,overlap_stop,count,*blocks,*blocks_shock;
  proc_info* myproc = &proc[myrank];
  proc_info* nproc = NULL;

  neighbor_send_type0      = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_type0);
  neighbor_recv_type0      = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_type0);
  neighbor_send_type1      = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_type1);
  neighbor_recv_type1      = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_type1);
  neighbor_send_shock_type = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_shock_type);
  neighbor_recv_shock_type = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_shock_type);

  for (i=0; i<myproc->nneigh; i++) {

    // processor ip in in direction idir-1 and has info nproc
    idir  = myproc->neigh_dir[i];
    ip    = myproc->neigh_id[i];
    nproc = &proc[ip];

    switch (abs(idir)) {  // remember idir is the direction PLUS 1; sign encodes above/below, left/right, etc.

      case 1: // neighbor is in x0-direction
        // build 0-type
        // first figure out overlap in 1-direction
        overlap_start = MAX(myproc->istart[1],nproc->istart[1]);
        overlap_stop  = MIN(myproc->istop[1] ,nproc->istop[1] );

        count         = NG;
        blocks        = malloc_rank1(count, sizeof *blocks);
        offsets0      = malloc_rank1(count, sizeof *offsets0);
        offsets1      = malloc_rank1(count, sizeof *offsets1);
        blocks_shock  = malloc_rank1(count, sizeof *blocks_shock);
        offsets_shock = malloc_rank1(count, sizeof *offsets_shock);

        i0 = istart[0];
        j0 = istart[1]/DJS(i0);
        k0 = istart[2]/DKS(i0,j0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,i0,j0,k0,0), &base0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,i0,j0,k0,0), &base1);
        MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,i0,j0,k0), &base_shock);

        // Send types
        for (c=0; c<count; c++) {
          ix = (idir > 0) ? myproc->istop[0] - NG + c : myproc->istart[0] + c;
          iy = JS(ix,overlap_start);
          blocks[c] = (JS(ix,overlap_stop) - JS(ix,overlap_start)) * nvars;
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,0,0), &offsets0[c]);
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,0,0), &offsets1[c]);
          offsets0[c] -= base0;
          offsets1[c] -= base1;
          blocks_shock[c] = JS(ix,overlap_stop) - JS(ix,overlap_start);
          MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,0), &offsets_shock[c]);
          offsets_shock[c] -= base_shock;
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
        MPI_Type_commit(&neighbor_send_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
        MPI_Type_commit(&neighbor_send_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_send_shock_type[i]);
        MPI_Type_commit(&neighbor_send_shock_type[i]);

        // Recv types
        for (c=0; c<count; c++) {
          ix = (idir > 0) ? myproc->istop[0] + c : myproc->istart[0] - NG + c;
          iy = JS(ix,overlap_start);
          blocks[c] = (JS(ix,overlap_stop) - JS(ix,overlap_start)) * nvars;
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,0,0), &offsets0[c]);
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,0,0), &offsets1[c]);
          offsets0[c] -= base0;
          offsets1[c] -= base1;
          blocks_shock[c] = JS(ix,overlap_stop) - JS(ix,overlap_start);
          MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,0), &offsets_shock[c]);
          offsets_shock[c] -= base_shock;
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
        MPI_Type_commit(&neighbor_recv_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
        MPI_Type_commit(&neighbor_recv_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_recv_shock_type[i]);
        MPI_Type_commit(&neighbor_recv_shock_type[i]);

        free(blocks);
        free(offsets0);
        free(offsets1);
        free(blocks_shock);
        free(offsets_shock);
        break;

      case 2: // neighbor is in x1-direction
        // build 1-type, this should be easy
        // figure out overlap in 0-direction
        overlap_start = MAX(myproc->istart[0],nproc->istart[0]);
        overlap_stop  = MIN(myproc->istop[0] ,nproc->istop[0] );

        count         = overlap_stop - overlap_start;
        blocks        = malloc_rank1(count, sizeof *blocks);
        offsets0      = malloc_rank1(count, sizeof *offsets0);
        offsets1      = malloc_rank1(count, sizeof *offsets1);
        blocks_shock  = malloc_rank1(count, sizeof *blocks_shock);
        offsets_shock = malloc_rank1(count, sizeof *offsets_shock);

        i0 = istart[0];
        j0 = istart[1]/DJS(i0);
        k0 = istart[2]/DKS(i0,j0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,i0,j0,k0,0), &base0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,i0,j0,k0,0), &base1);
        MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,i0,j0,k0), &base_shock);

        // Send types
        for (c=0; c<count; c++) {
          // y is the fastest-changing index, hence the need for small blocks...
          ix = overlap_start + c;
          iy = (idir > 0) ? JS(ix,istop[1]) - NG : JS(ix,istart[1]);
          blocks[c] = NG * nvars;
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,0,0), &offsets0[c]);
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,0,0), &offsets1[c]);
          offsets0[c] -= base0;
          offsets1[c] -= base1;
          blocks_shock[c] = NG;
          MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,0), &offsets_shock[c]);
          offsets_shock[c] -= base_shock;
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
        MPI_Type_commit(&neighbor_send_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
        MPI_Type_commit(&neighbor_send_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_send_shock_type[i]);
        MPI_Type_commit(&neighbor_send_shock_type[i]);

        // Recv types
        for (c=0; c<count; c++) {
          ix = overlap_start + c;
          iy = (idir > 0) ? JS(ix,istop[1]) : JS(ix,istart[1]) - NG;
          blocks[c] = NG * nvars;
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,0,0), &offsets0[c]);
          MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,0,0), &offsets1[c]);
          offsets0[c] -= base0;
          offsets1[c] -= base1;
          blocks_shock[c] = NG;
          MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,0), &offsets_shock[c]);
          offsets_shock[c] -= base_shock;
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
        MPI_Type_commit(&neighbor_recv_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
        MPI_Type_commit(&neighbor_recv_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_recv_shock_type[i]);
        MPI_Type_commit(&neighbor_recv_shock_type[i]);

        free(blocks);
        free(offsets0);
        free(offsets1);
        free(blocks_shock);
        free(offsets_shock);
        break;

      default:
        // should never be here
        fprintf(stderr,"Invalid idir in build_types2\n");
        exit(3);
    }  /* end switch */
  }  /* end for */
  #endif  /* USE_MPI==TRUE */
}
#endif


#if (NDIM==3)
void build_types3()
{
  #if (USE_MPI==TRUE)
  int i,idir,ip,ii,jj,kk,ix,iy,iz;
  int i0,j0,k0;
  MPI_Aint base0,base1,*offsets0,*offsets1,*offsets_shock,base_shock;
  int overlap_start[NDIM],overlap_stop[NDIM],count,*blocks,c,*blocks_shock;
  proc_info* myproc = &proc[myrank];
  proc_info* nproc = NULL;

  neighbor_send_type0      = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_type0);
  neighbor_recv_type0      = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_type0);
  neighbor_send_type1      = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_type1);
  neighbor_recv_type1      = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_type1);
  neighbor_send_shock_type = malloc_rank1(myproc->nneigh, sizeof *neighbor_send_shock_type);
  neighbor_recv_shock_type = malloc_rank1(myproc->nneigh, sizeof *neighbor_recv_shock_type);

  for (i=0; i<myproc->nneigh; i++) {

    // processor ip in in direction idir-1 and has info nproc
    idir  = myproc->neigh_dir[i];
    ip    = myproc->neigh_id[i];
    nproc = &proc[ip];

    switch (abs(idir)) {  // remember idir is the direction times +/- 1; sign encodes above/below, left/right, etc.

      case 1: // neighbor is in x0-direction
        // build 0-type
        // first figure out overlap in 1/2-directions
        overlap_start[1] = MAX(myproc->istart[1],nproc->istart[1]);
        overlap_stop[1]  = MIN(myproc->istop[1] ,nproc->istop[1] );

        overlap_start[2] = MAX(myproc->istart[2],nproc->istart[2]);
        overlap_stop[2]  = MIN(myproc->istop[2] ,nproc->istop[2] );

        #if (NDIM==3 && GEOM==SPHERICAL)
        // If procs are neighbors in phi across the origin, do something slightly different here
        if (idir < 0 && myproc->istart[0]==0  && nproc->istart[0]==0) {
          overlap_start[2] = MAX(myproc->istart[2],(nproc->istart[2]+n3/2)%n3);
          overlap_stop[2]  = MIN(myproc->istop[2],(nproc->istop[2]+n3/2-1)%n3+1);
        }
        #endif

        i0 = istart[0];
        j0 = JS(i0,istart[1]);
        k0 = KS(i0,j0,istart[2]);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,i0,j0,k0,0), &base0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,i0,j0,k0,0), &base1);
        MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,i0,j0,k0), &base_shock);

        // Send types
        count = 0;
        for (ii=0; ii<NG; ii++) {
          ix = (idir > 0) ? myproc->istop[0] + ii - NG : myproc->istart[0] + ii;
          count += JS(ix,overlap_stop[1]) - JS(ix,overlap_start[1]);
        }
        blocks        = malloc_rank1(count, sizeof *blocks);
        offsets0      = malloc_rank1(count, sizeof *offsets0);
        offsets1      = malloc_rank1(count, sizeof *offsets1);
        blocks_shock  = malloc_rank1(count, sizeof *blocks_shock);
        offsets_shock = malloc_rank1(count, sizeof *offsets_shock);

        c = 0;
        for (ii=0; ii<NG; ii++) {
          ix = (idir > 0) ? myproc->istop[0] + ii - NG : myproc->istart[0] + ii;
          for (iy = JS(ix,overlap_start[1]); iy < JS(ix,overlap_stop[1]); iy++) {
            iz = KS(ix,iy,overlap_start[2]);
            blocks[c] = (KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2])) * nvars;
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,iz,0), &offsets0[c]);
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,iz,0), &offsets1[c]);
            offsets0[c] -= base0;
            offsets1[c] -= base1;
            blocks_shock[c] = KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2]);
            MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,iz), &offsets_shock[c]);
            offsets_shock[c] -= base_shock;
            c++;
          }
        }

        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
        MPI_Type_commit(&neighbor_send_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
        MPI_Type_commit(&neighbor_send_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_send_shock_type[i]);
        MPI_Type_commit(&neighbor_send_shock_type[i]);

        free(blocks);
        free(offsets0);
        free(offsets1);
        free(blocks_shock);
        free(offsets_shock);

        // Recv types
        count = 0;
        for (ii=0; ii<NG; ii++) {
          int ix = (idir > 0) ? myproc->istop[0] + ii : myproc->istart[0] + ii - NG;
          count += JS(ix,overlap_stop[1]) - JS(ix,overlap_start[1]);
        }
        blocks        = malloc_rank1(count, sizeof *blocks);
        offsets0      = malloc_rank1(count, sizeof *offsets0);
        offsets1      = malloc_rank1(count, sizeof *offsets1);
        blocks_shock  = malloc_rank1(count, sizeof *blocks_shock);
        offsets_shock = malloc_rank1(count, sizeof *offsets_shock);

        c = 0;
        for (ii=0; ii<NG; ii++) {
          ix = (idir > 0) ? myproc->istop[0] + ii : myproc->istart[0] + ii - NG;
          #if (NDIM==3 && GEOM==SPHERICAL)
          // If procs are neighbors in phi across the origin, reflect iy in theta
          if (idir < 0 && myproc->istart[0]==0  && nproc->istart[0]==0) {
            ix = myproc->istart[0] - ii - 1;
          }
          #endif
          for (iy = JS(ix,overlap_start[1]); iy < JS(ix,overlap_stop[1]); iy++) {
            #if (NDIM==3 && GEOM==SPHERICAL)
            // If procs are neighbors in phi across the origin, reflect iy in theta
            if (idir < 0 && myproc->istart[0]==0  && nproc->istart[0]==0) {
              iy = my_grid_dims[1]/DJS(ix) - 1 - (iy - JS(ix,istart[1])) + JS(ix,istart[1]);
            }
            #endif
            iz = KS(ix,iy,overlap_start[2]);
            blocks[c] = (KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2])) * nvars;
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,iz,0), &offsets0[c]);
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,iz,0), &offsets1[c]);
            offsets0[c] -= base0;
            offsets1[c] -= base1;
            blocks_shock[c] = KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2]);
            MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,iz), &offsets_shock[c]);
            offsets_shock[c] -= base_shock;
            c++;
            #if (NDIM==3 && GEOM==SPHERICAL)
            // If procs are neighbors in phi across the origin, reflect iy in theta
            if (idir < 0 && myproc->istart[0]==0  && nproc->istart[0]==0) {
              iy = my_grid_dims[1]/DJS(ix) - 1 - (iy - JS(ix,istart[1])) + JS(ix,istart[1]);
            }
            #endif
          }
        }

        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
        MPI_Type_commit(&neighbor_recv_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
        MPI_Type_commit(&neighbor_recv_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_recv_shock_type[i]);
        MPI_Type_commit(&neighbor_recv_shock_type[i]);

        free(blocks);
        free(offsets0);
        free(offsets1);
        free(blocks_shock);
        free(offsets_shock);
        break;

      case 2: // neighbor is in x1-direction
        // build 1-type
        // figure out overlap in 0/2-directions
        overlap_start[0] = MAX(myproc->istart[0],nproc->istart[0]);
        overlap_stop[0]  = MIN(myproc->istop[0] ,nproc->istop[0] );

        overlap_start[2] = MAX(myproc->istart[2],nproc->istart[2]);
        overlap_stop[2]  = MIN(myproc->istop[2] ,nproc->istop[2] );

        #if (NDIM==3 && GEOM==SPHERICAL)
        // If procs are neighbors in phi across the polar axis, do something slightly different here
        if ((idir < 0 && myproc->istart[1]==0  && nproc->istart[1]==0 )
         || (idir > 0 && myproc->istop[1] ==n2 && nproc->istop[1] ==n2)) {
          overlap_start[2] = MAX(myproc->istart[2],(nproc->istart[2]+n3/2)%n3);
          overlap_stop[2]  = MIN(myproc->istop[2],(nproc->istop[2]+n3/2-1)%n3+1);
        }
        #endif

        count         = (overlap_stop[0]-overlap_start[0])*NG;
        blocks        = malloc_rank1(count, sizeof *blocks);
        offsets0      = malloc_rank1(count, sizeof *offsets0);
        offsets1      = malloc_rank1(count, sizeof *offsets1);
        blocks_shock  = malloc_rank1(count, sizeof *blocks_shock);
        offsets_shock = malloc_rank1(count, sizeof *offsets_shock);

        i0 = istart[0];
        j0 = JS(i0,istart[1]);
        k0 = KS(i0,j0,istart[2]);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,i0,j0,k0,0), &base0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,i0,j0,k0,0), &base1);
        MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,i0,j0,k0), &base_shock);

        // Send types
        c = 0;
        for (ix=overlap_start[0]; ix<overlap_stop[0]; ix++) {
          for (jj=0; jj<NG; jj++) {
            iy = (idir > 0) ? JS(ix,istop[1]) + jj - NG : JS(ix,istart[1]) + jj;
            #if (NDIM==3 && GEOM==SPHERICAL)
            // If procs are neighbors in phi across the polar axis, do something slightly different here
            if (idir < 0 && myproc->istart[1]==0  && nproc->istart[1]==0 ) {
              iy = JS(ix,istart[1]) - jj + NG - 1;
            }
            if (idir > 0 && myproc->istop[1] ==n2 && nproc->istop[1] ==n2) {
              iy = JS(ix,istop[1]) - jj - 1;
            }
            #endif
            iz = KS(ix,iy,overlap_start[2]);
            blocks[c] = (KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2])) * nvars;
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,iz,0), &offsets0[c]);
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,iz,0), &offsets1[c]);
            offsets0[c] -= base0;
            offsets1[c] -= base1;
            blocks_shock[c] = KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2]);
            MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,iz), &offsets_shock[c]);
            offsets_shock[c] -= base_shock;
            c++;
          }
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
        MPI_Type_commit(&neighbor_send_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
        MPI_Type_commit(&neighbor_send_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_send_shock_type[i]);
        MPI_Type_commit(&neighbor_send_shock_type[i]);

        // Recv types
        c = 0;
        for (ix=overlap_start[0]; ix<overlap_stop[0]; ix++) {
          for (jj=0; jj<NG; jj++) {
            iy = (idir > 0) ? JS(ix,istop[1]) + jj : JS(ix,istart[1]) + jj - NG;
            iz = KS(ix,iy,overlap_start[2]);
            blocks[c] = (KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2])) * nvars;
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,iz,0), &offsets0[c]);
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,iz,0), &offsets1[c]);
            offsets0[c] -= base0;
            offsets1[c] -= base1;
            blocks_shock[c] = KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2]);
            MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,iz), &offsets_shock[c]);
            offsets_shock[c] -= base_shock;
            c++;
          }
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
        MPI_Type_commit(&neighbor_recv_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
        MPI_Type_commit(&neighbor_recv_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_recv_shock_type[i]);
        MPI_Type_commit(&neighbor_recv_shock_type[i]);

        free(blocks);
        free(offsets0);
        free(offsets1);
        free(blocks_shock);
        free(offsets_shock);
        break;

      case 3: // neighbor is in x2-direction
        // build 2-type
        // figure out overlap in 0/1-directions
        overlap_start[0] = MAX(myproc->istart[0],nproc->istart[0]);
        overlap_stop[0]  = MIN(myproc->istop[0] ,nproc->istop[0] );

        overlap_start[1] = MAX(myproc->istart[1],nproc->istart[1]);
        overlap_stop[1]  = MIN(myproc->istop[1] ,nproc->istop[1] );

        i0 = istart[0];
        j0 = JS(i0,istart[1]);
        k0 = KS(i0,j0,istart[2]);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,i0,j0,k0,0), &base0);
        MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,i0,j0,k0,0), &base1);
        MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,i0,j0,k0), &base_shock);

        count = 0;
        for (ix = overlap_start[0]; ix < overlap_stop[0]; ix++) {
          count += JS(ix,overlap_stop[1]) - JS(ix,overlap_start[1]);
        }
        blocks        = malloc_rank1(count, sizeof *blocks);
        offsets0      = malloc_rank1(count, sizeof *offsets0);
        offsets1      = malloc_rank1(count, sizeof *offsets1);
        blocks_shock  = malloc_rank1(count, sizeof *blocks_shock);
        offsets_shock = malloc_rank1(count, sizeof *offsets_shock);

        // Send types
        c = 0;
        for (ix = overlap_start[0]; ix < overlap_stop[0]; ix++) {
          for (iy = JS(ix,overlap_start[1]); iy < JS(ix,overlap_stop[1]); iy++) {
            iz = (idir > 0) ? KS(ix,iy,istop[2]) - NG : KS(ix,iy,istart[2]);
            blocks[c] = NG * nvars;
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,iz,0), &offsets0[c]);
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,iz,0), &offsets1[c]);
            offsets0[c] -= base0;
            offsets1[c] -= base1;
            blocks_shock[c] = NG;
            MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag,ix,iy,iz), &offsets_shock[c]);
            offsets_shock[c] -= base_shock;
            c++;
          }
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_send_type0[i]);
        MPI_Type_commit(&neighbor_send_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_send_type1[i]);
        MPI_Type_commit(&neighbor_send_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_send_shock_type[i]);
        MPI_Type_commit(&neighbor_send_shock_type[i]);

        // Recv types
        c = 0;
        for (ix = overlap_start[0]; ix < overlap_stop[0]; ix++) {
          for (iy = JS(ix,overlap_start[1]); iy < JS(ix,overlap_stop[1]); iy++) {
            iz = (idir > 0) ? KS(ix,iy,istop[2]) : KS(ix,iy,istart[2]) - NG;
            blocks[c] = NG * nvars;
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_p ,ix,iy,iz,0), &offsets0[c]);
            MPI_Get_address(&NDP_ELEM_LINEAR(sim_ph,ix,iy,iz,0), &offsets1[c]);
            offsets0[c] -= base0;
            offsets1[c] -= base1;
            blocks_shock[c] = NG;
            MPI_Get_address(&ND_ELEM_LINEAR(sim_shock_flag, ix,iy,iz), &offsets_shock[c]);
            offsets_shock[c] -= base_shock;
            c++;
          }
        }
        MPI_Type_create_hindexed(count, blocks, offsets0, MPI_DOUBLE, &neighbor_recv_type0[i]);
        MPI_Type_commit(&neighbor_recv_type0[i]);
        MPI_Type_create_hindexed(count, blocks, offsets1, MPI_DOUBLE, &neighbor_recv_type1[i]);
        MPI_Type_commit(&neighbor_recv_type1[i]);
        MPI_Type_create_hindexed(count, blocks_shock, offsets_shock, MPI_INT, &neighbor_recv_shock_type[i]);
        MPI_Type_commit(&neighbor_recv_shock_type[i]);

        free(blocks);
        free(offsets0);
        free(offsets1);
        free(blocks_shock);
        free(offsets_shock);
        break;

      default:
        // should never be here
        fprintf(stderr,"Invalid idir in build_types2\n");
        exit(3);
    }
  }
  #endif /* USE_MPI==TRUE */
}
#endif /* NDIM==3 */


void init_consistent_tracer_id()
{
#if (USE_MPI==TRUE)
    if(consistent_tracer_ids) return;

    unsigned int n_accum, n_send;
    MPI_Status stat;
    n_accum = 0;
    if(myrank == 0 && numprocs > 1) {
        MPI_Send(&n_tracer_current, 1, MPI_UNSIGNED, 1, 0, MPI_COMM_WORLD);
    }
    for(int i = 1; i < numprocs-1; i++) {
        if(myrank == i) {
            MPI_Recv(&n_accum, 1, MPI_UNSIGNED, myrank-1, 0, MPI_COMM_WORLD, &stat);
            n_send = n_accum + n_tracer_current;
            MPI_Send(&n_send, 1, MPI_UNSIGNED, myrank+1, 0, MPI_COMM_WORLD);
        }
    }
    if(myrank == numprocs-1 && numprocs > 1) {
        MPI_Recv(&n_accum, 1, MPI_UNSIGNED, myrank-1, 0, MPI_COMM_WORLD, &stat);
    }

    tracer_t *tr = tracers;
    while(tr != NULL) {
        tr->id += n_accum;
        tr = tr->next;
    }
#endif
}

static tracer_mpi_packet_t **tagged_tracers;
static int n_per_tracer;
static int send_tracers_outstanding=0;
void mpi_send_tracers(int stage)
{
#if (USE_MPI==TRUE)
    static int firstc = 1;
    static int *nbuff, *nsend;
    int dd;
    tracer_t *tr = tracers;
    tracer_t *prev = NULL;
    tracer_t *next = NULL;
    n_per_tracer = NDIM+1;
    if(firstc) {
        nbuff = malloc_rank1(n_tracer_neighbors, sizeof(int));
        nsend = malloc_rank1(n_tracer_neighbors, sizeof(int));
        firstc = 0;
    }
    tagged_tracers = malloc_rank1(n_tracer_neighbors, sizeof(tracer_mpi_packet_t *));

    for(int i = 0; i < n_tracer_neighbors; i++) {
        nbuff[i] = 10;
        nsend[i] = 0;
        tagged_tracers[i] = malloc_rank1(nbuff[i], sizeof(tracer_mpi_packet_t));
    }

    int removed = 0;
    while(tr != NULL) {
        if(!tr->active) {
            tr = tr->next;
            continue;
        }
        int mpi_rank;
        if(stage == 0) {
            mpi_rank = x_to_rank(tr->xh, myrank);
        } else {
            mpi_rank = x_to_rank(tr->x, myrank);
        }
        if(mpi_rank == -1) { // tracer has left the building
            fprintf(stderr,"TRACER PROBLEM!!!   stage = %d\n", stage);
            // dump it's current history

            // and remove it
            /*next = tr->next;
            if(prev != NULL) {
                prev->next = next;
            } else {
                tracers = next;
            }
            //if(tr->hist != NULL) free(tr->hist);
            free(tr);
            tr = next;
            n_tracer_current--;
            */
            if(stage==0) {
              #if (NDIM==3)
                fprintf(stderr,"%g %g %g\n", tr->xh[0], tr->xh[1], tr->xh[2]);
              #endif
              #if (NDIM==2)
                fprintf(stderr,"%g %g\n", tr->xh[0], tr->xh[1]);
              #endif
              #if (NDIM==1)
                fprintf(stderr,"%g\n", tr->xh[0]);
              #endif
            } else {
              #if (NDIM==3)
                fprintf(stderr,"%g %g %g\n", tr->x[0], tr->x[1], tr->x[2]);
              #endif
              #if (NDIM==2)
                fprintf(stderr,"%g %g\n", tr->x[0], tr->x[1]);
              #endif
              #if (NDIM==1)
                fprintf(stderr,"%g\n", tr->x[0]);
              #endif
            }

            prev = tr;
            tr->active = 0;
            tr = tr->next;
        } else if(mpi_rank != myrank) { // deal with a tracer that has moved off myrank
            // which neighbor is this?
            int ineigh;
            for(ineigh=0;ineigh<n_tracer_neighbors;ineigh++) if(tracer_neighbors[ineigh]==mpi_rank) break;

            // now copy it for sending
            if(nsend[ineigh] >= nbuff[ineigh]) {
                nbuff[ineigh] += 10;
                tagged_tracers[ineigh] = realloc(tagged_tracers[ineigh], nbuff[ineigh]*sizeof(tracer_mpi_packet_t));
            }
            int iadd = nsend[ineigh];
            DLOOP {
                tagged_tracers[ineigh][iadd].d[dd] = tr->x[dd];
                tagged_tracers[ineigh][iadd].d[NDIM+dd] = tr->xh[dd];
            }
            tagged_tracers[ineigh][iadd].d[2*NDIM] = tr->mass;
            tagged_tracers[ineigh][iadd].id = tr->id;
            //tagged_tracers[ineigh][iadd].active = tr->active;
            nsend[ineigh]++;

            // and remove it
            next = tr->next;
            if(prev != NULL) {
                prev->next = next;
            } else {
                tracers = next;
            }
            //if(tr->hist != NULL) free(tr->hist);
//            if(myrank == 0 && mpi_rank == 1) {
//                fprintf(stderr,"!!!!!\n");
//                fprintf(stderr,"%d %d\n", ineigh, iadd);
//                fprintf(stderr,"%g %g %d\n", tr->x[0], tr->mass, tr->id);
//            }
            free(tr);
            tr = next;
            n_tracer_current--;
            removed++;

        } else {
            prev = tr;
            tr = tr->next;
        }
    }

/*    if(myrank == 0 && nsend[0] > 0) {
        tr = tracers;
        while(tr != NULL) {
            fprintf(stderr,"%g %g %d\n", tr->x[0], tr->mass, tr->id);
            tr = tr->next;
        }
    }
*/
    // send the data
    for(int i = 0; i < n_tracer_neighbors; i++) {
        int err = MPI_Isend(tagged_tracers[i], nsend[i], mpi_tracer_type, tracer_neighbors[i], tracer_neighbors[i]+995, MPI_COMM_WORLD, &mpi_tracer_requests[i]);
        if(err != MPI_SUCCESS) {
            fprintf(stderr,"MPI_Isend error %d\n", err);
        }
    }
    send_tracers_outstanding = 1;

    //mpi_recv_tracers();
#endif
}

void mpi_recv_tracers()
{
#if (USE_MPI==TRUE)
    if(!send_tracers_outstanding) return;
    int is_mesg,count,dd;
    int received = 0;
    tracer_mpi_packet_t *recv_buff;
    MPI_Status tracer_status;

    while(received < n_tracer_neighbors) {
        MPI_Iprobe(MPI_ANY_SOURCE, myrank+995, MPI_COMM_WORLD, &is_mesg, &tracer_status);
        if(is_mesg) {
            MPI_Get_count(&tracer_status, mpi_tracer_type, &count);
            //if(count > 0) fprintf(stderr,"proc %d received %d tracers from %d\n", myrank, count, tracer_status.MPI_SOURCE);
            recv_buff = malloc_rank1(count, sizeof(tracer_mpi_packet_t));
            int source = tracer_status.MPI_SOURCE;
            MPI_Recv(recv_buff, count, mpi_tracer_type, source, myrank+995, MPI_COMM_WORLD, &tracer_status);
/*            if(myrank > tracer_status.MPI_SOURCE) {
                for(int i = 0; i < count; i++) {
                    DLOOP fprintf(stderr,"%d %g\n", i, recv_buff[i].d[dd]);
                    fprintf(stderr,"%d %g\n", i, recv_buff[i].d[NDIM]);
                    fprintf(stderr,"%d %d\n", i, recv_buff[i].id);
                }
            }
*/
            for(int i = 0; i < count; i++) {
                tracer_t *head = tracers;
                tracers = malloc_rank1(1, sizeof(tracer_t));

                DLOOP tracers->x[dd] = recv_buff[i].d[dd];
                DLOOP tracers->xh[dd] = recv_buff[i].d[NDIM+dd];
                tracers->mass = recv_buff[i].d[2*NDIM];
                tracers->id = recv_buff[i].id;
                tracers->active = 1;
                if(tracers->id < 0 || tracers->id >= n_tracer_global) {
                    fprintf(stderr,"tracer id out of bounds!!!   %d\n", tracers->id);
                }
                tracers->next = head;
            }
            n_tracer_current += count;
            if(count > 0 && recv_buff != NULL) free(recv_buff);
            received++;
        }
    }

    MPI_Waitall(n_tracer_neighbors, mpi_tracer_requests, mpi_tracer_status);
    //fprintf(stderr,"proc %d has %u\n", myrank, n_tracer_current);
    send_tracers_outstanding = 0;
    for(int i = 0; i < n_tracer_neighbors; i++) free(tagged_tracers[i]);
    free(tagged_tracers);
#endif
}

#define TEST_VAR 0

#if (USE_LINEAR_ALLOCATION==TRUE)
void sync_mpi_boundaries(double ** p)
#else
void sync_mpi_boundaries(double NDP_PTR p)
#endif
{
  #if (USE_MPI==TRUE)
  int stride,i;
  int i0,j0,k0;
  proc_info *myproc = &proc[myrank];

  i0 = istart[0];
  j0 = JS(i0,istart[1]);
  k0 = KS(i0,j0,istart[2]);

  if (half_step_sync) {
    #if (NDIM==1)
    stride = 2;
    #else
    stride = 4;
    #endif

    for (i=0; i<myproc->nneigh; i++) {
      MPI_Isend(&NDP_ELEM_LINEAR(p,i0,j0,k0,0),
        1, neighbor_send_type1[i], myproc->neigh_id[i], 0, MPI_COMM_WORLD, &mpi_requests[stride*i]);
      MPI_Irecv(&NDP_ELEM_LINEAR(p,i0,j0,k0,0),
        1, neighbor_recv_type1[i], myproc->neigh_id[i], 0, MPI_COMM_WORLD, &mpi_requests[stride*i+1]);
      #if (NDIM>1)
      MPI_Isend(&ND_ELEM_LINEAR(sim_shock_flag, i0,j0,k0),
        1, neighbor_send_shock_type[i], myproc->neigh_id[i], 0, MPI_COMM_WORLD, &mpi_requests[stride*i+2]);
      MPI_Irecv(&ND_ELEM_LINEAR(sim_shock_flag, i0,j0,k0),
        1, neighbor_recv_shock_type[i], myproc->neigh_id[i], 0, MPI_COMM_WORLD, &mpi_requests[stride*i+3]);
      #endif
    }

  } else {

    for (i=0; i<myproc->nneigh; i++) {
      MPI_Isend(&NDP_ELEM_LINEAR(p,i0,j0,k0,0),
        1, neighbor_send_type0[i], myproc->neigh_id[i], 0, MPI_COMM_WORLD, &mpi_requests[2*i]);
      MPI_Irecv(&NDP_ELEM_LINEAR(p,i0,j0,k0,0),
        1, neighbor_recv_type0[i], myproc->neigh_id[i], 0, MPI_COMM_WORLD, &mpi_requests[2*i+1]);
    }

  }
  #endif /* USE_MPI==TRUE */
}

void complete_mpi_communication(int half_step)
{
  int i,nperproc;

  #if (USE_MPI==TRUE)
  if (proc[myrank].nneigh > 0) {
    #if (NDIM>1)
    nperproc = (half_step) ? 4 : 2;
    #else
    nperproc = 2;
    #endif

    MPI_Waitall(nperproc*proc[myrank].nneigh, mpi_requests, mpi_status);
    for (i=0; i<nperproc*proc[myrank].nneigh; i++) {
      if (mpi_status[i].MPI_ERROR != MPI_SUCCESS) {
        fprintf(stdout,"MPI request[%d] failed\n", i);
        fflush(stdout);
        exit(1);
      }
    }
  }
  #endif
}


void complete_mpi_bc(int half_step) {
  int i,nperproc;
  #if (USE_MPI==TRUE)
  if (proc[myrank].nneigh>0){
    #if (GEOM==SPHERICAL && NDIM==3)
    // Apply spherical origin boundary conditions where applicable
    int j,k,vv,jstart,jstop;
    int self_periodic = (istart[2] == 0 && istop[2] == global_grid_dims[2]);
#if (USE_LINEAR_ALLOCATION==TRUE)
    double ** p = (half_step) ? sim_ph : sim_p;
#else
    double NDP_PTR p = (half_step) ? sim_ph : sim_p;
#endif

    GPU_PRAGMA(omp target teams distribute parallel for)
    VLOOP {
      // Only do this if at origin and NOT self-periodic
      if (istart[0]==0 && !self_periodic && bc[vv].lo[0]==SPHERICAL_ORIGIN) {
        for (i=istart[0]-NG; i<istart[0]; i++) {
          #if (NDIM>1)
          JSLOOP(i,j) {
          #endif
            #if (NDIM==3)
            KSLOOP(i,j,k) {
            #endif
              if (vv==U1 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
              #if (DO_RADIATION)
              if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==0 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
              #endif
            #if (NDIM==3)
            }
            #endif
          #if (NDIM>1)
          }
          #endif
        }
      }

      // Only do this if at the north pole and NOT self-periodic
      if (istart[1]==0 && !self_periodic && bc[vv].lo[1]==SPHERICAL_ORIGIN) {
        ISLOOP(i) {
          jstart = JS(i,istart[1]);
          for (j=jstart-NG; j<jstart; j++) {
            #if (NDIM==3)
            KSLOOP(i,j,k) {
            #endif
              if (vv==U2 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
              #if (DO_RADIATION)
              if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==1 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
              #endif
            #if (NDIM==3)
            }
            #endif
          }
        }
      }

      // Only do this if at the south pole and NOT self-periodic
      if (istop[1]==global_grid_dims[1] && !self_periodic && bc[vv].hi[1]==SPHERICAL_ORIGIN) {
        ISLOOP(i) {
          jstop = JS(i,istop[1]);
          for (j=jstop; j<jstop+NG; j++) {
            #if (NDIM==3)
            KSLOOP(i,j,k) {
            #endif
              if (vv==U2 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
              #if (DO_RADIATION)
              if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==1 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
              #endif
            #if (NDIM==3)
            }
            #endif
          }
        }
      }
    }
    #endif
  }
  #endif
}

double mpi_min(double val)
{
  double min_val;

  #if (USE_MPI==TRUE)
  MPI_Allreduce(&val, &min_val, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  #else
  min_val = val;
  #endif

  return min_val;
}

double mpi_max(double val)
{
  double max_val;

  #if (USE_MPI==TRUE)
  MPI_Allreduce(&val, &max_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  #else
  max_val = val;
  #endif

  return max_val;
}

void mpi_global_reduce(double *arr, int len)
{
  int i;
  double *glob_arr;

  #if (USE_MPI==TRUE)
  glob_arr = malloc_rank1(len, sizeof *glob_arr);
  MPI_Allreduce(arr, glob_arr, len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (i=0; i<len; i++) arr[i] = glob_arr[i];
  free(glob_arr);
  #endif

  return;
}

double **iglob_arr;
int iglob_max_tag;
#if (USE_MPI==TRUE)
static MPI_Request *iglob_reduce_req;
#endif


void mpi_global_reduce_uint(unsigned int *arr, int len)
{
  int i;
  unsigned int *glob_arr;

  #if (USE_MPI==TRUE)
  glob_arr = malloc_rank1(len, sizeof *glob_arr);
  MPI_Allreduce(arr, glob_arr, len, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  for (i=0; i<len; i++) arr[i] = glob_arr[i];
  free(glob_arr);
  #endif

  return;
}

void mpi_global_ireduce_start(double *arr, int len, int tag)
{
  #if (USE_MPI==TRUE)
  if (tag == 0) {
    iglob_arr        = malloc_rank1(1,   sizeof *iglob_arr       );
    iglob_arr[0]     = malloc_rank1(len, sizeof *iglob_arr[0]    );
    iglob_reduce_req = malloc_rank1(1,   sizeof *iglob_reduce_req);
    iglob_max_tag    = 0;
  } else {
    iglob_arr        = realloc(iglob_arr, (tag+1)*sizeof *iglob_arr);
    iglob_arr[tag]   = malloc_rank1(len, sizeof *iglob_arr[tag]);
    iglob_reduce_req = realloc(iglob_reduce_req, (tag+1)*sizeof *iglob_reduce_req);
    iglob_max_tag    = tag;
  }
  MPI_Iallreduce(arr, iglob_arr[tag], len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &iglob_reduce_req[tag]);
  #endif
}

void mpi_global_ireduce_finish(double *arr, int len, int tag)
{
  #if (USE_MPI==TRUE)
  MPI_Status status;
  int err = MPI_Wait(&iglob_reduce_req[tag], &status);
  if (err != MPI_SUCCESS) {
    fprintf(stderr,"rank %d failed on MPI_Wait in mpi_global_ireduce_finish with error %d\n", myrank, err);
    fflush(stderr);
    exit(2);
  }
  memcpy(arr, iglob_arr[tag], len*sizeof(double));
  free(iglob_arr[tag]);
  if (tag == iglob_max_tag) {
    free(iglob_arr);
    free(iglob_reduce_req);
  }
  #endif
}

void mpi_bcast(double *arr, int len)
{
  #if (USE_MPI==TRUE)
  MPI_Bcast(arr, len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
}

int mpi_io_proc()
{
  #if (USE_MPI==TRUE)
  if (myrank==0) return 1;
  else return 0;
  #else
  return 1;
  #endif
}
