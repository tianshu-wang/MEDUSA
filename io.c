#include <err.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#ifndef u_short
typedef unsigned short int u_short;
#endif
#include <fts.h>

#include <unistd.h>

#include "decs.h"
#include "constants.h"
#include <hdf5.h>

static int dump_id=0;
static int pdump_id=0;
extern double Eerr;

#define RADIATION_NORM  (1e25)

#if (OUTPUT_PRECISION==SINGLE)
typedef float  Real;
#define HDF5PREC H5T_NATIVE_FLOAT
#elif (OUTPUT_PRECISION==DOUBLE)
typedef double Real;
#define HDF5PREC H5T_NATIVE_DOUBLE
#else
#error Must set OUTPUT_PRECISION=SINGLE or OUTPUT_PRECISION=DOUBLE in Makefile!
#endif

void run_log()
{
  fprintf(stderr,"# step = %d, dt = %1.5e, t = %1.5e, wallclock = %1.5e", istep, dt, t, telapsed);
  #if (DO_RADIATION==TRUE)
  fprintf(stderr,", implicit_iter = %d, Eerr = %1.2e", implicit_iter, Eerr);
  #endif
  if (detect_tbounce) fprintf(stderr,", central rho = %1.2e",rho_cent_max);
  if (use_chat) fprintf(stderr,", chat/CLIGHT = %1.2e",chat/CLIGHT);
  fprintf(stderr,"\n");

  implicit_iter = 0;

  return;
}

void output_timers() {
#if DO_TIMING==FALSE
  return;
#else
  static int firstc = 1;
  char const * timers_fname = "timing.csv";
  FILE * fp = NULL;
  double times[MAX_NTIMERS];
  double tmax[MAX_NTIMERS];
  double tmin[MAX_NTIMERS];
  double tavg[MAX_NTIMERS];

  if (firstc) {
    if (mpi_io_proc()) {
      if (0 != access(timers_fname, F_OK)) {
        fprintf(stderr, "[output_timers]:  Creating %s\n", timers_fname);
        fp = fopen(timers_fname, "w");
        fprintf(fp, "timer,time,step,ncycle,max,min,avg\n");
      } else {
        fprintf(stderr,"[output_timers]:  %s exists, appending...\n",timers_fname);
        fp = fopen(timers_fname, "a");
      }
    }
    firstc = 0;
  } else {
    if (mpi_io_proc()) {
      fp = fopen(timers_fname, "a");
    }
  }

  for (int tidx = 0; tidx < ntimers; ++tidx) {
    Timer * timer = get_timer(tidx);
    timer_get_total(timer, &times[tidx]);
  }
#if(USE_MPI==TRUE)
  MPI_Reduce(&times[0], &tmax[0], ntimers, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&times[0], &tmin[0], ntimers, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&times[0], &tavg[0], ntimers, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  for(int itemp=0;itemp<ntimers;itemp++){
    tmax[itemp] = times[0];
    tmin[itemp] = times[0];
    tavg[itemp] = times[0];
  }
#endif

  if (mpi_io_proc()) {
    for (int tidx = 0; tidx < ntimers; ++tidx) {
      Timer * timer = get_timer(tidx);
      fprintf(fp, "%s,%f,%d,%lu,%f,%f,%f\n",
          timer->name, t, istep, timer->ncycles, tmax[tidx], tmin[tidx],
          tavg[tidx]/numprocs);
    }
    fclose(fp);
  }
#endif // TIMING==FALSE
}

#if (NDIM==1)
/* Order of dumped variables:  x[0],RHO,UU,U1,U2,U3,{YE},{Er[0],F1[0],...},PRESS,CS,TEMP,ENT,GAMMA,dM,gr_lapse,gr_lapse_edge */
void dump_ascii1(int is_restart)
{
  int ii,jj=0,kk=0,vv,dd,g;
  double x[SPACEDIM];
  double *scale;
  double dm,mass;
  char name[256];
  char* str=NULL;
  FILE *fp;
  #if (USE_MPI)
  MPI_Status status;
  #endif

  TIMER_START("dump_ascii1");

  //GPU_PRAGMA(omp target update to(sim_p[:cell_count_all][:nvars]))
  //update_eos(eostable,sim_p,sim_eos);
  //GPU_PRAGMA(omp target update from(sim_eos[:cell_count_all][:NEOS]))

  scale = malloc_rank1(nvars, sizeof *scale);

  #if (PRINT_OPAC_DATA==TRUE)
  double *kappa = malloc_rank1(ngroups, sizeof *kappa);
  double *jg    = malloc_rank1(ngroups, sizeof *jg   );
  double *sc    = malloc_rank1(ngroups, sizeof *sc   );
  double *delta = malloc_rank1(ngroups, sizeof *delta);
  #endif

  if (is_restart) {
    if (mpi_io_proc()) mkdir("restarts", 0777);
    str = strrchr(model_name,'/');
    if (str)
      str++;
    else
      str = model_name;
    // sprintf(name,"restarts/%s.%d.%dg.%dms.restart_1d",model_name,n1,nr1,(int)(1000*perturb_delay));
    sprintf(name,"restarts/%s.%d.%dg.%dms.restart_1d",str,n1,nr1,(int)(1000*perturb_delay));
  } else {
    sprintf(name,"dumps/dump_%05d", dump_id);
  }
  for (int ioproc=0; ioproc<numprocs; ioproc++) {
    if (myrank==ioproc) {
      if (ioproc==0) {
        fp = fopen(name,"w");
      } else {
        fp = fopen(name,"a");
      }
      if (fp==NULL) {
        fprintf(stderr,"Failed to open %s, aborting\n", name);
        exit(1);
      }

      if (ioproc==0) {
        fprintf(fp,"# %g\n", t);
        mass = 0.0;
      } else {
        #if (USE_MPI)
        MPI_Recv(&mass, 1, MPI_DOUBLE, ioproc-1, 0, MPI_COMM_WORLD, &status);
        #endif
      }
      for (ii=istart[0]; ii<istop[0]; ii++) {
        ijk_to_x(ii,jj,kk,x);
        #if (GEOM==SPHERICAL)
        x[0] = r_of_x(rx_info,x[0]);
        #endif
        VLOOP { scale[vv] = 1.0; }
        DLOOP {
          fprintf(fp, "%16.14e ", x[dd]);
          scale[U1+dd] = ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][dd];
        }
        #if (DO_RADIATION==TRUE)
        for (vv=ifrad1; vv<ifrad1+ngroups; vv++) {
          scale[vv] = ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0];
        }
        #endif
        VLOOP { fprintf(fp, "%16.14e ", NDP_ELEM_LINEAR(sim_p,ii,jj,kk,vv)*scale[vv]); }
        #if (PRINT_OPAC_DATA==TRUE)
        opac_emis(NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO), NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP), NDP_ELEM_LINEAR(sim_p,ii,jj,kk,YE), kappa, jg, sc, delta, NULL, NULL, NULL);
        GLOOP { fprintf(fp, "%16.14e %16.14e %16.14e %16.14e ", kappa[g], jg[g], sc[g], delta[g]); }
        #endif

        EOSLOOP { fprintf(fp, "%16.15e ", NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,vv)); }
        dm = ND_ELEM_LINEAR(geom,ii,jj,kk).volume*NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
        fprintf(fp, "%16.15e ", mass+0.5*dm);
        mass += dm;

        #if (GR_MONOPOLE==TRUE)
        fprintf(fp, "%16.15e %16.15e ", ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0],ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1]);
        #else
        fprintf(fp, "%16.15e %16.15e ", 1.0,1.0);
        #endif

        fprintf(fp,"\n");
      }
      fclose(fp);
    }
    #if (USE_MPI==TRUE)
    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank==ioproc && ioproc<numprocs-1) {
      MPI_Send(&mass, 1, MPI_DOUBLE, ioproc+1, 0, MPI_COMM_WORLD);
    }
    #endif
  }

  if (!is_restart) dump_id++;

  #if (PRINT_OPAC_DATA==TRUE)
  free(kappa);
  free(jg);
  free(sc);
  free(delta);
  #endif

  free(scale);
  TIMER_STOP;
  return;
}
#endif


/*#################################
  #################################

    Generic XDMF output functions

  #################################
  #################################*/

void write_xml_closing(FILE *xml)
{
  fprintf(xml, "   </Grid>\n");
  fprintf(xml, " </Domain>\n");
  fprintf(xml, "</Xdmf>\n");

  fclose(xml);
}

void write_scalar(Real *data, hid_t file_id, char *name, hid_t filespace, hid_t memspace, FILE *xml)
{
  char fname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if (sname != NULL) sname++;
  else sname = fname;

  if (mpi_io_proc()) {
    fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", name);
    #if (NDIM==1)
    fprintf(xml, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", n1,sizeof(Real));
    #elif (NDIM==2)
    fprintf(xml, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", n1,n2,sizeof(Real));
    #elif (NDIM==3)
    fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", n1,n2,n3,sizeof(Real));
    #endif
    fprintf(xml, "        %s:%s\n", sname, name);
    fprintf(xml, "       </DataItem>\n");
    fprintf(xml, "     </Attribute>\n");
    fflush(xml);
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, HDF5PREC, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  herr_t status = H5Dwrite(dset_id, HDF5PREC, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  return;
}

void write_vector(Real *data, hid_t file_id, char *name, hid_t filespace, hid_t memspace, FILE *xml)
{
  // it would be nice to use xdmf Vector Attribute, but visit doesn't seem able to handle it properly

  int dd,total_cells;
  char fname[80];
  char cname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if (sname != NULL) sname++;
  else sname = fname;

  total_cells = 1;
  DLOOP { total_cells *= my_grid_dims[dd]; }

  DLOOP {
    sprintf(cname, "%s%d", name, dd);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, cname, HDF5PREC, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    herr_t status = H5Dwrite(dset_id, HDF5PREC, memspace, filespace, plist_id, &data[dd*total_cells]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    if (mpi_io_proc()) {
      fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", cname);
      #if (NDIM==1)
      fprintf(xml, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", n1,sizeof(Real));
      #elif (NDIM==2)
      fprintf(xml, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", n1,n2,sizeof(Real));
      #elif (NDIM==3)
      fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"%zu\" Format=\"HDF\">\n", n1,n2,n3,sizeof(Real));
      #endif
      fprintf(xml, "        %s:%s\n", sname, cname);
      fprintf(xml, "       </DataItem>\n");
      fprintf(xml, "     </Attribute>\n");
    }
  }

  return;
}

#if (NDIM==3)
#if (USE_LINEAR_ALLOCATION==TRUE)
void read_scalar(hid_t file_id, char const *name, int vv, double ** p)
#else
void read_scalar(hid_t file_id, char const *name, int vv, double NDP_PTR p)
#endif
{
  int i = 0, j = 0, k = 0;
  int ii = 0, jj = 0, kk = 0;
  int ind = 0;

  double * data = malloc(n1*n2*n3*sizeof(double));

  hid_t dset_id = H5Dopen(file_id, name, H5P_DEFAULT);
  if (dset_id < 0) {
    fprintf(stdout, "Could not open dataset \"%s\"\n", name);
    fflush(stdout);
    abort();
  }
  herr_t ierr = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (ierr < 0) {
    fprintf(stdout, "Could not read dataset \"%s\"\n", name);
    fflush(stdout);
    abort();
  }
  H5Dclose(dset_id);

  for (i = istart[0]; i < istop[0]; ++i)
  for (j = istart[1]; j < istop[1]; ++j)
  for (k = istart[2]; k < istop[2]; ++k) {
    ii = i;
    jj = JS(i,j);
    kk = KS(i,jj,k);
    ind = k + n3*(j + i*n2);
    NDP_ELEM_LINEAR(p,ii,jj,kk,vv) = data[ind];
  }

  free(data);
}
#endif

void dump_grid()
{
  int i=0,j=0,k=0,dd,jj,kk,ind,total_size;
  int nkp,kp,kkk;
  float *x[SPACEDIM];
  int *mask;
  double xp[NDIM];
  double scale;
  FILE *fp;
  hid_t dset_id,plist_id,file_id,filespace,memspace;
  herr_t status;
  hsize_t file_start[NDIM], file_count[NDIM];
  hsize_t mem_start[NDIM+1];
  hsize_t fdims[] = {n1+1, n2+1, n3+1};
  hsize_t dims[] = {my_grid_dims[0], my_grid_dims[1], my_grid_dims[2]};
  hsize_t rdims[1];
  const char *coordNames[] = {"/X", "/Y", "/Z"};
  char name[20];

  TIMER_START("dump_grid");

  if (mpi_io_proc())
    mkdir("dumps", 0777);

  DLOOP {
    if (istop[dd] == global_grid_dims[dd]) dims[dd]++;
  }

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  #if (USE_MPI==TRUE)
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  file_id = H5Fcreate("dumps/grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  filespace = H5Screate_simple(NDIM, fdims, NULL);
  DLOOP {
    file_start[dd] = istart[dd];
    file_count[dd] = dims[dd];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(NDIM, dims, NULL);
  DLOOP { mem_start[dd] = 0; }
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  total_size = 1;
  DLOOP { total_size *= dims[dd]; }
  SLOOP {
    x[dd] = malloc_rank1(total_size, sizeof(float));
  }
  total_size = 1;
  DLOOP { total_size *= my_grid_dims[dd]; }
  mask = malloc_rank1(total_size, sizeof(int));

  // first do the coordinates
  if (strcmp(grid_units, "km") == 0) {
    scale = 1.0e-5;
  } else {
    scale = 1.0;
  }

  ind = 0;
  for (i=istart[0]; i<istart[0]+dims[0]; i++) {
    #if (NDIM>1)
    for (j=istart[1]; j<istart[1]+dims[1]; j++) {
    #endif
      #if (NDIM==3)
      for (k=istart[2]; k<istart[2]+dims[2]; k++) {
      #endif
        ijk_to_Cart(i,j,k,xp);
        x[0][ind] = xp[0]*scale;
        #if (NDIM>1)
        x[1][ind] = xp[1]*scale;
        #if (NDIM==3)
        x[2][ind] = xp[2]*scale;
        #endif
        #endif
        ind++;
      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  DLOOP {
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, coordNames[dd], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x[dd]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }

  // now write areas
  ind = 0;
  for (i=istart[0]; i<istart[0]+dims[0] ;i++) {
    #if (NDIM>1)
    for (j=istart[1]; j<istart[1]+dims[1]; j++) {
    #endif
      jj = JS(i,j);
      #if (NDIM==3)
      for (k=istart[2]; k<istart[2]+dims[2]; k++) {
      #endif
        kk = KS(i,jj,k);
        x[0][ind] = ND_ELEM_LINEAR(geom,i,jj,kk).area[0];
        #if (NDIM>1)
        nkp = DKS(i,jj)/MIN(DKS(i,jj),DKS(i,jj-1));
        kp = kk*nkp;
        x[1][ind] = 0.0;
        for (kkk=kp; kkk<kp+nkp; kkk++) {
          x[1][ind] += ND_ELEM_LINEAR(geom,i,jj,kp).area[1];
        }
        #if (NDIM==3)
        x[2][ind] = ND_ELEM_LINEAR(geom,i,jj,kk).area[2];
        #endif
        #endif
        ind++;
      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }

  DLOOP {
    sprintf(name, "area%d", dd);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x[dd]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }

  H5Sclose(filespace);
  H5Sclose(memspace);
  // done with writing edge-based grid data


  // write "centered" grid data
  DLOOP {
    fdims[dd]--;
    dims[dd] = my_grid_dims[dd];
  }

  filespace = H5Screate_simple(NDIM, fdims, NULL);
  DLOOP {
    file_start[dd] = istart[dd];
    file_count[dd] = dims[dd];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(NDIM, dims, NULL);
  DLOOP { mem_start[dd] = 0; }
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  // now write volumes
  ind = 0;
  for (i=istart[0]; i<istart[0]+dims[0]; i++) {
    #if (NDIM>1)
    for (j=istart[1]; j<istart[1]+dims[1]; j++) {
    #endif
      jj = JS(i,j);
      #if (NDIM==3)
      for (k=istart[2]; k<istart[2]+dims[2]; k++) {
      #endif
        kk = KS(i,jj,k);
        x[0][ind++] = ND_ELEM_LINEAR(geom,i,jj,kk).volume;
      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }
  strcpy(name,"volume");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x[0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  // now write metric
  ind = 0;
  for (i=istart[0]; i<istart[0]+dims[0]; i++) {
    #if (NDIM>1)
    for (j=istart[1]; j<istart[1]+dims[1]; j++) {
    #endif
      jj = JS(i,j);
      #if (NDIM==3)
      for (k=istart[2]; k<istart[2]+dims[2]; k++) {
      #endif
        kk = KS(i,jj,k);
        x[0][ind] = ND_ELEM_LINEAR(geom,i,jj,kk).gcov[0][0];
        x[1][ind] = ND_ELEM_LINEAR(geom,i,jj,kk).gcov[0][1];
        x[2][ind] = ND_ELEM_LINEAR(geom,i,jj,kk).gcov[0][2];
        ind++;
      #if (NDIM==3)
      }
      #endif
    #if (NDIM>1)
    }
    #endif
  }
  SLOOP {
    sprintf(name, "gcov%d", dd);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x[dd]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }

  // now write a mask to identify non-active cells
  #if (NDIM>1)
  ind = 0;
  for (i=istart[0]; i<istart[0]+dims[0]; i++) {
    for (j=istart[1]; j<istart[1]+dims[1]; j++) {
      jj = JS(i,j);
      #if (NDIM==3)
      for (k=istart[2]; k<istart[2]+dims[2]; k++) {
      #endif
        mask[ind++] = ((j%DJS(i) == 0) && (k%DKS(i,jj) == 0));
      #if (NDIM==3)
      }
      #endif
    }
  }
  strcpy(name,"mask");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, mask);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  #endif /* NDIM>1 */

  H5Sclose(filespace);
  H5Sclose(memspace);

  SLOOP free(x[dd]);

  #if (DO_RADIATION)
  rdims[0] = 1;
  filespace = H5Screate_simple(1, rdims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? 1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(1, rdims, NULL);
  mem_start[0] = 0;
  file_count[0] = mpi_io_proc() ? 1 : 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  strcpy(name,"nr1");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &nr1);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  strcpy(name,"nr2");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &nr2);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  strcpy(name,"nr3");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &nr3);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(filespace);
  H5Sclose(memspace);

  rdims[0] = ngroups;
  filespace = H5Screate_simple(1, rdims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? ngroups : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(1, rdims, NULL);
  mem_start[0] = 0;
  file_count[0] = mpi_io_proc() ? ngroups : 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  strcpy(name,"egroup");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, egroup);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  strcpy(name,"degroup");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, degroup);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  #endif  /* DO_RADIATION */

  H5Fclose(file_id);
  TIMER_STOP;
  return;
}

void dump_tracers()
{
    static int firstc = 1;
    float *data;
    int64_t *pid;
    int *idata;
    char name[80];
    int index[NDIM];
    hid_t file_id,filespace,memspace,plist_id,group_id,dset_id;
    hid_t vec_filespace, vec_memspace;
    hsize_t mem_start[NDIM+1], file_start[NDIM+1], one, zero;
    hsize_t file_count[NDIM+1], file_grid_dims[NDIM+1];
    hsize_t *f_coord;
    herr_t status;

    mpi_recv_tracers();
    //fprintf(stderr,"proc %d dumping %d tracers\n", myrank, n_tracer_current);

    data = malloc_rank1(n_tracer_current, sizeof *data);
    f_coord = malloc_rank1(n_tracer_current, sizeof(hsize_t)) ;
    pid = malloc_rank1(n_tracer_current,sizeof(int64_t));

    if(firstc) {
        firstc = 0;
    }

    zero = 0;
    one = 1;

    sprintf(name,"dumps/part_%05d.h5part", pdump_id);

    #if (USE_MPI==TRUE)
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    #else
    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    #endif

    file_grid_dims[0] = 1;
    filespace = H5Screate_simple(1, file_grid_dims, NULL);
    file_start[0] = 0;
    file_count[0] = (mpi_io_proc()) ? one : zero;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Sclose(memspace);

    group_id = H5Gcreate(file_id, "/Step#0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);


    file_grid_dims[0] = n_tracer_global;
    hsize_t ntrace = n_tracer_current;
    filespace = H5Screate_simple(1, file_grid_dims, NULL);
    file_start[0] = 0;
    file_count[0] = n_tracer_current;
    //H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, (const hsize_t *)&f_coord);
    //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    memspace  = H5Screate_simple(1, &ntrace, NULL);

    int i=0;
    tracer_t *tr = tracers;
    while(tr != NULL) {
        f_coord[i] = tr->id;
        pid[i] = tr->id;
        data[i] = r_of_x(rx_info,tr->x[0]);
        #if(NDIM>1)
        data[i] *= sin(th_of_x(thx_info,tr->x[1]));
        #endif
        #if(NDIM==3)
        data[i] *= cos(tr->x[2]);
        #endif
        data[i] /= 1.e5;
        i++;
        tr = tr->next;
    }


    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/id", H5T_NATIVE_INT64, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_INT64, memspace, filespace, plist_id, pid);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/x", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    i=0;
    tr = tracers;
    while(tr != NULL) {
        data[i] = tr->active;
        i++;
        tr = tr->next;
    }
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/active", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    #if(NDIM>1)
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/y", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    tr = tracers;
    i = 0;
    while(tr != NULL) {
        #if(NDIM==2)
        data[i] = r_of_x(rx_info,tr->x[0])*cos(th_of_x(thx_info,tr->x[1]));
        #else
        data[i] = r_of_x(rx_info,tr->x[0])*sin(th_of_x(thx_info,tr->x[1]))*sin(tr->x[2]);
        #endif
        data[i] /= 1.e5;
        i++;
        tr = tr->next;
    }
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
    #endif

    #if(NDIM==3)
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/z", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    tr = tracers;
    i = 0;
    while(tr != NULL) {
        data[i] = r_of_x(rx_info,tr->x[0])*cos(th_of_x(thx_info,tr->x[1]));
        data[i] /= 1.e5;
        i++;
        tr = tr->next;
    }
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
    #endif


    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/mass", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    tr = tracers;
    i = 0;
    while(tr != NULL) {
        data[i] = tr->mass;
        i++;
        tr = tr->next;
    }
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    // now fill in useful quantities to keep track of

    i=0;
    tr = tracers;
    while(tr != NULL) {
        if(!tr->active) {
            tr = tr->next;
            data[i] = 0.;
            i++;
            continue;
        }
        x_to_ijk(tr->x, index);
        data[i] = NDP_ELEM_LINEAR(sim_p,index[0],index[1],index[2],RHO);
        tr = tr->next;
        i++;
    }
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/rho", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    i=0;
    tr = tracers;
    while(tr != NULL) {
        if(!tr->active) {
            tr = tr->next;
            data[i] = 0.;
            i++;
            continue;
        }
        x_to_ijk(tr->x, index);
        data[i] = NDP_ELEM_LINEAR(sim_p,index[0],index[1],index[2],YE);
        tr = tr->next;
        i++;
    }
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/ye", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    i=0;
    tr = tracers;
    while(tr != NULL) {
        if(!tr->active) {
            tr = tr->next;
            data[i] = 0.;
            i++;
            continue;
        }
        x_to_ijk(tr->x, index);
        data[i] = NDP_ELEM_LINEAR(sim_eos,index[0],index[1],index[2],TEMP);
        tr = tr->next;
        i++;
    }
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/T_MeV", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    i=0;
    tr = tracers;
    while(tr != NULL) {
        if(!tr->active) {
            tr = tr->next;
            data[i] = 0.;
            i++;
            continue;
        }
        x_to_ijk(tr->x, index);
        data[i] = NDP_ELEM_LINEAR(sim_eos,index[0],index[1],index[2],ENT);
        tr = tr->next;
        i++;
    }
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, "/Step#0/entropy", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
    else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    H5Sclose(filespace);
    H5Sclose(memspace);

    H5Fclose(file_id);

    if(data != NULL) free(data);
    if(f_coord != NULL) free(f_coord);
    if(pid != NULL) free(pid);

    pdump_id++;

}

/*#################################
  #################################

    ?-dim XDMF output functions

  #################################
  #################################*/


#if (NDIM==1)
FILE *write_xml_head1() {

    char name[80];
    FILE *fp;

    mkdir("dumps", 0777);

    sprintf(name,"dumps/dump_%05d.xmf", dump_id);
    fp = fopen(name,"w");
    if (fp == NULL) {
        fprintf(stderr,"Failed to open xmf file...perhaps there is no \"dumps\" directory???\n");
        exit(123);
    }

    fprintf(fp, "<?xml version=\"1.0\" ?>\n");
    fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(fp, "<Xdmf Version=\"2.0\">\n");
    fprintf(fp, " <Domain>\n");
    fprintf(fp, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
    fprintf(fp, "     <Time Value=\"%16.14e\"/>\n", t);
    fprintf(fp, "     <Topology TopologyType=\"1DMesh\" NumberOfElements=\"%d\"/>\n", n1+1);
    fprintf(fp, "     <Geometry Units=\"km\" GeometryType=\"X\">\n");
    fprintf(fp, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", n1+1, n2+1);
    fprintf(fp, "        grid.h5:/X\n");
    fprintf(fp, "       </DataItem>\n");
    fprintf(fp, "     </Geometry>\n");

    return fp;
}

void dump1()
{
  static int firstc = 1;
  static Real *data;
  static Real *data1,*data2,*data3,*data4,*data5,*data6;
  static double *kappa,*jg,*sc,*delta;
  double heat[3],cool[3];
  int i,j,v,dd,g,ind;
  int spec,ig;
  char name[80];
  FILE *xml=NULL, *fp;
  double vec[SPACEDIM];
  hid_t file_id,filespace,memspace,plist_id,group_id,dset_id;
  hid_t vec_filespace, vec_memspace;
  hsize_t mem_start[NDIM+1], file_start[NDIM+1], one, zero;
  hsize_t file_count[NDIM+1], file_grid_dims[NDIM+1];
  herr_t status;

  if (!dump_hdf) {
    if (mpi_io_proc()) mkdir("dumps", 0777);
    dump_ascii1(0);
    return;
  }

  TIMER_START("dump1");

  one = 1;
  zero = 0;

  //GPU_PRAGMA(omp target update to(sim_p[:cell_count_all][:nvars]))
  //update_eos(eostable,sim_p,sim_eos);
  //GPU_PRAGMA(omp target update from(sim_eos[:cell_count_all][:NEOS]))

  if (mpi_io_proc()) {
    xml = write_xml_head1();
    sprintf(name,"dumps/dumps.visit");
    if (dump_id == 0) fp = fopen(name,"w");
    else fp = fopen(name,"a");
    fprintf(fp,"dump_%05d.xmf\n", dump_id);
    fclose(fp);
  }

  if (firstc) {
    //if (mpi_io_proc()) {
    //  write_xml_expressions();
    //}
    dump_grid();
    data  = malloc_rank1(my_grid_dims[0], sizeof *data    );

    #if DO_RADIATION
    // AS:  Allocate additional data structures for dump heating/cooling rates
    kappa = malloc_rank1(ngroups,              sizeof *kappa);
    jg    = malloc_rank1(ngroups,              sizeof *jg   );
    sc    = malloc_rank1(ngroups,              sizeof *sc   );
    delta = malloc_rank1(ngroups,              sizeof *delta);
    data1 = malloc_rank1(NDIM*my_grid_dims[0], sizeof *data1);
    data2 = malloc_rank1(NDIM*my_grid_dims[0], sizeof *data2);
    #if NEUTRINO
    data3 = malloc_rank1(NDIM*my_grid_dims[0], sizeof *data3);
    data4 = malloc_rank1(NDIM*my_grid_dims[0], sizeof *data4);
    data5 = malloc_rank1(NDIM*my_grid_dims[0], sizeof *data5);
    data6 = malloc_rank1(NDIM*my_grid_dims[0], sizeof *data6);
    #endif /* NEUTRINO */
    #endif /* DO_RADIATION */

    firstc = 0;
  }

  sprintf(name,"dumps/dump_%05d.h5", dump_id);

  #if (USE_MPI==TRUE)
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  #else
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  #endif

  file_grid_dims[0] = 1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = (mpi_io_proc()) ? one : zero;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  DLOOP { file_grid_dims[dd] = global_grid_dims[dd]; }
  file_grid_dims[NDIM] = NDIM;    // for vectors
  filespace = H5Screate_simple(NDIM, file_grid_dims, NULL);
  DLOOP {
    file_start[dd] = istart[dd];
    file_count[dd] = my_grid_dims[dd];
  }
  file_start[NDIM] = 0;       // for vectors
  file_count[NDIM] = NDIM;    // "    "
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  vec_filespace = H5Screate_simple(NDIM+1, file_grid_dims, NULL);
  H5Sselect_hyperslab(vec_filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(NDIM, file_count, NULL);
  DLOOP { mem_start[dd] = 0; }
  mem_start[NDIM] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);
  vec_memspace = H5Screate_simple(NDIM+1, file_count, NULL);
  H5Sselect_hyperslab(vec_memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

#if 1
  // Extra output stuff for debugging...

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = ND_ELEM_LINEAR(sim_shock_flag,i,jj,kk);
  }
  write_scalar(data, file_id, "shock_flag", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = myrank;
  }
  write_scalar(data, file_id, "myrank", filespace, memspace, xml);
#endif

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
  }
  write_scalar(data, file_id, "rho", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,UU);
  }
  write_scalar(data, file_id, "u", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U1)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
  }
  write_scalar(data, file_id, "u1", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U2)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
  }
  write_scalar(data, file_id, "u2", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U3)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
  }
  write_scalar(data, file_id, "u3", filespace, memspace, xml);

  for (v=0; v<ncomp; v++) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U3+1+v);
    }
    sprintf(name, "comp%d", v);
    write_scalar(data, file_id, name, filespace, memspace, xml);
  }

  for (v=0; v<NEOS; v++) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind++] = NDP_ELEM_LINEAR(sim_eos,i,jj,kk,v);
    }
    sprintf(name, "eos%d", v);
    write_scalar(data, file_id, name, filespace, memspace, xml);
  }

  #if (GR_MONOPOLE==TRUE)
  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = ND_ELEM_LINEAR(geom,i,jj,kk).lapse[0];
  }
  write_scalar(data, file_id, "lapse", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) {
    data[ind++] = ND_ELEM_LINEAR(geom,i,jj,kk).lapse[1];
  }
  write_scalar(data, file_id, "lapse_edge", filespace, memspace, xml);
  #endif

  #if (DO_RADIATION)
  if (dump_rad_vars) {
    group_id = H5Gcreate(file_id, "/Erad0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Erad1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Erad2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, "/Frad0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, "/Frad0/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    double efac;
    if (strcmp(freq_type,"mev")==0) {
      efac = MEV_TO_ERG;
    }
    if (strcmp(freq_type,"hz")==0) {
      efac = HPLANCK;
    }

    GLOOP {

      if (g==0) {
        spec = 0;
        ig = g;
      }
      if (g==nr1) {
        spec = 1;
        ig = g - nr1;
      }
      if (g==nr1+nr2) {
        spec = 2;
        ig = g - nr1 - nr2;
      }

      ind = 0;
      for (i=istart[0]; i<istop[0]; i++) {
        data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)*efac/degroup[g]/RADIATION_NORM;
      }
      sprintf(name, "Erad%d/g%02d", spec, ig);
      write_scalar(data, file_id, name, filespace, memspace, xml);

      DLOOP {
        ind = 0;
        for (i=istart[0]; i<istop[0]; i++) {
          data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+NDIM*g+dd)/RADIATION_NORM*efac/degroup[g] * ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][dd];
        }
        sprintf(name, "Frad%d/dir%d/g%02d", spec, dd, ig);
        write_scalar(data, file_id, name, filespace, memspace, xml);
      }
      ig++;
    }

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind] = 0.0;
      for (g=0; g<nr1; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
      }
      ind++;
    }
    write_scalar(data, file_id, "Erad0/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind] = 0.0;
      for (g=0; g<nr1; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad0/dir0/total", filespace, memspace, xml);

    #if (NEUTRINO)
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
      }
      ind++;
    }
    write_scalar(data, file_id, "Erad1/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad1/dir0/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind] = 0.0;
      for (g=nr1+nr2; g<nr1+nr2+nr3; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
      }
      ind++;
    }
    write_scalar(data, file_id, "Erad2/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) {
      data[ind] = 0.0;
      for (g=nr1+nr2; g<nr1+nr2+nr3; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad2/dir0/total", filespace, memspace, xml);
    #endif

    // AS:  Dump local heating/cooling rates
    ind = 0;
    for (i=istart[0];i<istop[0];i++) {
      opac_emis(NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO), NDP_ELEM_LINEAR(sim_eos,i,jj,kk,TEMP), NDP_ELEM_LINEAR(sim_p,i,jj,kk,YE), kappa, jg, sc, delta, NULL, NULL, NULL);
      heat[0] = 0.0;
      cool[0] = 0.0;
      for (int g=0;g<nr1;g++) {
        heat[0] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[0] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }

      data1[ind] = heat[0];
      data2[ind] = cool[0];

      #if NEUTRINO
      heat[1] = 0.0;
      cool[1] = 0.0;
      for (int g=nr1;g<nr1+nr2;g++) {
        heat[1] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[1] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }
      heat[2] = 0.0;
      cool[2] = 0.0;
      for (int g=nr1+nr2;g<nr1+nr2+nr3;g++) {
        heat[2] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[2] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }

      data3[ind] = heat[1];
      data4[ind] = cool[1];
      data5[ind] = heat[2];
      data6[ind] = cool[2];
      #endif /* NEUTRINO */
      ind++;
    }
    write_scalar(data1, file_id, "Erad0/heat", filespace, memspace, xml);
    write_scalar(data2, file_id, "Erad0/cool", filespace, memspace, xml);
    #if NEUTRINO
    write_scalar(data3, file_id, "Erad1/heat", filespace, memspace, xml);
    write_scalar(data4, file_id, "Erad1/cool", filespace, memspace, xml);
    write_scalar(data5, file_id, "Erad2/heat", filespace, memspace, xml);
    write_scalar(data6, file_id, "Erad2/cool", filespace, memspace, xml);
    #endif /* NEUTRINO */
  } /* dump_rad_vars */
  #endif /* DO_RADIATION */

  if (mpi_io_proc()) {
    write_xml_closing(xml);
  }

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);

  dump_id++;
  TIMER_STOP;
  return;
}

#endif



#if (NDIM==2)
FILE *write_xml_head2()
{
  char name[80];
  FILE *fp;

  mkdir("dumps", 0777);

  sprintf(name,"dumps/dump_%05d.xmf", dump_id);
  fp = fopen(name,"w");
  if (fp == NULL) {
      fprintf(stderr,"Failed to open xmf file...perhaps there is no \"dumps\" directory???\n");
      exit(123);
  }

  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"2.0\">\n");
  fprintf(fp, " <Domain>\n");
  fprintf(fp, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
  fprintf(fp, "     <Time Value=\"%16.14e\"/>\n", t);
  fprintf(fp, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", n1+1, n2+1);
  fprintf(fp, "     <Geometry Units=\"km\" GeometryType=\"X_Y\">\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", n1+1, n2+1);
  fprintf(fp, "        grid.h5:/X\n");
  fprintf(fp, "       </DataItem>\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", n1+1, n2+1);
  fprintf(fp, "        grid.h5:/Y\n");
  fprintf(fp, "       </DataItem>\n");
  fprintf(fp, "     </Geometry>\n");

  return fp;
}

void write_xml_expressions()
{
  FILE *fp = NULL;

  fp = fopen("dumps/expressions.xml","w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");

  // convenience
  fprintf(fp, "<Object name=\"ExpressionList\">\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">rsq</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"coord(mesh)[0]^2+coord(mesh)[1]^2\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">r</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"sqrt(rsq)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">sth</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"coord(mesh)[0]/r\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">cth</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"coord(mesh)[1]/r\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  // spherical unit vectors
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">nr</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{sth, cth}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">nth</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{cth, -sth}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  // velocity vector
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">velocity</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{velocity0, velocity1}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  #if (DO_RADIATION)
  // flux vectors/luminosities/etc
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Flux0</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{Flux0_0, Flux0_1}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Fr0</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"dot(Flux0,nr)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Fth0</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"dot(Flux0,nth)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">L0</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"12.56637061436*rsq*Fr0\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">reduced_flux0</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"magnitude(Flux0)/(2.99792458e10*\\<Erad0/total\\>)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  #if (NEUTRINO)
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Flux1</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{Flux1_0, Flux1_1}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Fr1</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"dot(Flux1,nr)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Fth1</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"dot(Flux1,nth)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">L1</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"12.56637061436*rsq*Fr1\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">reduced_flux1</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"magnitude(Flux1)/(2.99792458e10*\\<Erad1/total\\>)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Flux2</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{Flux2_0, Flux2_1}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Fr2</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"dot(Flux2,nr)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Fth2</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"dot(Flux2,nth)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">L2</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"12.56637061436*rsq*Fr2\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">reduced_flux2</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"magnitude(Flux2)/(2.99792458e10*\\<Erad2/total\\>)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  #endif /* NEUTRINO */
  #endif /* DO_RADIATION */
  fprintf(fp, "</Object>\n");

  fclose(fp);

  return;
}

void dump2()
{
  static int firstc = 1;
  static Real *data;
  static Real *data1,*data2,*data3,*data4,*data5,*data6;
  static double *kappa,*jg,*sc,*delta;
  double heat[3],cool[3];
  double efac;
  int i,j,v,dd,jj,g;
  int ind,spec,ig,total_cells;
  char name[80];
  FILE *xml=NULL, *fp;
  double vec[SPACEDIM];
  hid_t file_id,filespace,memspace,plist_id,group_id,dset_id;
  hid_t vec_filespace, vec_memspace;
  hsize_t mem_start[NDIM+1], file_start[NDIM+1], one, zero;
  hsize_t file_count[NDIM+1], file_grid_dims[NDIM+1];
  herr_t status;

  if (!dump_hdf) return;

  TIMER_START("dump2");

  one = 1;
  zero = 0;

  //GPU_PRAGMA(omp target update to(sim_p[:cell_count_all][:nvars]))
  //update_eos(eostable,sim_p,sim_eos);
  //GPU_PRAGMA(omp target update from(sim_eos[:cell_count_all][:NEOS]))

  if (mpi_io_proc()) {
    xml = write_xml_head2();
    sprintf(name,"dumps/dumps.visit");
    if (dump_id == 0) fp = fopen(name,"w");
    else fp = fopen(name,"a");
    fprintf(fp,"dump_%05d.xmf\n", dump_id);
    fclose(fp);
  }

  if (firstc) {
    if (mpi_io_proc()) {
      write_xml_expressions();
    }
    dump_grid();
    data      = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data     );
    firstc = 0;

    #if DO_RADIATION
    // AS:  Allocate additional data structures for dump heating/cooling rates
    kappa = malloc_rank1(ngroups,                              sizeof *kappa);
    jg    = malloc_rank1(ngroups,                              sizeof *jg   );
    sc    = malloc_rank1(ngroups,                              sizeof *sc   );
    delta = malloc_rank1(ngroups,                              sizeof *delta);
    data1 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data1);
    data2 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data2);
    #if NEUTRINO
    data3 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data3);
    data4 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data4);
    data5 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data5);
    data6 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1], sizeof *data6);
    #endif /* NEUTRINO */
    #endif /* DO_RADIATION */
  }

  sprintf(name,"dumps/dump_%05d.h5", dump_id);

  #if (USE_MPI==TRUE)
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  #else
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  #endif

  file_grid_dims[0] = 1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? one : zero;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  DLOOP { file_grid_dims[dd] = global_grid_dims[dd]; }
  file_grid_dims[NDIM] = NDIM;    // for vectors
  filespace = H5Screate_simple(NDIM, file_grid_dims, NULL);
  DLOOP {
    file_start[dd] = istart[dd];
    file_count[dd] = my_grid_dims[dd];
  }
  file_start[NDIM] = 0;       // for vectors
  file_count[NDIM] = NDIM;    // "    "
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  vec_filespace = H5Screate_simple(NDIM+1, file_grid_dims, NULL);
  H5Sselect_hyperslab(vec_filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(NDIM, file_count, NULL);
  DLOOP { mem_start[dd] = 0; }
  mem_start[NDIM] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);
  vec_memspace = H5Screate_simple(NDIM+1, file_count, NULL);
  H5Sselect_hyperslab(vec_memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

#if 1
  // Extra output stuff for debugging...

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    data[ind++] = ND_ELEM_LINEAR(sim_shock_flag,i,jj,kk);
  }
  write_scalar(data, file_id, "shock_flag", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    data[ind++] = myrank;
  }
  write_scalar(data, file_id, "myrank", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    data[ind++] = DJS(i);
  }
  write_scalar(data, file_id, "dj", filespace, memspace, xml);

  // ind = 0;
  // for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
  //   jj = JS(i,j);
  //   data[ind++] = ND_ELEM(sim_Phi,i,jj,kk);
  // }
  // write_scalar(data, file_id, "Phi", filespace, memspace, xml);
#endif

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
  }
  write_scalar(data, file_id, "rho", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,UU);
  }
  write_scalar(data, file_id, "u", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U1)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
  }
  write_scalar(data, file_id, "u1", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U2)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
  }
  write_scalar(data, file_id, "u2", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U3)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
  }
  write_scalar(data, file_id, "u3", filespace, memspace, xml);

  for (v=0; v<ncomp; v++) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U3+1+v);
    }
    sprintf(name, "comp%d", v);
    write_scalar(data, file_id, name, filespace, memspace, xml);
  }

  for (v=0; v<NEOS; v++) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind++] = NDP_ELEM_LINEAR(sim_eos,i,jj,kk,v);
    }
    sprintf(name, "eos%d", v);
    write_scalar(data, file_id, name, filespace, memspace, xml);
  }

  #if (GR_MONOPOLE==TRUE)
  ind = 0;
  for (i=istart[0];i<istop[0];i++) for (j=istart[1];j<istop[1];j++) {
    jj = JS(i,j);
    data[ind] = ND_ELEM_LINEAR(geom,i,jj,kk).lapse[0];
    ind++;
  }
  write_scalar(data, file_id, "lapse", filespace, memspace, xml);
  ind = 0;
  for (i=istart[0];i<istop[0];i++) for (j=istart[1];j<istop[1];j++) {
    jj = JS(i,j);
    data[ind] = ND_ELEM_LINEAR(geom,i,jj,kk).lapse[1];
    ind++;
  }
  write_scalar(data, file_id, "lapse_edge", filespace, memspace, xml);
  #endif

  #if (DO_RADIATION)
  if (dump_rad_vars) {
    group_id = H5Gcreate(file_id, "/Erad0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Erad1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Erad2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, "/Frad0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    group_id = H5Gcreate(file_id, "/Frad0/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad0/dir1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1/dir1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2/dir1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    if (strcmp(freq_type,"mev")==0) {
      efac = MEV_TO_ERG;
    }
    if (strcmp(freq_type,"hz")==0) {
      efac = HPLANCK;
    }

    GLOOP {
      if (g==0) {
        spec = 0;
        ig = g;
      }
      if (g==nr1) {
        spec = 1;
        ig = g - nr1;
      }
      if (g==nr1+nr2) {
        spec = 2;
        ig = g - nr1 - nr2;
      }

      ind = 0;
      for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
        jj = JS(i,j);
        data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM*efac/degroup[g];
      }
      sprintf(name, "Erad%d/g%02d", spec, ig);
      write_scalar(data, file_id, name, filespace, memspace, xml);

      DLOOP {
        ind = 0;
        for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
          jj = JS(i,j);
          data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+NDIM*g+dd)/RADIATION_NORM*efac/degroup[g] * ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][dd];
        }
        sprintf(name, "Frad%d/dir%d/g%02d", spec, dd, ig);
        write_scalar(data, file_id, name, filespace, memspace, xml);
      }
      ig++;
    }

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=0; g<nr1; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
      }
      ind++;
    }
    write_scalar(data, file_id, "Erad0/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=0; g<nr1; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g*NDIM)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad0/dir0/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=0; g<nr1; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g*NDIM+1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad0/dir1/total", filespace, memspace, xml);

    #if (NEUTRINO==TRUE)
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
      }
      ind++;
    }
    write_scalar(data, file_id, "Erad1/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g*NDIM)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad1/dir0/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g*NDIM+1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad1/dir1/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=nr1+nr2; g<ngroups; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
      }
      ind++;
    }
    write_scalar(data, file_id, "Erad2/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=nr1+nr2; g<ngroups; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g*NDIM)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad2/dir0/total", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      data[ind] = 0.0;
      for (g=nr1+nr2; g<ngroups; g++) {
        data[ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+g*NDIM+1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
      }
      ind++;
    }
    write_scalar(data, file_id, "Frad2/dir1/total", filespace, memspace, xml);
    #endif /* NEUTRINO */
  }  /* if (dump_rad_vars) */
  #endif /* DO_RADIATION */

  // now write the vectors
  ind = 0;
  total_cells = my_grid_dims[0]*my_grid_dims[1];
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
    jj = JS(i,j);
    memcpy(vec, &NDP_ELEM_LINEAR(sim_p,i,jj,kk,U1), 2*sizeof(double));
    vec[0] *= ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
    vec[1] *= ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
    vec_transform_to_xyz(vec,i,j,0);
    data[ind + 0*total_cells] = vec[0];
    data[ind + 1*total_cells] = vec[1];
    ind++;
    //data[ind] = vec[2]; ind++;
  }
  write_vector(data, file_id, "velocity", filespace, memspace, xml);

  #if (DO_RADIATION)
  if (dump_rad_vars) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = 0.0;
      for (g=0; g<nr1; g++) {
        vec[0] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 0)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
        vec[1] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
      }
      vec_transform_to_xyz(vec,i,j,0);
      data[ind] = vec[0];
      data[ind+total_cells] = vec[1];
      ind++;
    }
    write_vector(data, file_id, "Flux0_", filespace, memspace, xml);

    #if (NEUTRINO)
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        vec[0] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 0)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
        vec[1] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
      }
      vec_transform_to_xyz(vec,i,j,0);
      data[ind] = vec[0];
      data[ind+total_cells] = vec[1];
      ind++;
    }
    write_vector(data, file_id, "Flux1_", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = 0.0;
      for (g=nr1+nr2; g<nr1+nr2+nr3; g++) {
        vec[0] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 0)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
        vec[1] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
      }
      vec_transform_to_xyz(vec,i,j,0);
      data[ind] = vec[0];
      data[ind+total_cells] = vec[1];
      ind++;
    }
    write_vector(data, file_id, "Flux2_", filespace, memspace, xml);
    #endif /* NEUTRINO */

    // AS:  Dump local heating/cooling rates
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) {
      jj = JS(i,j);
      opac_emis(NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO), NDP_ELEM_LINEAR(sim_eos,i,jj,kk,TEMP), NDP_ELEM_LINEAR(sim_p,i,jj,kk,YE), kappa, jg, sc, delta, NULL, NULL, NULL);
      heat[0] = 0.0;
      cool[0] = 0.0;
      for (int g=0; g<nr1; g++) {
        heat[0] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[0] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }
      heat[1] = 0.0;
      cool[1] = 0.0;
      for (int g=nr1; g<nr1+nr2; g++) {
        heat[1] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[1] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }
      heat[2] = 0.0;
      cool[2] = 0.0;
      for (int g=nr1+nr2; g<nr1+nr2+nr3; g++) {
        heat[2] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[2] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }

      data1[ind] = heat[0];
      data2[ind] = cool[0];
      data3[ind] = heat[1];
      data4[ind] = cool[1];
      data5[ind] = heat[2];
      data6[ind] = cool[2];
      ind++;
    }
    write_scalar(data1, file_id, "Erad0/heat", filespace, memspace, xml);
    write_scalar(data2, file_id, "Erad0/cool", filespace, memspace, xml);
    write_scalar(data3, file_id, "Erad1/heat", filespace, memspace, xml);
    write_scalar(data4, file_id, "Erad1/cool", filespace, memspace, xml);
    write_scalar(data5, file_id, "Erad2/heat", filespace, memspace, xml);
    write_scalar(data6, file_id, "Erad2/cool", filespace, memspace, xml);
  }  /* if (dump_rad_vars) */
  #endif /* DO_RADIATION */

  if (mpi_io_proc()) {
    write_xml_closing(xml);
  }

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);

  dump_id++;
  TIMER_STOP;
  return;
}


#endif /* NDIM==2 */


#if (NDIM==3)
FILE *write_xml_head3()
{
  char name[80];
  FILE *fp;

  mkdir("dumps", 0777);

  sprintf(name,"dumps/dump_%05d.xmf", dump_id);
  fp = fopen(name,"w");
  if (fp == NULL) {
    fprintf(stderr,"Failed to open xmf file...perhaps there is no \"dumps\" directory???\n");
    exit(123);
  }

  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"2.0\">\n");
  fprintf(fp, " <Domain>\n");
  fprintf(fp, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
  fprintf(fp, "     <Time Value=\"%16.14e\"/>\n", t);
  fprintf(fp, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", n1+1, n2+1, n3+1);
  fprintf(fp, "     <Geometry GeometryType=\"X_Y_Z\">\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", n1+1, n2+1, n3+1);
  fprintf(fp, "        grid.h5:/X\n");
  fprintf(fp, "       </DataItem>\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", n1+1, n2+1, n3+1);
  fprintf(fp, "        grid.h5:/Y\n");
  fprintf(fp, "       </DataItem>\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", n1+1, n2+1, n3+1);
  fprintf(fp, "        grid.h5:/Z\n");
  fprintf(fp, "       </DataItem>\n");

  fprintf(fp, "     </Geometry>\n");

  fflush(fp);

  return fp;
}

void write_xml_expressions()
{
  FILE *fp = NULL;

  fp = fopen("dumps/expressions.xml","w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");

  // convenience
  fprintf(fp, "<Object name=\"ExpressionList\">\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">rsq</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"coord(mesh)[0]^2+coord(mesh)[1]^2+coord(mesh)[2]^2\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">r</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"sqrt(rsq)\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">ScalarMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  // velocity vector
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">velocity</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{velocity0, velocity1, velocity2}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");

  #if (DO_RADIATION)
  // flux vectors
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Flux0</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{Flux0_0, Flux0_1, Flux0_2}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  #if (NEUTRINO)
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Flux1</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{Flux1_0, Flux1_1, Flux1_2}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  fprintf(fp, "  <Object name=\"Expression\">\n");
  fprintf(fp, "    <Field name=\"name\" type=\"string\">Flux2</Field>\n");
  fprintf(fp, "    <Field name=\"definition\" type=\"string\">\"{Flux2_0, Flux2_1, Flux2_2}\"</Field>\n");
  fprintf(fp, "    <Field name=\"type\" type=\"string\">VectorMeshVar</Field>\n");
  fprintf(fp, "  </Object>\n");
  #endif
  #endif
  fprintf(fp, "</Object>\n");

  fclose(fp);

  return;
}

void dump3()
{
  static int firstc = 1;
  static Real *data;
  static Real *data1,*data2,*data3,*data4,*data5,*data6;
  static double *kappa,*jg,*sc,*delta;
  double heat[3],cool[3];
  double efac;
  int i,j,k,v,dd,jj,kk,g;
  int ind,spec,ig,total_cells;
  char name[80];
  FILE *xml=NULL, *fp;
  double vec[SPACEDIM];
  hid_t file_id,filespace,memspace,plist_id,group_id,dset_id;
  hid_t vec_filespace, vec_memspace;
  hsize_t mem_start[NDIM+1], file_start[NDIM+1], one, zero;
  hsize_t file_count[NDIM+1], file_grid_dims[NDIM+1];
  herr_t status;

  if (!dump_hdf) return;

  TIMER_START("dump3");

  one = 1;
  zero = 0;

  //GPU_PRAGMA(omp target update to(sim_p[:cell_count_all][:nvars]))
  //update_eos(eostable,sim_p,sim_eos);
  //GPU_PRAGMA(omp target update from(sim_eos[:cell_count_all][:NEOS]))

  if (mpi_io_proc()) {
    xml = write_xml_head3();
    sprintf(name,"dumps/dumps.visit");
    if (dump_id == 0) fp = fopen(name,"w");
    else fp = fopen(name,"a");
    fprintf(fp,"dump_%05d.xmf\n", dump_id);
    fclose(fp);
  }

  if (firstc) {
    if (mpi_io_proc()) {
      write_xml_expressions();
    }
    dump_grid();
    data = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data);
    firstc = 0;

    #if DO_RADIATION
    // AS:  Allocate additional data structures for dump heating/cooling rates
    kappa = malloc_rank1(ngroups,                                              sizeof *kappa);
    jg    = malloc_rank1(ngroups,                                              sizeof *jg   );
    sc    = malloc_rank1(ngroups,                                              sizeof *sc   );
    delta = malloc_rank1(ngroups,                                              sizeof *delta);
    data1 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data1);
    data2 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data2);
    #if NEUTRINO
    data3 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data3);
    data4 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data4);
    data5 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data5);
    data6 = malloc_rank1(NDIM*my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2], sizeof *data6);
    #endif /* NEUTRINO */
    #endif /* DO_RADIATION */
  }

  sprintf(name,"dumps/dump_%05d.h5", dump_id);

  #if (USE_MPI==TRUE)
  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  if (!file_id && mpi_io_proc()) {
    fprintf(stderr,"Failed to create dump file %s\n", name);
    exit(2);
  }
  H5Pclose(plist_id);
  #else
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  #endif

  file_grid_dims[0] = 1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = (mpi_io_proc()) ? one : zero;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace  = H5Screate_simple(1, (mpi_io_proc()) ? &one : &zero, NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  DLOOP { file_grid_dims[dd] = global_grid_dims[dd]; }
  file_grid_dims[NDIM] = NDIM;    // for vectors
  filespace = H5Screate_simple(NDIM, file_grid_dims, NULL);
  DLOOP {
    file_start[dd] = istart[dd];
    file_count[dd] = my_grid_dims[dd];
  }
  file_start[NDIM] = 0;       // for vectors
  file_count[NDIM] = NDIM;    // "    "
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(NDIM, file_count, NULL);
  DLOOP { mem_start[dd] = 0; }
  mem_start[NDIM] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

#if 0
  // Extra output stuff for debugging...

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = ND_ELEM_LINEAR(sim_shock_flag,i,jj,kk);
  }
  write_scalar(data, file_id, "shock_flag", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    data[ind++] = myrank;
  }
  write_scalar(data, file_id, "myrank", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    data[ind++] = DJS(i);
  }
  write_scalar(data, file_id, "dj", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    data[ind++] = DKS(i,jj);
  }
  write_scalar(data, file_id, "dk", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = ND_ELEM(sim_Phi,i,jj,kk);
  }
  write_scalar(data, file_id, "Phi", filespace, memspace, xml);
#endif


  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
  }
  write_scalar(data, file_id, "rho", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,UU);
  }
  write_scalar(data, file_id, "u", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U1)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
  }
  write_scalar(data, file_id, "u1", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U2)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
  }
  write_scalar(data, file_id, "u2", filespace, memspace, xml);

  ind=0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U3)*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
  }
  write_scalar(data, file_id, "u3", filespace, memspace, xml);

  for (v=0; v<ncomp; v++) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
      jj = JS(i,j);
      kk = KS(i,jj,k);
      data[ind++] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,U3+1+v);
    }
    sprintf(name, "comp%d", v);
    write_scalar(data, file_id, name, filespace, memspace, xml);
  }

  for (v=0; v<NEOS; v++) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
      jj = JS(i,j);
      kk = KS(i,jj,k);
      data[ind++] = NDP_ELEM_LINEAR(sim_eos,i,jj,kk,v);
    }
    sprintf(name, "eos%d", v);
    write_scalar(data, file_id, name, filespace, memspace, xml);
  }

  #if (GR_MONOPOLE==TRUE)
  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind] = ND_ELEM_LINEAR(geom,i,jj,kk).lapse[0];
    ind++;
  }
  write_scalar(data, file_id, "lapse", filespace, memspace, xml);

  ind = 0;
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    data[ind] = ND_ELEM_LINEAR(geom,i,jj,kk).lapse[1];
    ind++;
  }
  write_scalar(data, file_id, "lapse_edge", filespace, memspace, xml);
  #endif

  #if (DO_RADIATION)
  if (dump_rad_vars) {
    group_id = H5Gcreate(file_id, "/Erad0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    #if (NEUTRINO)
    group_id = H5Gcreate(file_id, "/Erad1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Erad2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    #endif

    group_id = H5Gcreate(file_id, "/Frad0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    #if (NEUTRINO)
    group_id = H5Gcreate(file_id, "/Frad1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    #endif

    group_id = H5Gcreate(file_id, "/Frad0/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad0/dir1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad0/dir2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    #if (NEUTRINO)
    group_id = H5Gcreate(file_id, "/Frad1/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1/dir1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad1/dir2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2/dir0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2/dir1", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    group_id = H5Gcreate(file_id, "/Frad2/dir2", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    #endif

    double efac;
    if (strcmp(freq_type,"mev")==0) {
      efac = MEV_TO_ERG;
    }
    if (strcmp(freq_type,"hz")==0) {
      efac = HPLANCK;
    }

    total_cells = my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2];
    memset(data1, 0, NDIM*total_cells*sizeof(Real));
    memset(data2, 0, NDIM*total_cells*sizeof(Real));
    #if (NEUTRINO==TRUE)
    memset(data3, 0, NDIM*total_cells*sizeof(Real));
    memset(data4, 0, NDIM*total_cells*sizeof(Real));
    memset(data5, 0, NDIM*total_cells*sizeof(Real));
    memset(data6, 0, NDIM*total_cells*sizeof(Real));
    #endif
    Real * data_Etot[3] = {data1, data3, data5};
    Real * data_Ftot[3][3] = {
      {&data2[0*total_cells], &data2[1*total_cells], &data2[2*total_cells]},
      {&data4[0*total_cells], &data4[1*total_cells], &data4[2*total_cells]},
      {&data6[0*total_cells], &data6[1*total_cells], &data6[2*total_cells]},
    };

    GLOOP {
      if (g==0) {
        spec = 0;
        ig = g;
      }
      if (g==nr1) {
        spec = 1;
        ig = g - nr1;
      }
      if (g==nr1+nr2) {
        spec = 2;
        ig = g - nr1 - nr2;
      }

      ind = 0;
      for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
        jj = JS(i,j);
        kk = KS(i,jj,k);
        data[ind] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM*efac/degroup[g];
        data_Etot[spec][ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/RADIATION_NORM;
        ind++;
      }
      sprintf(name, "Erad%d/g%02d", spec, ig);
      write_scalar(data, file_id, name, filespace, memspace, xml);

      DLOOP {
        ind = 0;
        for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
          jj = JS(i,j);
          kk = KS(i,jj,k);
          data[ind] = NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+NDIM*g+dd)/RADIATION_NORM*efac/degroup[g]*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][dd];
          data_Ftot[spec][dd][ind] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1+NDIM*g+dd)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][dd];
          ind++;
        }
        sprintf(name, "Frad%d/dir%d/g%02d", spec, dd, ig);
        write_scalar(data, file_id, name, filespace, memspace, xml);
      }
      ig++;
    }
    #if (NEUTRINO==TRUE)
    int const nspec = 3;
    #else
    int const nspec = 1;
    #endif
    for(int spec = 0; spec < nspec; ++spec) {
      sprintf(name, "Erad%d/total", spec);
      write_scalar(data_Etot[spec], file_id, name, filespace, memspace, xml);
      DLOOP {
        sprintf(name, "Frad%d/dir%d/total", spec, dd);
        write_scalar(data_Ftot[spec][dd], file_id, name, filespace, memspace, xml);
      }
    }
  }  /* if (dump_rad_vars) */
  #endif /* DO_RADIATION */

  // now write the vectors
  ind = 0;
  total_cells = my_grid_dims[0]*my_grid_dims[1]*my_grid_dims[2];
  for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
    jj = JS(i,j);
    kk = KS(i,jj,k);
    memcpy(vec, &NDP_ELEM_LINEAR(sim_p,i,jj,kk,U1), 3*sizeof(double));
    vec[0] *= ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
    vec[1] *= ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
    vec[2] *= ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
    vec_transform_to_xyz(vec,i,j,k);
    data[ind + 0*total_cells] = vec[0];
    data[ind + 1*total_cells] = vec[1];
    data[ind + 2*total_cells] = vec[2];
    ind++;
  }
  write_vector(data, file_id, "velocity", filespace, memspace, xml);

  #if (DO_RADIATION)
  if (dump_rad_vars) {
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
      jj = JS(i,j);
      kk = KS(i,jj,k);
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = 0.0;
      for (g=0; g<nr1; g++) {
        vec[0] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 0)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
        vec[1] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
        vec[2] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 2)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
      }
      vec_transform_to_xyz(vec,i,j,k);
      data[ind + 0*total_cells] = vec[0];
      data[ind + 1*total_cells] = vec[1];
      data[ind + 2*total_cells] = vec[2];
      ind++;
    }
    write_vector(data, file_id, "Flux0_", filespace, memspace, xml);

    #if (NEUTRINO)
    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
      jj = JS(i,j);
      kk = KS(i,jj,k);
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = 0.0;
      for (g=nr1; g<nr1+nr2; g++) {
        vec[0] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM    )/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
        vec[1] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
        vec[2] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 2)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
      }
      vec_transform_to_xyz(vec,i,j,k);
      data[ind + 0*total_cells] = vec[0];
      data[ind + 1*total_cells] = vec[1];
      data[ind + 2*total_cells] = vec[2];
      ind++;
    }
    write_vector(data, file_id, "Flux1_", filespace, memspace, xml);

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
      jj = JS(i,j);
      kk = KS(i,jj,k);
      vec[0] = 0.0;
      vec[1] = 0.0;
      vec[2] = 0.0;
      for (g=nr1+nr2; g<nr1+nr2+nr3; g++) {
        vec[0] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM    )/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][0];
        vec[1] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 1)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][1];
        vec[2] += NDP_ELEM_LINEAR(sim_p,i,jj,kk,ifrad1 + g*NDIM + 2)/RADIATION_NORM*ND_ELEM_LINEAR(geom,i,jj,kk).scale[0][2];
      }
      vec_transform_to_xyz(vec,i,j,k);
      data[ind + 0*total_cells] = vec[0];
      data[ind + 1*total_cells] = vec[1];
      data[ind + 2*total_cells] = vec[2];
      ind++;
    }
    write_vector(data, file_id, "Flux2_", filespace, memspace, xml);
    #endif  /* NEUTRINO */

    ind = 0;
    for (i=istart[0]; i<istop[0]; i++) for (j=istart[1]; j<istop[1]; j++) for (k=istart[2]; k<istop[2]; k++) {
      jj = JS(i,j);
      kk = KS(i,jj,k);

      opac_emis(NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO), NDP_ELEM_LINEAR(sim_eos,i,jj,kk,TEMP), NDP_ELEM_LINEAR(sim_p,i,jj,kk,YE), kappa, jg, sc, delta, NULL, NULL, NULL);
      heat[0] = 0.0;
      cool[0] = 0.0;
      for (int g=0; g<nr1; g++) {
        heat[0] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[0] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }
      #if NEUTRINO
      heat[1] = 0.0;
      cool[1] = 0.0;
      for (int g=nr1; g<nr1+nr2; g++) {
        heat[1] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[1] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }
      heat[2] = 0.0;
      cool[2] = 0.0;
      for (int g=nr1+nr2; g<nr1+nr2+nr3; g++) {
        heat[2] += CLIGHT*kappa[g]*NDP_ELEM_LINEAR(sim_p,i,jj,kk,irad1+g)/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
        cool[2] += jg[g]/NDP_ELEM_LINEAR(sim_p,i,jj,kk,RHO);
      }
      #endif

      data1[ind] = heat[0];
      data2[ind] = cool[0];
      #if NEUTRINO
      data3[ind] = heat[1];
      data4[ind] = cool[1];
      data5[ind] = heat[2];
      data6[ind] = cool[2];
      #endif
      ind++;
    }
    write_scalar(data1, file_id, "Erad0/heat", filespace, memspace, xml);
    write_scalar(data2, file_id, "Erad0/cool", filespace, memspace, xml);
    #if NEUTRINO
    write_scalar(data3, file_id, "Erad1/heat", filespace, memspace, xml);
    write_scalar(data4, file_id, "Erad1/cool", filespace, memspace, xml);
    write_scalar(data5, file_id, "Erad2/heat", filespace, memspace, xml);
    write_scalar(data6, file_id, "Erad2/cool", filespace, memspace, xml);
    #endif
  }  /* if (dump_rad_vars) */
  #endif  /* RADIATION */

  if (mpi_io_proc()) {
    write_xml_closing(xml);
  }

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);

  dump_id++;
  TIMER_STOP;
  return;
}

void read3(char const * fname)
{
  int ii,jj,kk;
  hid_t file_id;

  TIMER_START("read3");

  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0) {
    fprintf(stdout, "Could not open file: \"%s\"\n", fname);
    fflush(stdout);
    abort();
  }
  read_scalar(file_id, "rho", RHO, sim_p);
  read_scalar(file_id, "comp0", YE, sim_p);
  read_scalar(file_id, "u1", U1, sim_p);
  read_scalar(file_id, "u2", U2, sim_p);
  read_scalar(file_id, "u3", U3, sim_p);
  read_scalar(file_id, "eos2", TEMP, sim_eos);
  H5Fclose(file_id);

  ZLOOP {
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,UU) = u_given_rtx(eostable,
        NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO),
        NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP),
        NDP_ELEM_LINEAR(sim_p,ii,jj,kk,YE));
  }

  reset_boundaries(sim_p,0);
  complete_mpi_communication(0);

  ZGLOOP {
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1) /= ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0];
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U2) /= ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][1];
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U3) /= ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][2];
  }

  TIMER_STOP;
}
#endif


void restart_dump()
{
  restart_dump_v2();
}

void restart_dump_v1()
{
  int i,j,k,dd,vv,p,jstart,jstop,kstart;
  hsize_t file_grid_dims[NDIM+1],file_start[NDIM+1],file_count[NDIM+1];
  hsize_t mem_grid_dims[NDIM+1],mem_start[NDIM+1],one,zero;
  herr_t status;
  char name[256];
  double *pstart=NULL;

  #if (NDIM>1)
  // If not dumping HDF, then disable restarts in 2-d and 3-d
  if (!dump_hdf) return;
  #endif

  TIMER_START("restart_dump");

  one = 1;
  zero = 0;
  if (mpi_io_proc())
    mkdir("restarts", 0777);

  sprintf(name,"restarts/restart_%07d.h5", istep);
  if (mpi_io_proc()) {
    fprintf(stderr,"writing restart file %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  #if (USE_MPI==TRUE)
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // first write some info about the current state of the run (e.g. time, istep, ...)
  file_grid_dims[0] = 1;
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? one : zero;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  hid_t memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, "NDIM", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  int ndim = NDIM;
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ndim);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Step", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &istep);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time_Next_Dump", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_dump);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time_Next_Pdump", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_pdump);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time_Next_Restart", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_restart);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Step_Next_Restart", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &i_next_restart);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Next_Dump_ID", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &dump_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Next_Pdump_ID", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &pdump_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  dset_id = H5Dcreate(file_id, "N_Tracer", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  int ntracer = n_tracer_global;
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ntracer);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(filespace);
  H5Sclose(memspace);

  #if (GR_MONOPOLE==TRUE)
  /* dump lapse */
  file_grid_dims[0] = n1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "lapse", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse[0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  /* dump lapse edge */
  file_grid_dims[0] = n1+1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1+1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1+1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "lapse_edge", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse_edge[0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  #endif /* GR_MONOPOLE==TRUE */

  #if (NDIM==1)
  p = 0;
  for (i=0; i<n1; i++) {
    file_grid_dims[0] = nvars;
    filespace = H5Screate_simple(1, file_grid_dims, NULL);
    file_start[0] = 0;
    file_count[0] = nvars;
    if (i >= istart[0] && i < istop[0]) {
      jstart = JS(i,istart[1]);
      kstart = KS(i,jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,i,jstart,kstart,0);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    } else {
      jstart = JS(istart[0],istart[1]);
      kstart = KS(istart[0],jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,istart[0],jstart,kstart,0);
      file_count[0] = 0;
      H5Sselect_none(filespace);
    }

    mem_grid_dims[0] = nvars;
    memspace = H5Screate_simple(NDIM, mem_grid_dims, NULL);
    mem_start[0] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

    sprintf(name, "p%d", p++);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, pstart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
  }
  #endif /* NDIM==1 */

  #if (NDIM==2)
  p = 0;
  for (i=0; i<n1; i++) {
    file_grid_dims[0] = global_grid_dims[1]/dj[i];
    file_grid_dims[1] = nvars;
    filespace = H5Screate_simple(2, file_grid_dims, NULL);
    if (i >= istart[0] && i < istop[0]) {
      jstart = JS(i,istart[1]);
      kstart = KS(i,jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,i,jstart,kstart,0);
      file_count[0] = my_grid_dims[1]/dj[i];
      file_count[1] = nvars;
      file_start[0] = jstart;
      file_start[1] = 0;
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    } else {
      jstart = JS(istart[0],istart[1]);
      kstart = KS(istart[0],jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,istart[0],jstart,kstart,0);
      file_count[0] = 0;
      file_count[1] = 0;
      H5Sselect_none(filespace);
    }

    mem_grid_dims[0] = my_grid_dims[1]/dj[i] + 2*NG;
    mem_grid_dims[1] = nvars;
    memspace = H5Screate_simple(NDIM, mem_grid_dims, NULL);
    mem_start[0] = 0;
    mem_start[1] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

    sprintf(name,"p%d", p++);
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, pstart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
  }
  #endif /* NDIM==2 */

  #if (NDIM==3)
  p = 0;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2/dj[i]; j++) {
      file_grid_dims[0] = global_grid_dims[2]/dk[i][j];
      file_grid_dims[1] = nvars;
      filespace = H5Screate_simple(2, file_grid_dims, NULL);
      jstart = JS(i,istart[1]);
      jstop = JS(i,istop[1]);
      if (i >= istart[0] && i < istop[0] && j >= jstart && j < jstop) {
        kstart = KS(i,j,istart[2]);
        pstart = &NDP_ELEM_LINEAR(sim_p,i,j,kstart,0);
        file_count[0] = my_grid_dims[2]/dk[i][j];
        file_count[1] = nvars;
        file_start[0] = kstart;
        file_start[1] = 0;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
      } else {
        jstart = JS(istart[0],istart[1]);
        kstart = KS(istart[0],jstart,istart[2]);
        pstart = &NDP_ELEM_LINEAR(sim_p,istart[0],jstart,kstart,0);
        file_count[0] = 0;
        file_count[1] = 0;
        H5Sselect_none(filespace);
      }

      mem_grid_dims[0] = my_grid_dims[2]/dk[i][j];
      mem_grid_dims[1] = nvars;
      memspace = H5Screate_simple(2, mem_grid_dims, NULL);
      mem_start[0] = 0;
      mem_start[1] = 0;
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

      sprintf(name,"p%d", p++);
      plist_id = H5Pcreate(H5P_DATASET_CREATE);
      dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
      H5Pclose(plist_id);

      plist_id = H5Pcreate(H5P_DATASET_XFER);
      #if (USE_MPI==TRUE)
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      #endif
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, pstart);
      H5Dclose(dset_id);
      H5Pclose(plist_id);

      H5Sclose(memspace);
      H5Sclose(filespace);
    }
  }
  #endif /* NDIM==3 */

  hid_t group_id = H5Gcreate(file_id, "/Tracers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
  file_grid_dims[0] = n_tracer_global;
  hsize_t ntrace = n_tracer_current;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = n_tracer_current;
  //H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, (const hsize_t *)&f_coord);
  //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace  = H5Screate_simple(1, &ntrace, NULL);

  hsize_t *f_coord = malloc_rank1(n_tracer_current, sizeof(hsize_t));
  double *data = malloc_rank1(n_tracer_current, sizeof(double));
  char *flag = malloc_rank1(n_tracer_current, sizeof(char));

  i=0;
  tracer_t *tr = tracers;
  while(tr != NULL) {
      f_coord[i] = tr->id;
      data[i] = tr->x[0];
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/x", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  #if(NDIM>1)
  i=0;
  tr = tracers;
  while(tr != NULL) {
      data[i] = tr->x[1];
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/y", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  #endif
  #if(NDIM==3)
  i=0;
  tr = tracers;
  while(tr != NULL) {
      data[i] = tr->x[2];
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/z", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  #endif

  i=0;
  tr = tracers;
  while(tr != NULL) {
      data[i] = tr->mass;
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/mass", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  i=0;
  tr = tracers;
  while(tr != NULL) {
      flag[i] = tr->active;
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/active", H5T_NATIVE_CHAR, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_CHAR, memspace, filespace, plist_id, flag);;
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(f_coord!=NULL) free(f_coord);
  if(data!=NULL) free(data);
  if(flag!=NULL) free(flag);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  TIMER_STOP;

  return;
}

void restart_dump_v2()
{
  int i,j,k,dd,vv,p,jstart,jstop,kstart;
  hsize_t file_grid_dims[NDIM+1],file_start[NDIM+1],file_count[NDIM+1];
  hsize_t mem_grid_dims[NDIM+1],mem_start[NDIM+1],one,zero;
  herr_t status;
  char name[256];
  double *pstart=NULL;

  #if (NDIM>1)
  // If not dumping HDF, then disable restarts in 2-d and 3-d
  if (!dump_hdf) return;
  #endif

  TIMER_START("restart_dump");

  one = 1;
  zero = 0;
  if (mpi_io_proc())
    mkdir("restarts", 0777);

  sprintf(name,"restarts/restart_%07d.h5", istep);
  if (mpi_io_proc()) {
    fprintf(stderr,"writing restart file %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  #if (USE_MPI==TRUE)
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // first write some info about the current state of the run (e.g. time, istep, ...)
  file_grid_dims[0] = 1;
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? one : zero;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  hid_t memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, "NDIM", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  int ndim = NDIM;
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ndim);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Step", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &istep);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "NumProcs", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &numprocs);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time_Next_Dump", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_dump);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time_Next_Pdump", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_pdump);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Time_Next_Restart", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_restart);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Step_Next_Restart", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &i_next_restart);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Next_Dump_ID", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &dump_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "Next_Pdump_ID", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &pdump_id);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_CREATE);

  dset_id = H5Dcreate(file_id, "N_Tracer", H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  int ntracer = n_tracer_global;
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ntracer);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(filespace);
  H5Sclose(memspace);

  #if (GR_MONOPOLE==TRUE)
  /* dump lapse */
  file_grid_dims[0] = n1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "lapse", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse[0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  /* dump lapse edge */
  file_grid_dims[0] = n1+1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1+1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1+1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "lapse_edge", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if (mpi_io_proc()) status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse_edge[0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
  #endif /* GR_MONOPOLE==TRUE */


  /* DR: output all of the primitives at once */
  int mycount = get_prims_count();
  int * counts = malloc(numprocs*sizeof(int));
#if (USE_MPI==TRUE)
  MPI_Allgather(&mycount, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
#endif

  double * prims = malloc(mycount*sizeof(double));
  pack_prims(sim_p, prims);

  file_grid_dims[0] = 0;
  for(int ip = 0; ip < numprocs; ++ip) {
    file_grid_dims[0] += counts[ip];
  }
  file_start[0] = 0;
  for(int ip = 0; ip < myrank; ++ip) {
    file_start[0] += counts[ip];
  }
  file_count[0] = mycount;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(1, file_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "prims", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, prims);

  H5Pclose(plist_id);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  free(counts);
  free(prims);



  hid_t group_id = H5Gcreate(file_id, "/Tracers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
  file_grid_dims[0] = n_tracer_global;
  hsize_t ntrace = n_tracer_current;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = n_tracer_current;
  //H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, (const hsize_t *)&f_coord);
  //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  memspace  = H5Screate_simple(1, &ntrace, NULL);

  hsize_t *f_coord = malloc_rank1(n_tracer_current, sizeof(hsize_t));
  double *data = malloc_rank1(n_tracer_current, sizeof(double));
  char *flag = malloc_rank1(n_tracer_current, sizeof(char));

  i=0;
  tracer_t *tr = tracers;
  while(tr != NULL) {
      f_coord[i] = tr->id;
      data[i] = tr->x[0];
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/x", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  #if(NDIM>1)
  i=0;
  tr = tracers;
  while(tr != NULL) {
      data[i] = tr->x[1];
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/y", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  #endif
  #if(NDIM==3)
  i=0;
  tr = tracers;
  while(tr != NULL) {
      data[i] = tr->x[2];
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/z", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  #endif

  i=0;
  tr = tracers;
  while(tr != NULL) {
      data[i] = tr->mass;
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/mass", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  i=0;
  tr = tracers;
  while(tr != NULL) {
      flag[i] = tr->active;
      i++;
      tr = tr->next;
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "/Tracers/active", H5T_NATIVE_CHAR, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, f_coord);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  status = H5Dwrite(dset_id, H5T_NATIVE_CHAR, memspace, filespace, plist_id, flag);;
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(f_coord!=NULL) free(f_coord);
  if(data!=NULL) free(data);
  if(flag!=NULL) free(flag);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  TIMER_STOP;

  return;
}

int find_last_restart(char * const path, char filename[256])
{
  FTS * ftsp;
  FTSENT *p, *chp;
  char * const argv[2] = {path, NULL};
  int fts_options = FTS_COMFOLLOW | FTS_LOGICAL | FTS_NOCHDIR;
  char iter_str[10];
  int iter_max = -1;
  int iter;

  if((ftsp = fts_open(argv, fts_options, NULL)) == NULL) {
    warn("fts_open");
    return -1;
  }
  chp = fts_children(ftsp, 0);
  if(chp == NULL) {
    err(EINVAL, "no restart file found");
    return -1;
  }
  while((p = fts_read(ftsp)) != NULL) {
    if(FTS_F == p->fts_info) {
      if(sscanf(p->fts_name, "restart_%7[0-9].h5", iter_str)) {
        iter = atoi(iter_str);
        if(iter > iter_max) {
          iter_max = iter;
          strncpy(filename, p->fts_path, 256);
        }
      }
    }
  }
  fts_close(ftsp);
  return 0;
}

void restart_read()
{
  restart_read_v2();
}

void restart_read_v1()
{
  int ierr;
  int i,j,k,dd,ndim,p,g,ii,jj,kk,jstart,jstop,kstart;
  hsize_t file_grid_dims[NDIM+1],file_start[NDIM+1],file_count[NDIM+1];
  hsize_t mem_grid_dims[NDIM+1],mem_start[NDIM+1],one;
  char *name = "p";
  char slice[20];
  hid_t dset_id;
  herr_t status;
  double *pstart=NULL;

  TIMER_START("restart_read");

  if (restart_from_last) {
    if(mpi_io_proc()) {
      ierr = find_last_restart(restart_dir, restart_file);
      if(ierr != 0) {
        fprintf(stderr,"could not find restart file");
        exit(1);
      }
    }
    #if USE_MPI==TRUE
    MPI_Bcast(restart_file, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    #endif
  }

  if (mpi_io_proc()) {
    fprintf(stderr,"restarting from %s...", restart_file);
  }

  one = 1;

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  #if (USE_MPI==TRUE)
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  hid_t file_id = H5Fopen(restart_file, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  if (file_id < 0) {
    fprintf(stderr,"\n\nrestart failed.  does %s exist?\n\n", restart_file);
    exit(2);
  }

  // first read some info about the current state of the run (e.g. time, istep, ...)
  file_grid_dims[0] = 1;
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = 1;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  hid_t memspace  = H5Screate_simple(1, &one, NULL);

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "NDIM", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ndim);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&ndim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Step", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &istep);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&istep, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time_Next_Dump", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_dump);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t_next_dump, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  if (t + dt_dump < t_next_dump) t_next_dump = t+dt_dump;

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time_Next_Restart", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_restart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t_next_restart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  if (t + dt_restart < t_next_restart) t_next_restart = t+dt_restart;


  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Step_Next_Restart", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &i_next_restart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&i_next_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Next_Dump_ID", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &dump_id);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&dump_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time_Next_Pdump", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_pdump);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t_next_pdump, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  if (t + dt_pdump < t_next_pdump) t_next_pdump = t+dt_pdump;

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Next_Pdump_ID", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &pdump_id);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&pdump_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "N_Tracer", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    int ntrace;
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ntrace);
    n_tracer_global = ntrace;
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&n_tracer_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  H5Sclose(filespace);
  H5Sclose(memspace);

  #if (GR_MONOPOLE==TRUE)
  /* read lapse */
  file_grid_dims[0] = n1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "lapse", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse[0]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(gr_lapse, n1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  H5Sclose(filespace);
  H5Sclose(memspace);

  /* read lapse edge */
  file_grid_dims[0] = n1+1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1+1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1+1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "lapse_edge", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse_edge[0]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(gr_lapse_edge, n1+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  H5Sclose(filespace);
  H5Sclose(memspace);

  /* copy 1-d to 3-d */
  ZLOOP {
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
  }
  ii = istop[0];
  #if (NDIM>1)
  JSLOOP(ii,jj) {
  #endif
    #if (NDIM==3)
    KSLOOP(ii,jj,kk) {
    #endif
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
    #if (NDIM==3)
    }
    #endif
  #if (NDIM>1)
  }
  #endif

  #if (NDIM==2)
  for (ii=istart[0]; ii<istop[0]; ii++) {
    jj = JS(ii,istop[1]);
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
  }
  #endif

  #if (NDIM==3)
  for (ii=istart[0]; ii<istop[0]; ii++) {
    jj = JS(ii,istop[1]);
    KSLOOP(ii,jj,kk) {
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
    JSLOOP(ii,jj) {
      kk = KS(ii,jj,istop[2]);
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
  }
  #endif
  #endif /* GR_MONOPOLE==TRUE */

  #if (NDIM==1)
  p = 0;
  for (i=0; i<n1; i++) {
    file_grid_dims[0] = nvars;
    filespace = H5Screate_simple(1, file_grid_dims, NULL);
    if (i >= istart[0] && i < istop[0]) {
      jstart = JS(i,istart[1]);
      kstart = KS(i,jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,i,jstart,kstart,0);
      file_start[0] = 0;
      file_count[0] = nvars;
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    } else {
      jstart = JS(istart[0],istart[1]);
      kstart = KS(istart[0],jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,istart[0],jstart,kstart,0);
      file_count[0] = 0;
      H5Sselect_none(filespace);
    }
    mem_grid_dims[0] = nvars;
    memspace = H5Screate_simple(1, mem_grid_dims, NULL);
    mem_start[0] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

    sprintf(slice,"p%d", p++);
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, slice, plist_id);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, pstart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
  }
  #endif

  #if (NDIM==2)
  p = 0;
  for (i=0; i<n1; i++) {
    file_grid_dims[0] = global_grid_dims[1]/dj[i];
    file_grid_dims[1] = nvars;
    filespace = H5Screate_simple(2, file_grid_dims, NULL);
    if (i >= istart[0] && i < istop[0]) {
      jstart = JS(i,istart[1]);
      kstart = KS(i,jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,i,jstart,kstart,0);
      file_count[0] = my_grid_dims[1]/dj[i];
      file_count[1] = nvars;
      file_start[0] = jstart;
      file_start[1] = 0;
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
    } else {
      jstart = JS(istart[0],istart[1]);
      kstart = KS(istart[0],jstart,istart[2]);
      pstart = &NDP_ELEM_LINEAR(sim_p,istart[0],jstart,kstart,0);
      file_count[0] = 0;
      file_count[1] = 0;
      H5Sselect_none(filespace);
    }

    mem_grid_dims[0] = my_grid_dims[1]/dj[i] + 2*NG;
    mem_grid_dims[1] = nvars;
    memspace = H5Screate_simple(2, mem_grid_dims, NULL);
    mem_start[0] = 0;
    mem_start[1] = 0;
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

    sprintf(slice,"p%d", p++);
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, slice, plist_id);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &NDP_ELEM_LINEAR(sim_p,MIN(MAX(i,istart[0]),istop[0]-1),jstart,kstart,0));
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
  }
  #endif

  #if (NDIM==3)
  p = 0;
  for (i=0; i<n1; i++) {
    for (j=0; j<n2/dj[i]; j++) {
      file_grid_dims[0] = global_grid_dims[2]/dk[i][j];
      file_grid_dims[1] = nvars;
      filespace = H5Screate_simple(2, file_grid_dims, NULL);
      jstart = JS(i,istart[1]);
      jstop = JS(i,istop[1]);
      if (i >= istart[0] && i < istop[0] && j >= jstart && j < jstop) {
        kstart = KS(i,j,istart[2]);
        pstart = &NDP_ELEM_LINEAR(sim_p,i,j,kstart,0);
        file_count[0] = my_grid_dims[2]/dk[i][j];
        file_count[1] = nvars;
        file_start[0] = kstart;
        file_start[1] = 0;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
      } else {
        jstart = JS(istart[0],istart[1]);
        kstart = KS(istart[0],jstart,istart[2]);
        pstart = &NDP_ELEM_LINEAR(sim_p,istart[0],jstart,kstart,0);
        file_count[0] = 0;
        file_count[1] = 0;
        H5Sselect_none(filespace);
      }

      mem_grid_dims[0] = my_grid_dims[2]/dk[i][j];
      mem_grid_dims[1] = nvars;
      memspace = H5Screate_simple(2, mem_grid_dims, NULL);
      mem_start[0] = 0;
      mem_start[1] = 0;
      H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

      sprintf(slice,"p%d", p++);
      plist_id = H5Pcreate(H5P_DATASET_ACCESS);
      dset_id = H5Dopen(file_id, slice, plist_id);
      H5Pclose(plist_id);

      plist_id = H5Pcreate(H5P_DATASET_XFER);
      #if (USE_MPI==TRUE)
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      #endif
      status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, pstart);
      H5Dclose(dset_id);
      H5Pclose(plist_id);

      H5Sclose(memspace);
      H5Sclose(filespace);
    }
  }
  #endif

  // now get tracers
#define N_TRACER_RESTART_BATCH 200000
  double xbatch[NDIM][N_TRACER_RESTART_BATCH];
  double x[NDIM];
  int index[NDIM];

  file_grid_dims[0] = n_tracer_global;
  int nread = 0;
  int nbatch;
  int nbuff = 100;
  hsize_t *pid = malloc_rank1(nbuff, sizeof(hsize_t));
  int current_id = 0;
  while(nread < n_tracer_global) {
      nbatch = MIN(N_TRACER_RESTART_BATCH, n_tracer_global-nread);
      if(mpi_io_proc()) {
        fprintf(stderr,"Reading %d %d\n", nread, current_id);
        hsize_t ntrace = n_tracer_current;
        filespace = H5Screate_simple(1, file_grid_dims, NULL);
        file_start[0] = nread;
        file_count[0] = nbatch;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

        memspace  = H5Screate_simple(1, file_count, NULL);
        mem_grid_dims[0] = file_count[0];
        memspace = H5Screate_simple(1, mem_grid_dims, NULL);
        mem_start[0] = 0;
        H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

        plist_id = H5Pcreate(H5P_DATASET_ACCESS);
        dset_id = H5Dopen(file_id, "/Tracers/x", plist_id);
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        #if (USE_MPI==TRUE)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        #endif
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &xbatch[0][0]);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        #if(NDIM>1)
        plist_id = H5Pcreate(H5P_DATASET_ACCESS);
        dset_id = H5Dopen(file_id, "/Tracers/y", plist_id);
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        #if (USE_MPI==TRUE)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        #endif
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &xbatch[1][0]);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        #endif
        #if(NDIM==3)
        plist_id = H5Pcreate(H5P_DATASET_ACCESS);
        dset_id = H5Dopen(file_id, "/Tracers/z", plist_id);
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        #if (USE_MPI==TRUE)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        #endif
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &xbatch[2][0]);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        #endif
      }
      MPI_Bcast(xbatch, NDIM*N_TRACER_RESTART_BATCH, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // now check for tracers on my rank
      for(i=0;i<nbatch;i++) {
        DLOOP x[dd] = xbatch[dd][i];
        x_to_ijk(x, index);
        if(!is_outside_rank(index, myrank)) { // this one lives on my rank
            if(current_id == nbuff) {
                nbuff += 100;
                pid = realloc(pid, nbuff*sizeof(hsize_t));
            }
            pid[current_id] = nread+i;
            current_id++;
        }
      }

      nread += nbatch;
  }
#undef N_TRACER_RESTART_BATCH
  n_tracer_current = current_id;

  // now make the linked list so we can read in and populate the tracer data
  tracers = NULL;
  for(i=current_id-1; i >=0; i--) {
    tracer_t *new = malloc_rank1(1, sizeof(tracer_t));
    new->next = tracers;
    new->id = pid[i];
    tracers = new;
  }

  // tracer id's are all as they should be now, so set the flag saying so
  consistent_tracer_ids = 1;

  // now have everyone read in the data they need
  file_grid_dims[0] = n_tracer_global;
  hsize_t ntrace = n_tracer_current;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = n_tracer_current;
  memspace  = H5Screate_simple(1, &ntrace, NULL);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, pid);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  double *data = malloc_rank1(n_tracer_current, sizeof(double));
  char *flag = malloc_rank1(n_tracer_current, sizeof(char));

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/x", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  if(status < 0) {
    fprintf(stderr,"Error on %d reading in %d Tracers/x\n", myrank, n_tracer_current);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tracer_t *tr = tracers;
  while(tr != NULL) {
    tr->x[0] = data[i];
    tr = tr->next;
    i++;
  }
  #if(NDIM>1)
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/y", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  if(status < 0) {
    fprintf(stderr,"Error on %d reading in %d Tracers/y\n", myrank, n_tracer_current);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->x[1] = data[i];
    tr = tr->next;
    i++;
  }
  #endif
  #if(NDIM==3)
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/z", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->x[2] = data[i];
    tr = tr->next;
    i++;
  }
  #endif
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/mass", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->mass = data[i];
    tr = tr->next;
    i++;
  }
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/active", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_CHAR, memspace, filespace, plist_id, flag);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->active = flag[i];
    tr = tr->next;
    i++;
  }

  free(pid);
  if(data!=NULL) free(data);
  if(flag!=NULL) free(flag);

  //check if all my tracers are really mine...
  tr = tracers;
  while(tr != NULL) {
    x_to_ijk(tr->x, index);
    if(is_outside_rank(index, myrank)) {
        fprintf(stderr,"Found a problem with read in tracer\n");
    }
    tr = tr->next;
  }

  H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_io_proc()) fprintf(stderr,"done\n");

  TIMER_STOP;
  return;
}

void restart_read_v2()
{
  int ierr;
  int i,j,k,dd,ndim,p,g,ii,jj,kk,jstart,jstop,kstart;
  hsize_t file_grid_dims[NDIM+1],file_start[NDIM+1],file_count[NDIM+1];
  hsize_t mem_grid_dims[NDIM+1],mem_start[NDIM+1],one;
  char *name = "p";
  char slice[20];
  hid_t dset_id;
  herr_t status;
  double *pstart=NULL;

  TIMER_START("restart_read");

  if (restart_from_last) {
    if(mpi_io_proc()) {
      ierr = find_last_restart(restart_dir, restart_file);
      if(ierr != 0) {
        fprintf(stderr,"could not find restart file");
        exit(1);
      }
    }
    #if USE_MPI==TRUE
    MPI_Bcast(restart_file, 256, MPI_CHAR, 0, MPI_COMM_WORLD);
    #endif
  }

  if (mpi_io_proc()) {
    fprintf(stderr,"restarting from %s...", restart_file);
  }

  one = 1;

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  #if (USE_MPI==TRUE)
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  #endif
  hid_t file_id = H5Fopen(restart_file, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  if (file_id < 0) {
    fprintf(stderr,"\n\nrestart failed.  does %s exist?\n\n", restart_file);
    exit(2);
  }

  // first read some info about the current state of the run (e.g. time, istep, ...)
  file_grid_dims[0] = 1;
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = 1;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
  hid_t memspace  = H5Screate_simple(1, &one, NULL);

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "NDIM", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ndim);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&ndim, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Step", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &istep);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&istep, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "NumProcs", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    int rnumprocs;
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &rnumprocs);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
    if (rnumprocs != numprocs) {
      fprintf(stdout, "Trying to restart using \"%d\" processes, but restart file was created using \"%d\" processes. This is not supported.", rnumprocs, numprocs);
      fflush(stdout);
      abort();
    }
  }

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time_Next_Dump", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_dump);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t_next_dump, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  if (t + dt_dump < t_next_dump) t_next_dump = t+dt_dump;

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time_Next_Restart", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_restart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t_next_restart, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  if (t + dt_restart < t_next_restart) t_next_restart = t+dt_restart;


  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Step_Next_Restart", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &i_next_restart);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&i_next_restart, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Next_Dump_ID", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &dump_id);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&dump_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Time_Next_Pdump", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t_next_pdump);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&t_next_pdump, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  if (t + dt_pdump < t_next_pdump) t_next_pdump = t+dt_pdump;

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "Next_Pdump_ID", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &pdump_id);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&pdump_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "N_Tracer", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    int ntrace;
    status = H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &ntrace);
    n_tracer_global = ntrace;
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&n_tracer_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
  #endif

  H5Sclose(filespace);
  H5Sclose(memspace);

  #if (GR_MONOPOLE==TRUE)
  /* read lapse */
  file_grid_dims[0] = n1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "lapse", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse[0]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(gr_lapse, n1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  H5Sclose(filespace);
  H5Sclose(memspace);

  /* read lapse edge */
  file_grid_dims[0] = n1+1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = mpi_io_proc() ? n1+1 : 0;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  mem_grid_dims[0] = mpi_io_proc() ? n1+1 : 0;
  memspace = H5Screate_simple(1, mem_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  if (mpi_io_proc()) {
    plist_id = H5Pcreate(H5P_DATASET_ACCESS);
    dset_id = H5Dopen(file_id, "lapse_edge", plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    #if (USE_MPI==TRUE)
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    #endif
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &gr_lapse_edge[0]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(gr_lapse_edge, n1+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  #endif
  H5Sclose(filespace);
  H5Sclose(memspace);

  /* copy 1-d to 3-d */
  ZLOOP {
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
  }
  ii = istop[0];
  #if (NDIM>1)
  JSLOOP(ii,jj) {
  #endif
    #if (NDIM==3)
    KSLOOP(ii,jj,kk) {
    #endif
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
    #if (NDIM==3)
    }
    #endif
  #if (NDIM>1)
  }
  #endif

  #if (NDIM==2)
  for (ii=istart[0]; ii<istop[0]; ii++) {
    jj = JS(ii,istop[1]);
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
  }
  #endif

  #if (NDIM==3)
  for (ii=istart[0]; ii<istop[0]; ii++) {
    jj = JS(ii,istop[1]);
    KSLOOP(ii,jj,kk) {
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
    JSLOOP(ii,jj) {
      kk = KS(ii,jj,istop[2]);
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
  }
  #endif
  #endif /* GR_MONOPOLE==TRUE */



  /* DR: read serialized 3D data and unpack */
  int mycount = get_prims_count();
  int * counts = malloc(numprocs*sizeof(int));
  MPI_Allgather(&mycount, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);

  double * prims = malloc(mycount*sizeof(double));

  file_grid_dims[0] = 0;
  for(int ip = 0; ip < numprocs; ++ip) {
    file_grid_dims[0] += counts[ip];
  }
  file_start[0] = 0;
  for(int ip = 0; ip < myrank; ++ip) {
    file_start[0] += counts[ip];
  }
  file_count[0] = mycount;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  memspace = H5Screate_simple(1, file_grid_dims, NULL);
  mem_start[0] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "prims", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, prims);

  H5Pclose(plist_id);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

  unpack_prims(prims, sim_p);

  free(counts);
  free(prims);



  // now get tracers
#define N_TRACER_RESTART_BATCH 200000
  double xbatch[NDIM][N_TRACER_RESTART_BATCH];
  double x[NDIM];
  int index[NDIM];

  file_grid_dims[0] = n_tracer_global;
  int nread = 0;
  int nbatch;
  int nbuff = 100;
  hsize_t *pid = malloc_rank1(nbuff, sizeof(hsize_t));
  int current_id = 0;
  while(nread < n_tracer_global) {
      nbatch = MIN(N_TRACER_RESTART_BATCH, n_tracer_global-nread);
      if(mpi_io_proc()) {
        fprintf(stderr,"Reading %d %d\n", nread, current_id);
        hsize_t ntrace = n_tracer_current;
        filespace = H5Screate_simple(1, file_grid_dims, NULL);
        file_start[0] = nread;
        file_count[0] = nbatch;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

        memspace  = H5Screate_simple(1, file_count, NULL);
        mem_grid_dims[0] = file_count[0];
        memspace = H5Screate_simple(1, mem_grid_dims, NULL);
        mem_start[0] = 0;
        H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

        plist_id = H5Pcreate(H5P_DATASET_ACCESS);
        dset_id = H5Dopen(file_id, "/Tracers/x", plist_id);
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        #if (USE_MPI==TRUE)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        #endif
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &xbatch[0][0]);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        #if(NDIM>1)
        plist_id = H5Pcreate(H5P_DATASET_ACCESS);
        dset_id = H5Dopen(file_id, "/Tracers/y", plist_id);
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        #if (USE_MPI==TRUE)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        #endif
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &xbatch[1][0]);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        #endif
        #if(NDIM==3)
        plist_id = H5Pcreate(H5P_DATASET_ACCESS);
        dset_id = H5Dopen(file_id, "/Tracers/z", plist_id);
        H5Pclose(plist_id);
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        #if (USE_MPI==TRUE)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
        #endif
        status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &xbatch[2][0]);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        #endif
      }
      MPI_Bcast(xbatch, NDIM*N_TRACER_RESTART_BATCH, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      // now check for tracers on my rank
      for(i=0;i<nbatch;i++) {
        DLOOP x[dd] = xbatch[dd][i];
        x_to_ijk(x, index);
        if(!is_outside_rank(index, myrank)) { // this one lives on my rank
            if(current_id == nbuff) {
                nbuff += 100;
                pid = realloc(pid, nbuff*sizeof(hsize_t));
            }
            pid[current_id] = nread+i;
            current_id++;
        }
      }

      nread += nbatch;
  }
#undef N_TRACER_RESTART_BATCH
  n_tracer_current = current_id;

  // now make the linked list so we can read in and populate the tracer data
  tracers = NULL;
  for(i=current_id-1; i >=0; i--) {
    tracer_t *new = malloc_rank1(1, sizeof(tracer_t));
    new->next = tracers;
    new->id = pid[i];
    tracers = new;
  }

  // tracer id's are all as they should be now, so set the flag saying so
  consistent_tracer_ids = 1;

  // now have everyone read in the data they need
  file_grid_dims[0] = n_tracer_global;
  hsize_t ntrace = n_tracer_current;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = n_tracer_current;
  memspace  = H5Screate_simple(1, &ntrace, NULL);
  if(n_tracer_current) H5Sselect_elements(filespace, H5S_SELECT_SET, n_tracer_current, pid);
  else H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  double *data = malloc_rank1(n_tracer_current, sizeof(double));
  char *flag = malloc_rank1(n_tracer_current, sizeof(char));

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/x", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  if(status < 0) {
    fprintf(stderr,"Error on %d reading in %d Tracers/x\n", myrank, n_tracer_current);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tracer_t *tr = tracers;
  while(tr != NULL) {
    tr->x[0] = data[i];
    tr = tr->next;
    i++;
  }
  #if(NDIM>1)
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/y", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  if(status < 0) {
    fprintf(stderr,"Error on %d reading in %d Tracers/y\n", myrank, n_tracer_current);
  }
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->x[1] = data[i];
    tr = tr->next;
    i++;
  }
  #endif
  #if(NDIM==3)
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/z", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->x[2] = data[i];
    tr = tr->next;
    i++;
  }
  #endif
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/mass", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->mass = data[i];
    tr = tr->next;
    i++;
  }
  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  dset_id = H5Dopen(file_id, "/Tracers/active", plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  #if (USE_MPI==TRUE)
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  #endif
  status = H5Dread(dset_id, H5T_NATIVE_CHAR, memspace, filespace, plist_id, flag);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  i=0;
  tr = tracers;
  while(tr != NULL) {
    tr->active = flag[i];
    tr = tr->next;
    i++;
  }

  free(pid);
  if(data!=NULL) free(data);
  if(flag!=NULL) free(flag);

  //check if all my tracers are really mine...
  tr = tracers;
  while(tr != NULL) {
    x_to_ijk(tr->x, index);
    if(is_outside_rank(index, myrank)) {
        fprintf(stderr,"Found a problem with read in tracer\n");
    }
    tr = tr->next;
  }

  H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);
  if (mpi_io_proc()) fprintf(stderr,"done\n");

  TIMER_STOP;
  return;
}


FILE *mpi_fopen(char *name, const char *mode)
{
  FILE *fp=NULL;

  if (mpi_io_proc()) {
    fp = fopen(name,mode);
    if (fp == NULL) {
        fprintf(stderr,"failed to open input file %s\n", name);
        exit(123);
    }
  }
  return fp;
}

void mpi_fclose(FILE *fp)
{
  if (mpi_io_proc()) fclose(fp);
}

int mpi_fgets(char *line, int num, FILE *fp, MPI_Comm comm)
{
  int status = 1;
  int nbuf;

  if (mpi_io_proc()) {
    if (fgets(line, num, fp)==NULL) status = 0;
    nbuf = strlen(line)+1;
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&status, 1, MPI_INT, 0, comm);
  MPI_Bcast(&nbuf, 1, MPI_INT, 0, comm);
  // if (status) MPI_Bcast(line, num, MPI_CHAR, 0, comm);
  if (status) MPI_Bcast(line, nbuf, MPI_CHAR, 0, comm);
  #endif

  return status;
}

#define VALUE_TO_STRING(x) #x
#define VALUE(x) VALUE_TO_STRING(x)
#define VAR_NAME_VALUE(var) #var "="  VALUE(var)

void parse_input(char *name)
{
  int pos,cpos;
  // char line[2048],lhs[512],rhs[512];
  int nbuf=2048;
  char *line=NULL, *lhs=NULL, *rhs=NULL;
  char *eq, *com;
  FILE *fp;

  TIMER_START("parse_input");

  line = malloc_rank1(nbuf+1, sizeof *line);
  lhs  = malloc_rank1(nbuf+1, sizeof *lhs);
  rhs  = malloc_rank1(nbuf+1, sizeof *rhs);
  for (int i=0; i<nbuf+1; i++) line[i] = '#';  // done to avoid an "uninitialised byte(s)" error in valgrind
  for (int i=0; i<nbuf+1; i++) lhs[i] = '#';  // done to avoid an "uninitialised byte(s)" error in valgrind
  for (int i=0; i<nbuf+1; i++) rhs[i] = '#';  // done to avoid an "uninitialised byte(s)" error in valgrind

  // set some (sensible?) defaults
  //tianshu
  rhocut = 1e13;
  use_kom = 1;
  kom_dtmin = 1e-12;
  kom_epsilon = 1e-2;
  kom_delta = 1e-2;
  Ecut1 = 100;
  Ecut2 = 50;
  Ecut3 = 50;
  Nskip = 1;
  rx_info.facta = 10*1e5;
  rx_info.factb = 15*1e5;
  rx_info.factc = 5*1e5;
  //thx_info.nice_alpha = 1.0007;  // the closer to 1, the closer to constant \Delta \mu
  thx_info.nice_alpha = 1.01;  // much below 1.0007, Romberg has trouble converging
                             // though that may not actually matter
  thx_info.poly_xt = 1.2;
  thx_info.poly_alpha = 4;
  //end tianshu

  cfl = 0.9;
  n1 = 1024;
  n2 = 128;
  n3 = 256;
  nr1 = 12;
  nr2 = 12;
  nr3 = 12;
  n_tracer_target = 0;
  ccsn_dr_min = 5.0e4;
  strcpy(grid_units,"cm");
  emin1 = 1.0;
  emin2 = 1.0;
  emin3 = 1.0;
  emax1 = 300.0;
  emax2 = 100.0;
  emax3 = 100.0;
  force_flux_limit = 0;
  dt_grow_safety = 1.1;
  dt_init_safety = 0.1;
  dtmax = 1.0e200;
  dtmin = 1.0e-200;
  strcpy(freq_type,"mev");
  ncomp = 0;
  restart_create_1d = 0;  // flag to create 1-d ASCII file for restart
  restart_from_1d = 0;  // flag to restart from 1-d ASCII file
  restart_from_3d = 0;  // flag to restart from a 3-d HDF5 dump
  restart_from_hdf = 0;  // flag to restart from HDF5 dump
  restart_from_last = 0;  // flag to ignore input file and use the most recent restart dump instead
  restart_delay = 0.010;  // post-bounce 1-d ASCII restart file delay [s]
  strcpy(restart_file, "");
  strcpy(restart_dir, ".");
  dump_hdf = 1;
  dump_rad_vars = 1;  // flag to dump radiation variables (E[g],F1[g],F2[g],F3[g]) in all groups
  detect_tbounce = 0;  // flag to detect and record the time of core bounce
  implicit_err_tol = 1.0e-10;
  max_implicit_iter = 20;
  hydro_interp_order = 3;
  rad_interp_order = 2;
  multipole_lmax = 12;  // maximum multipole order
  dT_tol = 1.0e-2;
  du_tol = 1.0e-2;
  dye_tol = 1.0e-2;
  initial_dt = 1.0e100;
  tmax = 0.0;
  dt_dump = 1.0e100;
  dt_pdump = dt_dump;
  dt_restart = 1.0e100;
  di_restart = INT_MAX/2;
  nstep_log = 100;  // number of steps between calls to run_log()
  nstep_analysis = 100;  // number of steps between calls to analysis_inloop()
  nstep_timing = INT_MAX;  // frequency of the timing output
  chat = CLIGHT;
  chat_safety = 10.0;
  use_chat = 0;
  dx[0] = 1.0;
  dx[1] = 1.0;
  dx[2] = 1.0;
  max_steps = 2147483647;
  max_wtime = 1e100;  // maximum wall time for run; when reached, will force restart dump
  mass_inside_tracers = 0.;
  mass_outside_tracers = 0.9 * DBL_MAX / MSUN;
  L_enu = 2.0e52;
  T_enu = 4.0;
  e_floor = -1.0;  // minimum value of internal (thermal) energy
  rho_floor = -1.0;  // minimum value of density
  perturb_r1m = 1000e5; // velocity perturbation, region 1, rmin [cm]
  perturb_r1p = 1300e5; // velocity perturbation, region 1, rmax [cm]
  perturb_l1 = 12; // velocity perturbation lmax, region 1
  perturb_m1 = 1;  // velocity perturbation mmax, region 1
  perturb_n1 = 1; // velocity perturbation nmax, region 1
  perturb_dv1 = 0.0; // velocity perturbation amplitude [cm/s], region 1
  perturb_r2m = 1600e5; // velocity perturbation, region 2, rmin [cm]
  perturb_r2p = 2200e5; // velocity perturbation, region 2, rmax [cm]
  perturb_l2 = 10; // velocity perturbation lmax, region 2
  perturb_m2 = 1;  // velocity perturbation mmax, region 2
  perturb_n2 = 1; // velocity perturbation nmax, region 2
  perturb_dv2 = 0.0; // velocity perturbation amplitude [cm/s], region 2
  perturb_r3m = 4000e5; // velocity perturbation, region 3, rmin [cm]
  perturb_r3p = 17000e5; // velocity perturbation, region 3, rmax [cm]
  perturb_l3 = 4; // velocity perturbation lmax, region 3
  perturb_m3 = 1;  // velocity perturbation mmax, region 3
  perturb_n3 = 1; // velocity perturbation nmax, region 3
  perturb_dv3 = 0.0; // velocity perturbation amplitude [cm/s], region 3
  perturb_level = 0.1;  // relative level [%] for random density perturbation
  perturb_delay = 0.010;  // post-bounce perturbation delay [s]; perturb immediately if <0
  rotate_Omega0 = 1.0;  // Omega0 [s^-1] for rotating model:  Omega(r) = Omega0/(1 + (r/A)^2)
  rotate_A = 1000.0;  // A [km] for rotating model:  Omega(r) = Omega0/(1 + (r/A)^2)
  outer_radius = 1.0e9; // outer boundary radius for the CCSN problem in cm
  decomp_only = 0;  // flag to test decomposition only without allocating or running
  decomp_npmin = 4;  // test processor decomposition in range npmin:npskip:npmax (shouldn't go below 4)
  decomp_npmax = 100;  // test processor decomposition in range npmin:npskip:npmax
  decomp_npskip = 1;  // test processor decomposition in range npmin:npskip:npmax
  decomp_from_file = 0;  // flag to load decomposition proc info from file ("proc_info.txt")
  strcpy(decomp_path, ".");
  include_inelastic = 1;
  strcpy(opac_param_file,"MUST SPECIFY OPAC PARAM FILE NAME VIA opac_param_file\n");
  strcpy(opac_file,"MUST SPECIFY OPAC FILE NAME VIA opac_file\n");
  strcpy(eos_file,"MUST SPECIFY EOS FILE NAME VIA eos_file\n");
  strcpy(inelastic_root,"");

  fp = mpi_fopen(name,"r");

  eq = NULL;
  com = NULL;

  // Dump Makefile configuration
  if (mpi_io_proc()) {
    #if (NDIM==1)
    fprintf(stderr,"NDIM=1\n");
    #elif (NDIM==2)
    fprintf(stderr,"NDIM=2\n");
    #elif (NDIM==3)
    fprintf(stderr,"NDIM=3\n");
    #else
    #error NDIM VALUE UNKNOWN!
    #endif

    #if (USE_MPI==TRUE)
    fprintf(stderr,"USE_MPI=TRUE\n");
    #elif (USE_MPI==FALSE)
    fprintf(stderr,"USE_MPI=FALSE\n");
    #else
    #error USE_MPI VALUE UNKNOWN!
    #endif

    #if (USE_OMP==TRUE)
    fprintf(stderr,"USE_OMP=TRUE\n");
    #elif (USE_OMP==FALSE)
    fprintf(stderr,"USE_OMP=FALSE\n");
    #else
    #error USE_OMP VALUE UNKNOWN!
    #endif

    #if (GEOM==CARTESIAN)
    fprintf(stderr,"GEOM=CARTESIAN\n");
    #elif (GEOM==SPHERICAL)
    fprintf(stderr,"GEOM=SPHERICAL\n");
    #elif (GEOM==CYLINDRICAL)
    fprintf(stderr,"GEOM=CYLINDRICAL\n");
    #else
    #error GEOM VALUE UNKNOWN!
    #endif

    #if (EOS==GAMMA_LAW)
    fprintf(stderr,"EOS=GAMMA_LAW\n");
    #elif (EOS==POLYTROPIC)
    fprintf(stderr,"EOS=POLYTROPIC\n");
    #elif (EOS==COLLAPSE)
    fprintf(stderr,"EOS=COLLAPSE\n");
    #else
    #error EOS VALUE UNKNOWN!
    #endif

    #if (GRAV==NO_GRAV)
    fprintf(stderr,"GRAV=NO_GRAV\n");
    #elif (GRAV==FIXED_GRAV)
    fprintf(stderr,"GRAV=FIXED_GRAV\n");
    #elif (GRAV==PRESCRIBED_GRAV)
    fprintf(stderr,"GRAV=PRESCRIBED_GRAV\n");
    #elif (GRAV==SPHERICAL_MONOPOLE_GRAV)
    fprintf(stderr,"GRAV=SPHERICAL_MONOPOLE_GRAV\n");
    #elif (GRAV==SPHERICAL_MULTIPOLE_GRAV)
    fprintf(stderr,"GRAV=SPHERICAL_MULTIPOLE_GRAV\n");
    #elif (GRAV==USER_GRAV)
    fprintf(stderr,"GRAV=USER_GRAV\n");
    #else
    #error GRAV VALUE UNKNOWN!
    #endif

    #if (GR_MONOPOLE==TRUE)
    fprintf(stderr,"GR_MONOPOLE=TRUE\n");
    #elif (GR_MONOPOLE==FALSE)
    fprintf(stderr,"GR_MONOPOLE=FALSE\n");
    #else
    #error GR_MONOPOLE VALUE UNKNOWN!
    #endif

    #if (PN_POTENTIAL==TRUE)
    fprintf(stderr,"PN_POTENTIAL=TRUE\n");
    #elif (PN_POTENTIAL==FALSE)
    fprintf(stderr,"PN_POTENTIAL=FALSE\n");
    #else
    #error PN_POTENTIAL VALUE UNKNOWN!
    #endif

    #if (USE_EXT_SRC==TRUE)
    fprintf(stderr,"USE_EXT_SRC=TRUE\n");
    #elif (USE_EXT_SRC==FALSE)
    fprintf(stderr,"USE_EXT_SRC=FALSE\n");
    #else
    #error USE_EXT_SRC VALUE UNKNOWN!
    #endif

    #if (USE_AUX_SRC==TRUE)
    fprintf(stderr,"USE_AUX_SRC=TRUE\n");
    #elif (USE_AUX_SRC==FALSE)
    fprintf(stderr,"USE_AUX_SRC=FALSE\n");
    #else
    #error USE_AUX_SRC VALUE UNKNOWN!
    #endif

    #if (DO_HYDRO==TRUE)
    fprintf(stderr,"DO_HYDRO=TRUE\n");
    #elif (DO_HYDRO==FALSE)
    fprintf(stderr,"DO_HYDRO=FALSE\n");
    #else
    #error DO_HYDRO VALUE UNKNOWN!
    #endif

    #if (PRINT_OPAC_DATA==TRUE)
    fprintf(stderr,"PRINT_OPAC_DATA=TRUE\n");
    #elif (PRINT_OPAC_DATA==FALSE)
    fprintf(stderr,"PRINT_OPAC_DATA=FALSE\n");
    #else
    #error PRINT_OPAC_DATA VALUE UNKNOWN!
    #endif

    #if (DO_RADIATION==TRUE)
    fprintf(stderr,"DO_RADIATION=TRUE\n");
    #elif (DO_RADIATION==FALSE)
    fprintf(stderr,"DO_RADIATION=FALSE\n");
    #else
    #error DO_RADIATION VALUE UNKNOWN!
    #endif

    #if (PHOTON==TRUE)
    fprintf(stderr,"PHOTON=TRUE\n");
    #elif (PHOTON==FALSE)
    fprintf(stderr,"PHOTON=FALSE\n");
    #else
    #error PHOTON VALUE UNKNOWN!
    #endif

    #if (NEUTRINO==TRUE)
    fprintf(stderr,"NEUTRINO=TRUE\n");
    #elif (NEUTRINO==FALSE)
    fprintf(stderr,"NEUTRINO=FALSE\n");
    #else
    #error NEUTRINO VALUE UNKNOWN!
    #endif

    #if (PERTURB==NONE)
    fprintf(stderr,"PERTURB=NONE\n");
    #elif (PERTURB==VELOCITY_SPHERICAL_HARMONIC)
    fprintf(stderr,"PERTURB=VELOCITY_SPHERICAL_HARMONIC\n");
    #elif (PERTURB==DENSITY_RANDOM)
    fprintf(stderr,"PERTURB=DENSITY_RANDOM\n");
    #else
    #error PERTURB VALUE UNKNOWN!
    #endif

    #if (DO_ROTATION==TRUE)
    fprintf(stderr,"DO_ROTATION=TRUE\n");
    #elif (DO_ROTATION==FALSE)
    fprintf(stderr,"DO_ROTATION=FALSE\n");
    #else
    #error DO_ROTATION VALUE UNKNOWN!
    #endif

    #if (ENFORCE_FLOORS==TRUE)
    fprintf(stderr,"ENFORCE_FLOORS=TRUE\n");
    #elif (ENFORCE_FLOORS==FALSE)
    fprintf(stderr,"ENFORCE_FLOORS=FALSE\n");
    #else
    #error ENFORCE_FLOORS VALUE UNKNOWN!
    #endif

    #if (DENDRITIC_GRID==TRUE)
    fprintf(stderr,"DENDRITIC_GRID=TRUE\n");
    #elif (DENDRITIC_GRID==FALSE)
    fprintf(stderr,"DENDRITIC_GRID=FALSE\n");
    #else
    #error DENDRITIC_GRID VALUE UNKNOWN!
    #endif

    #if (OUTPUT_PRECISION==SINGLE)
    fprintf(stderr,"OUTPUT_PRECISION=SINGLE\n");
    #elif (OUTPUT_PRECISION==DOUBLE)
    fprintf(stderr,"OUTPUT_PRECISION=DOUBLE\n");
    #else
    #error OUTPUT_PRECISION VALUE UNKNOWN!
    #endif

    #if (DO_TIMING==TRUE)
    fprintf(stderr,"DO_TIMING=TRUE\n");
    #elif (DO_TIMING==FALSE)
    fprintf(stderr,"DO_TIMING=FALSE\n");
    #else
    #error DO_TIMING VALUE UNKNOWN!
    #endif
  }

  while (mpi_fgets(line,nbuf,fp,MPI_COMM_WORLD)) {
    eq = strchr(line,'=');
    if (eq == NULL) continue;
    pos = eq - line;
    strncpy(lhs,line,pos);
    lhs[pos] = 0;
    char *trim_lhs = trim_whitespace(lhs);

    com = strchr(line,'#');
    if (com == NULL) {
      strcpy(rhs,&line[pos+1]);
      rhs[(int)strlen(line)-pos-1] = 0;
    } else {
      cpos = com-line;
      if (cpos < pos) continue;
      strncpy(rhs,&line[pos+1],cpos-pos-1);
      rhs[cpos-pos-1] = 0;
    }
    char *trim_rhs = trim_whitespace(rhs);

    if (strcmp(trim_lhs,"cfl")==0) {
      cfl = atof(trim_rhs);
      if (cfl >= 1.0) {
        if (mpi_io_proc()) fprintf(stderr,"WARNING: cfl >= 1, THIS WON'T WORK!\n");
        //exit(123);
      }
      if (mpi_io_proc()) fprintf(stderr,"Set cfl = %f\n", cfl);
    } else if (strcmp(trim_lhs,"restart_create_1d")==0) {
      restart_create_1d = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_create_1d = %d\n", restart_create_1d);
    } else if (strcmp(trim_lhs,"restart_from_1d")==0) {
      restart_from_1d = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_from_1d = %d\n", restart_from_1d);
    } else if (strcmp(trim_lhs,"restart_from_3d")==0) {
      restart_from_3d = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_from_3d = %d\n", restart_from_3d);
    } else if (strcmp(trim_lhs,"restart_from_hdf")==0) {
      restart_from_hdf = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_from_hdf = %d\n", restart_from_hdf);
    } else if (strcmp(trim_lhs,"restart_from_last")==0) {
      restart_from_last = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_from_last = %d\n", restart_from_last);
    } else if (strcmp(trim_lhs,"restart_delay")==0) {
      restart_delay = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_delay = %g\n", restart_delay);
    } else if (strcmp(trim_lhs,"multipole_lmax")==0) {
      multipole_lmax = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set multipole_lmax = %d\n", multipole_lmax);
    } else if (strcmp(trim_lhs,"n1")==0) {
      n1 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set n1 = %d\n", n1);
    } else if (strcmp(trim_lhs,"n2")==0) {
      n2 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set n2 = %d\n", n2);
    } else if (strcmp(trim_lhs,"n3")==0) {
      n3 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set n3 = %d\n", n3);
    } else if (strcmp(trim_lhs,"nr1")==0) {
      nr1 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set nr1 = %d\n", nr1);
    } else if (strcmp(trim_lhs,"n_tracer_target")==0) {
      n_tracer_target = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set n_tracer_target = %d\n", n_tracer_target);
    } else if(strcmp(trim_lhs,"ccsn_dr_min")==0) {
        ccsn_dr_min = atof(trim_rhs);
        if(mpi_io_proc()) fprintf(stderr,"Set ccsn_dr_min = %g\n", ccsn_dr_min);
    } else if (strcmp(trim_lhs,"nr2")==0) {
      nr2 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set nr2 = %d\n", nr2);
    } else if (strcmp(trim_lhs,"nr3")==0) {
      nr3 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set nr3 = %d\n", nr3);
    } else if (strcmp(trim_lhs,"emin1")==0) {
      emin1 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set emin1 = %g\n", emin1);
    } else if (strcmp(trim_lhs,"emin2")==0) {
      emin2 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set emin2 = %g\n", emin2);
    } else if (strcmp(trim_lhs,"emin3")==0) {
      emin3 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set emin3 = %g\n", emin3);
    } else if (strcmp(trim_lhs,"emax1")==0) {
      emax1 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set emax1 = %g\n", emax1);
    } else if (strcmp(trim_lhs,"emax2")==0) {
      emax2 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set emax2 = %g\n", emax2);
    } else if (strcmp(trim_lhs,"emax3")==0) {
      emax3 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set emax3 = %g\n", emax3);
    } else if (strcmp(trim_lhs,"force_flux_limit")==0) {
      force_flux_limit = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set force_flux_limit = %d\n", force_flux_limit);
    } else if (strcmp(trim_lhs,"dt_grow_safety")==0) {
      dt_grow_safety = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_grow_safety = %f\n", dt_grow_safety);
    } else if (strcmp(trim_lhs,"dt_init_safety")==0) {
      dt_init_safety = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_init_safety = %f\n", dt_init_safety);
    } else if (strcmp(trim_lhs,"dt_max")==0) {
      dtmax = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_max = %g\n", dtmax);
    } else if (strcmp(trim_lhs,"dt_min")==0) {
      dtmin = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_min = %g\n", dtmin);
    } else if (strcmp(trim_lhs,"freq_type")==0) {
      trim_lhs = trim_rhs;
      for ( ; *trim_lhs; ++trim_lhs) *trim_lhs = tolower(*trim_lhs); // convert to all lowercase
      strcpy(freq_type,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set freq_type = %s\n", freq_type);
    } else if (strcmp(trim_lhs,"grid_units")==0) {
      trim_lhs = trim_rhs;
      for ( ; *trim_lhs; ++trim_lhs) *trim_lhs = tolower(*trim_lhs); // convert to all lowercase
      strcpy(grid_units,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set grid_units = %s\n", grid_units);
    } else if (strcmp(trim_lhs,"model_name")==0) {
      strcpy(model_name,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set model_name = %s\n", model_name);
    } else if (strcmp(trim_lhs,"restart_file")==0) {
      strcpy(restart_file,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_file = %s\n", restart_file);
    } else if (strcmp(trim_lhs,"restart_dir")==0) {
      strcpy(restart_dir,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set restart_dir = %s\n", restart_dir);
    } else if (strcmp(trim_lhs,"num_comp")==0) {
      ncomp = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set num_comp = %d\n", ncomp);
    } else if (strcmp(trim_lhs,"implicit_err_tol")==0) {
      implicit_err_tol = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set implicit_err_tol = %g\n", implicit_err_tol);
    } else if (strcmp(trim_lhs,"dT_tol")==0) {
      dT_tol = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dT_tol = %g\n", dT_tol);
    } else if (strcmp(trim_lhs,"du_tol")==0) {
      du_tol = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set du_tol = %g\n", du_tol);
    } else if (strcmp(trim_lhs,"dye_tol")==0) {
      dye_tol = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dye_tol = %g\n", dye_tol);
    } else if (strcmp(trim_lhs,"initial_dt")==0) {
      initial_dt = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set initial_dt = %g\n", initial_dt);
    } else if (strcmp(trim_lhs,"tmax")==0) {
      tmax = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set tmax = %g\n", tmax);
    } else if (strcmp(trim_lhs,"dump_hdf")==0) {
      dump_hdf = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dump_hdf = %d\n", dump_hdf);
    } else if (strcmp(trim_lhs,"dump_rad_vars")==0) {
      dump_rad_vars = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dump_rad_vars = %d\n", dump_rad_vars);
    } else if (strcmp(trim_lhs,"detect_tbounce")==0) {
      detect_tbounce = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set detect_tbounce = %d\n", detect_tbounce);
    } else if (strcmp(trim_lhs,"dt_dump")==0) {
      dt_dump = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_dump = %g\n", dt_dump);
    } else if (strcmp(trim_lhs,"dt_pdump")==0) {
      dt_pdump = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_pdump = %g\n", dt_pdump);
    } else if (strcmp(trim_lhs,"dt_restart")==0) {
      dt_restart = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set dt_restart = %g\n", dt_restart);
    } else if (strcmp(trim_lhs,"di_restart")==0) {
      di_restart = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set di_restart = %d\n", di_restart);
    } else if (strcmp(trim_lhs,"nstep_log")==0) {
      nstep_log = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set nstep_log = %d\n", nstep_log);
    } else if (strcmp(trim_lhs,"nstep_analysis")==0) {
      nstep_analysis = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set nstep_analysis = %d\n", nstep_analysis);
    } else if (strcmp(trim_lhs,"nstep_timing")==0) {
      nstep_timing = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set nstep_timing = %d\n", nstep_timing);
    } else if (strcmp(trim_lhs,"chat")==0) {
      chat = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set chat = %g\n", chat);
    } else if (strcmp(trim_lhs,"chat_safety")==0) {
      chat_safety = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set chat_safety = %g\n", chat_safety);
    } else if (strcmp(trim_lhs,"hydro_interp_order")==0) {
      hydro_interp_order = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set hydro_interp_order = %d\n", hydro_interp_order);
    } else if (strcmp(trim_lhs,"rad_interp_order")==0) {
      rad_interp_order = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set rad_interp_order = %d\n", rad_interp_order);
    } else if (strcmp(trim_lhs,"use_chat")==0) {
      use_chat = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set use_chat = %d\n", use_chat);
    } else if (strcmp(trim_lhs,"max_implicit_iter")==0) {
      max_implicit_iter = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set max_implicit_iter = %d\n", max_implicit_iter);
    } else if (strcmp(trim_lhs,"max_steps")==0) {
      max_steps = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set max_steps = %d\n", max_steps);
    } else if (strcmp(trim_lhs,"max_wtime")==0) {
      max_wtime = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set max_wtime = %g\n", max_wtime);
    } else if (strcmp(trim_lhs,"include_inelastic")==0) {
      include_inelastic = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set include_inelastic = %d\n", include_inelastic);
    } else if (strcmp(trim_lhs,"opac_param_file")==0) {
      strcpy(opac_param_file,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set opac_param_file = %s\n", opac_param_file);
    } else if (strcmp(trim_lhs,"opac_file")==0) {
      strcpy(opac_file,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set opac_file = %s\n", opac_file);
    } else if (strcmp(trim_lhs,"inelastic_root")==0) {
      strcpy(inelastic_root,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set inelastic_root = %s\n", inelastic_root);
    } else if (strcmp(trim_lhs,"eos_file")==0) {
      strcpy(eos_file,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set eos_file = %s\n", eos_file);
    } else if (strcmp(trim_lhs,"L_enu")==0) {
      L_enu = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set L_enu = %g\n", L_enu);
    } else if (strcmp(trim_lhs,"mass_inside_tracers")==0) {
      mass_inside_tracers = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set mass_inside_tracers = %g\n", mass_inside_tracers);
    } else if (strcmp(trim_lhs,"mass_outside_tracers")==0) {
      mass_outside_tracers = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set mass_outside_tracers = %g\n", mass_outside_tracers);
    } else if (strcmp(trim_lhs,"T_enu")==0) {
      T_enu = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set T_enu = %g\n", T_enu);
    } else if (strcmp(trim_lhs,"e_floor")==0) {
      e_floor = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set e_floor = %g\n", e_floor);
    } else if (strcmp(trim_lhs,"rho_floor")==0) {
      rho_floor = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set rho_floor = %g\n", rho_floor);
    } else if (strcmp(trim_lhs, "perturb_r1m")==0) {
      perturb_r1m = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set perturb_r1m = %g\n", perturb_r1m);
    } else if (strcmp(trim_lhs, "perturb_r1p")==0) {
      perturb_r1p = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set perturb_r1p = %g\n", perturb_r1p);
    } else if (strcmp(trim_lhs, "perturb_l1")==0) {
      perturb_l1 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_l1 = %d\n", perturb_l1);
    } else if (strcmp(trim_lhs, "perturb_m1")==0) {
      perturb_m1 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_m1 = %d\n", perturb_m1);
    } else if (strcmp(trim_lhs,"perturb_n1")==0) {
      perturb_n1 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_n1 = %d\n", perturb_n1);
    } else if (strcmp(trim_lhs,"perturb_dv1")==0) {
      perturb_dv1 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_dv1 = %g\n", perturb_dv1);
    } else if (strcmp(trim_lhs, "perturb_r2m")==0) {
      perturb_r2m = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set perturb_r2m = %g\n", perturb_r2m);
    } else if (strcmp(trim_lhs, "perturb_r2p")==0) {
      perturb_r2p = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set perturb_r2p = %g\n", perturb_r2p);
    } else if (strcmp(trim_lhs, "perturb_l2")==0) {
      perturb_l2 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_l2 = %d\n", perturb_l2);
    } else if (strcmp(trim_lhs, "perturb_m2")==0) {
      perturb_m2 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_m2 = %d\n", perturb_m2);
    } else if (strcmp(trim_lhs,"perturb_n2")==0) {
      perturb_n2 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_n2 = %d\n", perturb_n2);
    } else if (strcmp(trim_lhs,"perturb_dv2")==0) {
      perturb_dv2 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_dv2 = %g\n", perturb_dv2);
    } else if (strcmp(trim_lhs, "perturb_r3m")==0) {
      perturb_r3m = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set perturb_r3m = %g\n", perturb_r3m);
    } else if (strcmp(trim_lhs, "perturb_r3p")==0) {
      perturb_r3p = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set perturb_r3p = %g\n", perturb_r3p);
    } else if (strcmp(trim_lhs, "perturb_l3")==0) {
      perturb_l3 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_l3 = %d\n", perturb_l3);
    } else if (strcmp(trim_lhs, "perturb_m3")==0) {
      perturb_m3 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_m3 = %d\n", perturb_m3);
    } else if (strcmp(trim_lhs,"perturb_n3")==0) {
      perturb_n3 = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_n3 = %d\n", perturb_n3);
    } else if (strcmp(trim_lhs,"perturb_dv3")==0) {
      perturb_dv3 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_dv3 = %g\n", perturb_dv3);
    } else if (strcmp(trim_lhs,"perturb_level")==0) {
      perturb_level = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_level = %g\n", perturb_level);
    } else if (strcmp(trim_lhs,"perturb_delay")==0) {
      perturb_delay = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set perturb_delay = %g\n", perturb_delay);
    } else if (strcmp(trim_lhs,"rotate_Omega0")==0) {
      rotate_Omega0 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set rotate_Omega0 = %g\n", rotate_Omega0);
    } else if (strcmp(trim_lhs,"rotate_A")==0) {
      rotate_A = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set rotate_A = %g\n", rotate_A);
    } else if (strcmp(trim_lhs,"outer_radius")==0) {
      outer_radius = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set outer_radius = %g\n", outer_radius);
    } else if (strcmp(trim_lhs,"decomp_only")==0) {
      decomp_only = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set decomp_only = %d\n", decomp_only);
    } else if (strcmp(trim_lhs,"decomp_npmin")==0) {
      decomp_npmin = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set decomp_npmin = %d\n", decomp_npmin);
    } else if (strcmp(trim_lhs,"decomp_npmax")==0) {
      decomp_npmax = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set decomp_npmax = %d\n", decomp_npmax);
    } else if (strcmp(trim_lhs,"decomp_npskip")==0) {
      decomp_npskip = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set decomp_npskip = %d\n", decomp_npskip);
    } else if (strcmp(trim_lhs,"decomp_from_file")==0) {
      decomp_from_file = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set decomp_from_file = %d\n", decomp_from_file);
    } else if (strcmp(trim_lhs,"decomp_path")==0) {
      strcpy(decomp_path,trim_rhs);
      if (mpi_io_proc()) fprintf(stderr,"Set decomp_path = %s\n", decomp_path);
    // tianshu
    } else if (strcmp(trim_lhs,"use_kom")==0) {
      use_kom = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set use_kom = %d\n", use_kom);
    } else if (strcmp(trim_lhs,"kom_dtmin")==0) {
      kom_dtmin = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set kom_dtmin = %g\n", kom_dtmin);
    } else if (strcmp(trim_lhs,"kom_epsilon")==0) {
      kom_epsilon = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set kom_epsilon = %g\n", kom_epsilon);
    } else if (strcmp(trim_lhs,"kom_delta")==0) {
      kom_delta = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set kom_delta = %g\n", kom_delta);
    } else if (strcmp(trim_lhs,"Ecut1")==0) {
      Ecut1 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set Ecut1 = %g\n", Ecut1);
    } else if (strcmp(trim_lhs,"Ecut2")==0) {
      Ecut2 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set Ecut2 = %g\n", Ecut2);
    } else if (strcmp(trim_lhs,"Ecut3")==0) {
      Ecut3 = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set Ecut3 = %g\n", Ecut3);
    } else if (strcmp(trim_lhs,"rhocut")==0) {
      rhocut = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set rhocut = %g\n", rhocut);
#if (USE_LARGER_STEP==TRUE)
    } else if (strcmp(trim_lhs,"Nskip")==0) {
      Nskip = atoi(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set Nskip = %d\n", Nskip);
#endif
#if (RCOORD==SINH_MODIFIED)
    } else if (strcmp(trim_lhs,"facta")==0) {
      rx_info.facta = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set facta = %g\n", rx_info.facta);
    } else if (strcmp(trim_lhs,"factb")==0) {
      rx_info.factb = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set factb = %g\n", rx_info.factb);
    } else if (strcmp(trim_lhs,"factc")==0) {
      rx_info.factc = atof(trim_rhs);
      if (mpi_io_proc()) fprintf(stderr, "Set factc = %g\n", rx_info.factc);
#endif
    // end tianshu
    } else {
      if (mpi_io_proc()) fprintf(stderr,"Warning: ignoring unknown parameter %s\n", trim_lhs);
    }
  }

  mpi_fclose(fp);

  #if (NDIM==1)
  n2 = n3 = 1;
  #elif (NDIM==2)
  n3 = 1;
  #endif

  nhydro = NPRIM + ncomp;
  nvars = nhydro;
  #if (DO_RADIATION==TRUE)
  #if (PHOTON==TRUE)
  nvars += nr1*(1+NDIM);
  #elif (NEUTRINO==TRUE)
  nvars += (nr1+nr2+nr3)*(1+NDIM);
  #endif
  #endif
  ninterp = nvars + 2;    // interpolate all hydro/radiation primitives + pressure and Gamma

  free(line);
  free(lhs);
  free(rhs);

  if (dump_hdf==0 && dump_rad_vars==1) {
    if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  dump_hdf=0, so resetting dump_rad_vars=0\n");
    fflush(stderr);
    dump_rad_vars=0;
  }

  #if (ENFORCE_FLOORS==TRUE)
  if (e_floor < 0.0) {
    if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  e_floor has not been set!\n");
    fflush(stderr);
    exit(1);
  }
  if (rho_floor < 0.0) {
    if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  rho_floor has not been set!\n");
    fflush(stderr);
    exit(1);
  }
  #endif

  if (restart_create_1d) {
    #if (NDIM>1)
    if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_create_1d=1, but NDIM>1!\n");
    fflush(stderr);
    exit(1);
    #endif
    if (dump_hdf) {
      // MAS:  IS THIS TOTALLY NECESSARY?  WHY DOES IT NOT WORK TO HAVE HDF *AND* 1D RESTART DUMPS??
       // if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_create_1d=1, so resetting dump_hdf=0\n");
       // fflush(stderr);
       // dump_hdf=0;
    }
    if (restart_from_1d) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_create_1d=1, so resetting restart_from_1d=0\n");
      fflush(stderr);
      restart_from_1d=0;
    }
    if (restart_from_last) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_create_1d=1, so resetting restart_from_last=0\n");
      fflush(stderr);
      restart_from_last=0;
    }
    if (restart_from_hdf) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_create_1d=1, so resetting restart_from_hdf=0\n");
      fflush(stderr);
      restart_from_hdf=0;
    }
    if (!detect_tbounce) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_create_1d=1, so resetting detect_tbounce=1\n");
      fflush(stderr);
      detect_tbounce=1;
    }
  }

  if (restart_from_1d) {
    // #if (NDIM==1)
    // if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_from_1d=1, but NDIM==1!\n");
    // fflush(stderr);
    // exit(1);
    // #endif
    if (restart_from_hdf) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_from_1d=1, so resetting restart_from_hdf=0\n");
      fflush(stderr);
      restart_from_hdf=0;
    }
    if (restart_from_last) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_from_1d=1, so resetting restart_from_last=0\n");
      fflush(stderr);
      restart_from_last=0;
    }
    if (detect_tbounce) {
      if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  restart_from_1d=1, so resetting detect_tbounce=0\n");
      fflush(stderr);
      detect_tbounce=0;
    }
  }

  if (decomp_only) {
    #if (NDIM==1)
    if (mpi_io_proc()) fprintf(stderr,"[parse_input]:  decomp_only=1, but NDIM==1!\n");
    fflush(stderr);
    exit(1);
    #endif
  }

  TIMER_STOP;
  return;
}
