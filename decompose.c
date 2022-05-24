#include <inttypes.h>
#include <stdarg.h>

#include "decs.h"

// Minimum size of grid in each direction
// If using dendritic grid, must be >= 2^NG
// If not using dendritic grid, must be >= NG
#define MIN_IBLOCK_SIZE  4
#define MIN_JBLOCK_SIZE  2
#define MIN_KBLOCK_SIZE  2
// Minimum number of zones in order to refine...
#define MIN_JREFINE_SIZE 2
#if POLAR_AVG==TRUE
#define MIN_KREFINE_SIZE 4
#else
#define MIN_KREFINE_SIZE 2
#endif // POLAR_AVG

// Set the level of debugging output (1 for a little, 5 for a lot)
const int out_level=0;

void fnx_perr(const int level, const char *fmt, ...)
{  
  va_list ap;
  FILE* fp=stderr;
  
  if (mpi_io_proc() && level <= out_level) {
    va_start(ap, fmt);
    vfprintf(fp, fmt, ap);
    fflush(fp);
    va_end(ap);
  }
}

void fnx_pout(const int level, const char *fmt, ...)
{  
  va_list ap;
  FILE* fp=stdout;

  if (mpi_io_proc() && level <= out_level) {
    va_start(ap, fmt);
    vfprintf(fp, fmt, ap);
    fflush(fp);
    va_end(ap);
  }
}

#if (NDIM>1)
typedef struct result_ {
  int nproc_dendritic;
  int nproc_cart;
  float ideal;
  uint64_t actual;
  int worst;
  uint32_t max_neigh;
  int cart_dims[NDIM];
} result;

static result res;

#if (NDIM==3 && GEOM==SPHERICAL)
void dump_grid_VTK(proc_info* proc)
{
  FILE* fp=NULL;
  int ind,ii,jj,kk,ip,j,k,k1,nkp,j0,njp,k0;
  int ND_PTR nodes=NULL;
  double x,y,z;
  const int imax = 608;

  my_grid_dims[0] = n1;
  my_grid_dims[1] = n2;
  my_grid_dims[2] = n3;

  nodes = dendritic_malloc_int();
  
  /* Compute nodes for decomposition */
  ind = 0;
  for (ii=0; ii<=imax; ii++) {
    for (jj=0; jj<=(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<=(n3/4)/dk[ii][jj]; kk++) {
        ND_ELEM(nodes,ii,jj,kk) = ind++;
      }
    }
  }

  fp = fopen("grid.vtk","w");
  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"Unstructured Grid Example\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
  fprintf(fp,"POINTS %d float\n",ind);
  ind = 0;
  for (ii=0; ii<=imax; ii++) {
    for (jj=0; jj<=(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<=(n3/4)/dk[ii][jj]; kk++) {
        x = r_of_x(rx_info,startx[0] + ii*dx[0]);
        // x /= 1.0e5;
        y = th_of_x(thx_info,startx[1] + jj*dj[ii]*dx[1]);
        z = startx[2] + kk*dk[ii][jj]*dx[2];
        fprintf(fp,"%f %f %f\n",x,y,z);
      }
    }
  }
  
  ind = 0;
  for (ii=0; ii<imax; ii++) {
    for (jj=0; jj<(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<(n3/4)/dk[ii][jj]; kk++) {
        ind++;
      }
    }
  }

  fprintf(fp,"CELLS %d %d\n",ind,ind*(8+1));

  for (ii=0; ii<imax; ii++) {
    for (jj=0; jj<(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<(n3/4)/dk[ii][jj]; kk++) {
        // corner
        fprintf(fp,"8 %d ",nodes[ii][jj][kk]);

        // neighbor in the k-direction at lower i
        fprintf(fp,"%d ",nodes[ii][jj][kk+1]);

        // neighbor(s) in the j-direction at lower i(check for k-refinement)
        k1 = (kk*dk[ii][jj])/dk[ii][jj+1];
        nkp = (dk[ii][jj+1] < dk[ii][jj]) ? dk[ii][jj]/dk[ii][jj+1] : 1;
        for (k=k1; k<k1+nkp+1; k+=nkp) {
          fprintf(fp,"%d ",nodes[ii][jj+1][k]);
        }

        // neighbor(s) in the j- and k-directions at upper i (check for k-refinement)
        j0 = (jj*dj[ii])/dj[ii+1];
        njp = (dj[ii+1] < dj[ii]) ? dj[ii]/dj[ii+1] : 1;
        for (j=j0; j<j0+njp+1; j+=njp) {
          k0 = (kk*dk[ii][jj])/dk[ii+1][j];
          nkp = (dk[ii+1][j] < dk[ii][jj]) ? dk[ii][jj]/dk[ii+1][j] : 1;
          for (k=k0; k<k0+nkp+1; k+=nkp) {        
            fprintf(fp,"%d ",nodes[ii+1][j][k]);
          }
        }
        fprintf(fp,"\n");
      }
    }
  }

  fprintf(fp,"CELL_TYPES %d\n",ind);

  for (ii=0; ii<imax; ii++) {
    for (jj=0; jj<(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<(n3/4)/dk[ii][jj]; kk++) {
        fprintf(fp,"11\n");
      }
    }
  }
  
  fprintf(fp,"CELL_DATA %d\n",ind);
  fprintf(fp,"SCALARS jrefine float\n");
  // fprintf(fp,"SCALARS krefine float\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (ii=0; ii<imax; ii++) {
    for (jj=0; jj<(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<(n3/4)/dk[ii][jj]; kk++) {
        fprintf(fp,"%d\n",dj[ii]);
        // fprintf(fp,"%d\n",dk[ii][jj]);
        // fprintf(fp,"%d\n",ND_ELEM(part,ii,jj,kk));
        // fprintf(fp,"%f\n",r_of_x(istart[0]+ii*dx[0]));
      }
    }
  }
    
  fprintf(fp,"SCALARS procid float\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for (ii=0; ii<imax; ii++) {
    for (jj=0; jj<(n2/2)/dj[ii]; jj++) {
      for (kk=0; kk<(n3/4)/dk[ii][jj]; kk++) {
        for (ip=0; ip<numprocs; ip++) {
          if (ii>=proc[ip].istart[0] && ii<proc[ip].istop[0] 
            && jj*dj[ii]>=proc[ip].istart[1] && jj*dj[ii]<proc[ip].istop[1] 
            && kk*dk[ii][jj]>=proc[ip].istart[2] && kk*dk[ii][jj]<proc[ip].istop[2]) {
            fprintf(fp,"%d\n",ip);
            break;
          } 
        }
      }
    }
  }  
  
  fclose(fp);
}
#endif /* NDIM==3 */

/* Find the integer factors (not necessarily prime) of a given number n, 
 * store factors in array f, return number of factors
 */
int find_factors(int n, int **f)
{
  int nf,i;

  nf = 2;
  *f = malloc(2*sizeof(int));
  (*f)[0] = 1;
  for (i=2; i<n/2+1; i++) {
    if (n%i == 0) {
      (*f)[nf-1] = i;
      nf++;
      *f = realloc(*f, nf*sizeof(int));
    }
  }
  (*f)[nf-1] = n;
  
  return nf;
}

/* Determine which pairs of factors in f multiply to give n, then
 * store the pairs in c and return nc, the number of combinations
 */
int comb_factors2(int n, int *f, int nf, int **c)
{
  int i,j,nc=0;

  *c = malloc(2*sizeof(int));

  for (i=0; i<nf; i++) {
    for (j=0; j<nf; j++) {
      if (f[i]*f[j] == n) {
        *c = realloc(*c, (nc+1)*2*sizeof(int));
        (*c)[2*nc  ] = f[i];
        (*c)[2*nc+1] = f[j];
        nc++;
      }
    }
  }
  
  return nc;
}

/* Determine which triples of factors in f multiply to give n, then
 * store the triples in c and return nc, the number of combinations
 */
int comb_factors3(int n, int *f, int nf, int **c)
{
  int i,j,k,nc=0;

  *c = malloc(3*sizeof(int));

  for (i=0; i<nf; i++) {
    for (j=0; j<nf; j++) {
      for (k=0; k<nf; k++) {
        if (f[i]*f[j]*f[k] == n) {
          *c = realloc(*c, (nc+1)*3*sizeof(int));
          (*c)[3*nc  ] = f[i];
          (*c)[3*nc+1] = f[j];
          (*c)[3*nc+2] = f[k];
          nc++;
        }
      }
    }
  }
  
  return nc;
}

// Compute a measure of how far from the target trial is
float split_obj(int trial[], float target[])
{
  int i;
  float obj = 0;
  for (i=0; i<NDIM; i++)  obj += SQR(trial[i]-target[i]);
  
  return obj;
}

// Split the Cartesian portion of the grid
double split_cart(int np, float np1, float np2, float np3, int best[], int istop_dendritic, int jstop_dendritic, int jstart_dendritic)
{
  int *factors=NULL;
  int *com=NULL;
  int i,j,nf,nc,iobj_min,ok;
  int narr[3];
  float targ[3],obj_min,obj;

  nf = find_factors(np, &factors);

  #if (NDIM==3)
  nc = comb_factors3(np, factors, nf, &com);
  #else
  nc = comb_factors2(np, factors, nf, &com);
  #endif

  targ[0] = np1;
  targ[1] = np2;
  targ[2] = np3;
  narr[0] = n1 - istop_dendritic;
  narr[1] = jstart_dendritic - jstop_dendritic;
  narr[2] = n3;

  if (nc == 0) fnx_perr(0,"nc == 0\n");

  // Find the feasible combination of factors that minimizes the objective function
  obj_min = 1.0e20;
  iobj_min = -1;
  for (i=0; i<nc; i++) {
    /* For each combination of factors check that 1) the combination gives a block of
     * at least the minimum size, and 2) for a 3D dendritic grid, that the third factor
     * both divides the number of k-cells and that the quotient is a factor of 4
     */
    ok  = (narr[0]/com[NDIM*i  ] >= MIN_IBLOCK_SIZE);
    ok *= (narr[1]/com[NDIM*i+1] >= MIN_JBLOCK_SIZE);
    ok *= (narr[1] % com[NDIM*i+1] == 0);
    ok *= (narr[1]/com[NDIM*i+1] % MIN_JBLOCK_SIZE == 0);
    #if (NDIM==3)
    ok *= (narr[2]/com[NDIM*i+2] >= MIN_KBLOCK_SIZE);
    #if (GEOM==SPHERICAL)
    ok *= (narr[2] % com[NDIM*i+2] == 0);
    ok *= ((narr[2]/com[NDIM*i+2]) % 4 == 0);
    ok *= (np + 2*com[NDIM*i] + 1 <= numprocs);
    // ok *= (com[NDIM*i] <= (numprocs - np - 1)/2);
    #endif
    #endif
        
    if (!ok) continue;
    // Measure distance to ideal target
    obj = split_obj(&com[NDIM*i], targ);
    if (obj < obj_min) {
      // New minimizer found
      obj_min = obj;
      iobj_min = i;
    }
  }

  // If no feasible combinations found (e.g., np is prime), return trivial splitting
  if (iobj_min < 0) {
    for (i=0; i<NDIM; i++) best[i] = 1;
    if (narr[0] > narr[1] && narr[0] > narr[2]) best[0] = np;
    else if (narr[1] > narr[2]) best[1] = np;
    else best[2] = np;
    free(factors);
    free(com);
    
    return obj_min;
  }

  // Otherwise, return the minimizer
  for (i=0; i<NDIM; i++) best[i] = com[NDIM*iobj_min+i];
  free(factors);
  free(com);

  return obj_min;
}

// Check factors for feasibility
int check_factors(int np, int* factors, int nf, int istop_dendritic, int jstop_dendritic, int jstart_dendritic)
{
  int *com=NULL;
  int nc,check,ok,i,j;
  int narr[NDIM];

  if (np !=2 && nf==2) return 0;

  narr[0] = n1 - istop_dendritic;
  narr[1] = jstart_dendritic - jstop_dendritic;

  // Determine the specific combinations of factors that give np
  #if (NDIM==2)
  nc = comb_factors2(np, factors, nf, &com);
  #elif (NDIM==3)
  nc = comb_factors3(np, factors, nf, &com);
  narr[2] = n3;
  #endif

  // For each combination of factors check that 1) the combination gives a block of
  // at least the minimum size, and 2) for a 3D dendritic grid, that the last factor
  // both divides the number of k-cells and that the quotient is a multiple of MIN_KBLOCK_SIZE
  check = 0;
  for (i=0; i<nc; i++) {
    ok  = (narr[0]/com[NDIM*i  ] >= MIN_IBLOCK_SIZE);
    ok *= (narr[1]/com[NDIM*i+1] >= MIN_JBLOCK_SIZE);
    ok *= (narr[1] % com[NDIM*i+1] == 0);
    ok *= (narr[1]/com[NDIM*i+1] % MIN_JBLOCK_SIZE == 0);
    #if (NDIM==3)
    ok *= (narr[2]/com[NDIM*i+2] >= MIN_KBLOCK_SIZE);
    #if (GEOM==SPHERICAL)
    ok *= (narr[2] % com[NDIM*i+2] == 0);
    ok *= (narr[2]/com[NDIM*i+2] % MIN_KBLOCK_SIZE == 0);
    // ok *= (np + 2*com[NDIM*i] + 1 <= numprocs);
    // ok *= (com[NDIM*i] <= (numprocs - np - 1)/2);
    #endif
    #endif
        
    if (ok) check = 1;
  }

  free(com);

  return check;
}

// Copy (some) proc info from p1 to p2
void copy_proc(proc_info *p1, proc_info *p2)
{
  int i;

  for (i=0; i<NDIM; i++) {
    p2->istart[i] = p1->istart[i];
    p2->istop[i]  = p1->istop[i];
  }
  p2->ncells = p1->ncells;
}

// Try to split the biggest proc transversely from np0 to np1 into two parts,
// create new proc info, recompute cell counts, and return 1 if successful
int split_proc(proc_info *proc, int np0, int np1)
{
  uint64_t ncnt1,ncnt2,biggest;
  int ibiggest,ip,idim,i,j,s,jstop;
  int i0,i1,j0,j1,k00,k10,k01,k11,ok1,ok2;

  // Find the proc with the most cells (provided it has at the required minimum j- or k-cells)
  ibiggest = -1;
  biggest = 0;
  for (ip=np0; ip<np1; ip++) {
    i0  = proc[ip].istart[0];
    i1  = proc[ip].istop[0];
    j0  = proc[ip].istart[1]/dj[i0];
    j1  = proc[ip].istop[1]/dj[i0];
    ok1 = (j1-j0 >= 2*MIN_JREFINE_SIZE);
    fnx_pout(5,"[split_proc]:  ip=%d, ncells=%d\n",ip,proc[ip].ncells);
    fnx_pout(5,"[split_proc]:  i0=%d, i1=%d, j0=%d, j1=%d, dj=%d\n",i0,i1,j0,j1,dj[i0]);
    #if (NDIM==3)
    k00 = proc[ip].istart[2]/dk[i0][j0];
    k10 = proc[ip].istop[2]/dk[i0][j0];
    k01 = proc[ip].istart[2]/dk[i0][j1-1];
    k11 = proc[ip].istop[2]/dk[i0][j1-1];
    ok2 = (MIN(k10-k00,k11-k01) >= 2*MIN_KREFINE_SIZE);
    fnx_pout(5,"[split_proc]:  k00=%d, k10=%d, k01=%d, k11=%d\n",k00,k10,k01,k11);
    #else
    ok2 = 0;
    #endif
    
    if ((ok1 || ok2) && proc[ip].ncells > biggest) {
      biggest  = proc[ip].ncells;
      ibiggest = ip;
    }
  }
  fnx_pout(4,"[split_proc]:  Done with loop from %d to %d\n",np0,np1);

  if (ibiggest == -1) {   // could not split proc
    return 0;
  }

  // Loop over all procs above the best one and copy proc info forward from ip-1 to ip
  for (ip=np1; ip>ibiggest; ip--)  copy_proc(&proc[ip-1], &proc[ip]);

  #if (NDIM==3)
  // If 3-d, determine which (transverse) dimension to split along
  idim = -1;
  i0  = proc[ibiggest].istart[0];
  i1  = proc[ibiggest].istop[0];
  j0  = proc[ibiggest].istart[1]/dj[i0];
  j1  = proc[ibiggest].istop[1]/dj[i0];
  if (j1-j0 >= 2*MIN_JREFINE_SIZE) idim = 1;
  fnx_pout(5,"[split_proc]:  Determine split, ibiggest=%d, i0=%d, i1=%d, j0=%d, j1=%d\n",ibiggest,i0,i1,j0,j1);
  k00 = proc[ibiggest].istart[2]/dk[i0][j0];
  k10 = proc[ibiggest].istop[2]/dk[i0][j0];
  k01 = proc[ibiggest].istart[2]/dk[i0][j1-1];
  k11 = proc[ibiggest].istop[2]/dk[i0][j1-1];
  fnx_pout(5,"[split_proc]:  Determine split, k00=%d, k10=%d, k01=%d, k11=%d\n",k00,k10,k01,k11);
  if (MIN(k10-k00,k11-k01) >= 2*MIN_KREFINE_SIZE) {
    if (idim<0 || (MIN(k10-k00,k11-k01) > j1-j0)) idim = 2;
  }
  fnx_pout(4,"[split_proc]:  Splitting proc %d of size %d x %d x %d with %llu cells in dimension %d!\n",
    ibiggest,i1-i0,j1-j0,MIN(k10-k00,k11-k01),proc[ibiggest].ncells,idim);
  #else
  // Otherwise, can only split in theta-direction
  idim = 1;
  #endif

  // Split the proc by bisecting in the idim direction, rounding integer division up in the southern hemisphere
  s = (proc[ibiggest].istart[idim] + proc[ibiggest].istop[idim])/2;
  if (s>=n2/2) s = (proc[ibiggest].istart[idim] + proc[ibiggest].istop[idim] + 1)/2;
  fnx_pout(4,"[split_proc]:  istart=%d, istop=%d, s=%d\n",proc[ibiggest].istart[idim],proc[ibiggest].istop[idim],s);
  proc[ibiggest  ].istop [idim] = s;
  proc[ibiggest+1].istart[idim] = s;

  // Recompute the number of cells in each part
  ncnt1 = 0;
  ncnt2 = 0;
  i0 = proc[ibiggest].istart[0];
  i1 = proc[ibiggest].istop[0];
  for (i=i0; i<i1; i++) {
    j0 = proc[ibiggest].istart[1]/dj[i];
    j1 = proc[ibiggest].istop[1]/dj[i];
    fnx_pout(5,"[split_proc]:  recompute ncnt1, i=%d, j0=%d, j1=%d\n",i,j0,j1);
    #if (NDIM==3)
    fnx_pout(5,"[split_proc]:  recompute ncnt1, dk0=%d, dk1=%d, k00=%d, k10=%d, k01=%d, k11=%d\n",dk[i][j0],dk[i][j1-1],
      proc[ibiggest].istart[2]/dk[i][j0  ],proc[ibiggest].istop[2]/dk[i][j0  ],
      proc[ibiggest].istart[2]/dk[i][j1-1],proc[ibiggest].istop[2]/dk[i][j1-1]);
    for (j=j0; j<j1; j++) ncnt1 += (proc[ibiggest].istop[2] - proc[ibiggest].istart[2])/dk[i][j];
    #else
    ncnt1 += (j1-j0);
    #endif

    j0 = proc[ibiggest+1].istart[1]/dj[i];
    j1 = proc[ibiggest+1].istop[1]/dj[i];
    fnx_pout(5,"[split_proc]:  recompute ncnt2, i=%d, j0=%d, j1=%d\n",i,j0,j1);
    #if (NDIM==3)
    fnx_pout(5,"[split_proc]:  recompute ncnt2, dk0=%d, dk1=%d, k00=%d, k10=%d, k01=%d, k11=%d\n",dk[i][j0],dk[i][j1-1],
      proc[ibiggest+1].istart[2]/dk[i][j0  ],proc[ibiggest+1].istop[2]/dk[i][j0  ],
      proc[ibiggest+1].istart[2]/dk[i][j1-1],proc[ibiggest+1].istop[2]/dk[i][j1-1]);
    for (j=j0; j<j1; j++) ncnt2 += (proc[ibiggest+1].istop[2] - proc[ibiggest+1].istart[2])/dk[i][j];
    #else
    ncnt2 += (j1-j0);
    #endif
  }
  if (ncnt1 + ncnt2 != biggest) {
    fnx_perr(0,"[split_proc]:  Error!  ncells changed from %llu to %llu + %llu = %llu\n",biggest,ncnt1,ncnt2,ncnt1+ncnt2);
    #if (NDIM==3)
    fnx_perr(0,"[split_proc]:  proc %d had (%d, %d) x (%d, %d) x (%d, %d)\n",ibiggest,
      proc[ibiggest].istart[0],proc[ibiggest].istop[0],
      proc[ibiggest].istart[1],proc[ibiggest].istop[1],
      proc[ibiggest].istart[2],proc[ibiggest].istop[2]);
    #else
    fnx_perr(0,"[split_proc]:  proc %d had (%d, %d) x (%d, %d)\n",ibiggest,
      proc[ibiggest].istart[0],proc[ibiggest].istop[0],
      proc[ibiggest].istart[1],proc[ibiggest].istop[1]);
    #endif
    exit(1);
  }
  proc[ibiggest  ].ncells = ncnt1;
  proc[ibiggest+1].ncells = ncnt2;
  fnx_pout(4,"[split_proc]:  Done splitting!  proc %d has %llu cells, proc %d has %llu cells\n",
    ibiggest,ncnt1,ibiggest+1,ncnt2);

  return 1;
}

// Take the largest chunk of indices in i and shrink it, adjusting all procs above it in i
// If we can't shrink any chunks, no choice but to remove the last one
int shrink_biggest(proc_info *proc, int npmin, int npmax)
{
  int ip,j,ibig,ishell;
  uint64_t big,nshell;

  // First, find the proc with the biggest total number of cells (and at least the required minimum)
  ibig = -1;
  big = 0;
  for (ip=npmin; ip<npmax; ip++) {
    if ((proc[ip].istop[0]-proc[ip].istart[0] > MIN_IBLOCK_SIZE) && (proc[ip].ncells > big)) {
      big  = proc[ip].ncells;
      ibig = ip;
    }
  }
  if (ibig < 0) {
    fnx_pout(1,"[shrink_biggest]:  Could find no chunk to shrink in range %d to %d...\n",npmin,npmax);
    return 1;
  }

  // Now loop over procs above (?) this one and adjust istart[0], istop[0], and ncells
  for (ip=ibig; ip<npmax; ip++) {
    ishell = proc[ip].istop[0]-1;
    #if (NDIM==2)
    nshell = n2/dj[ishell];
    #elif (NDIM==3)
    nshell = 0;
    for (j=0; j<n2/dj[ishell]; j++)  nshell += n3/dk[ishell][j];
    #endif
    proc[ip].istop[0]--;
    proc[ip].ncells -= nshell;
    if (ip+1 < npmax) {
      proc[ip+1].istart[0]--;
      proc[ip+1].ncells += nshell;
    }
  }
  return 0;
}

// Take the smallest chunk of indices in i and grow it, adjusting all procs above it in i
int grow_smallest(proc_info *proc, int npmin, int npmax)
{
  int ip,j,ismall,ishell;
  uint64_t small,nshell;

  // First, find the proc with the smallest number of cells on the next i-shell above it
  ismall = -1;
  small = n1*n2*n3;
  for (ip=npmin; ip<npmax; ip++) {
    ishell = proc[ip].istop[0];
    #if (NDIM==2)
    nshell = n2/dj[ishell];
    #elif (NDIM==3)
    nshell = 0;
    for (j=0; j<n2/dj[ishell]; j++)  nshell += n3/dk[ishell][j];
    #endif
    if (proc[ip].ncells + nshell < small) {
      small  = proc[ip].ncells + nshell;
      ismall = ip;
    }
  }
  if (ismall < 0) {
    fnx_pout(1,"[grow_smallest]:  Could find no chunk to grow in range %d to %d...\n",npmin,npmax);
    return 1;
  }

  // Now loop over procs above (?) this one and adjust istart[0], istop[0], and ncells
  for (ip=ismall; ip<npmax; ip++) {
    ishell = proc[ip].istop[0];
    #if (NDIM==2)
    nshell = n2/dj[ishell];
    #elif (NDIM==3)
    nshell = 0;
    for (j=0; j<n2/dj[ishell]; j++)  nshell += n3/dk[ishell][j];
    #endif
    proc[ip].istop[0]++;
    proc[ip].ncells += nshell;
    if (ip+1 < npmax) {
      proc[ip+1].istart[0]++;
      proc[ip+1].ncells -= nshell;
    }
  }
  return 0;
}

// Add proc jp as a neighbor of proc ip in direction idir; reallocate arrays as needed
void link_neighbor(proc_info *proc, int ip, int jp, int idir)
{
  if (proc[ip].nneigh == 0) {
    fnx_pout(4,"[link_neighbor]:  proc %d just got its first neighbor %d in dir %d\n",ip,jp,idir);
    proc[ip].nneigh++;
    proc[ip].neigh_id = malloc(sizeof(uint32_t));
    proc[ip].neigh_dir = malloc(sizeof(int));
    proc[ip].neigh_id[0] = jp;
    proc[ip].neigh_dir[0] = idir;
  } else {
    fnx_pout(4,"[link_neighbor]:  proc %d has a new neighbor %d in dir %d\n",ip,jp,idir);
    proc[ip].nneigh++;
    proc[ip].neigh_id = realloc(proc[ip].neigh_id, proc[ip].nneigh*sizeof(uint32_t));
    proc[ip].neigh_dir = realloc(proc[ip].neigh_dir, proc[ip].nneigh*sizeof(int));
    proc[ip].neigh_id[proc[ip].nneigh-1] = jp;
    proc[ip].neigh_dir[proc[ip].nneigh-1] = idir;
  }
}

void is_neighbor(proc_info *proc, int ip, int jp)
{
  int idir,d,overlap;
  int cheap_check = 0;

  if (ip==jp) fnx_perr(0,"[is_neighbor]:  ip==jp=%d!\n",ip);

  // Determine the direction in which procs ip and jp are neighbors
  idir = 0;
  for (d=0; d<NDIM; d++) {
    overlap = ((proc[ip].istop[d] - proc[jp].istart[d])*(proc[ip].istart[d] - proc[jp].istop[d]) < 0);
    if (!overlap) idir = d;
    cheap_check += overlap;
  }

  // If procs ip and jp do not share NDIM-1 overlapping dimensions, then they are not neighbors
  if (cheap_check != NDIM-1) return;

  // If procs ip and jp are neighbors, link their proc info to each other
  if (periodic[idir]) {    // treat this carefully
    if (proc[ip].istop[idir] == proc[jp].istart[idir]) {
      link_neighbor(proc,ip,jp,  idir+1);
      link_neighbor(proc,jp,ip,-(idir+1));
    }
    if (proc[jp].istop[idir] == proc[ip].istart[idir]) {
      link_neighbor(proc,ip,jp,-(idir+1));
      link_neighbor(proc,jp,ip,  idir+1);
    }
    if (proc[ip].istart[idir] == 0 && proc[jp].istop[idir] == global_grid_dims[idir]) {
      link_neighbor(proc,ip,jp,-(idir+1));
      link_neighbor(proc,jp,ip,  idir+1);
    }
    if (proc[jp].istart[idir] == 0 && proc[ip].istop[idir] == global_grid_dims[idir]) {
      link_neighbor(proc,ip,jp,  idir+1);
      link_neighbor(proc,jp,ip,-(idir+1));
    }
  } else {
    if (proc[ip].istop[idir] == proc[jp].istart[idir]) {
      link_neighbor(proc,ip,jp,  idir+1);
      link_neighbor(proc,jp,ip,-(idir+1));
    } else if (proc[jp].istop[idir] == proc[ip].istart[idir]) {
      link_neighbor(proc,ip,jp,-(idir+1));
      link_neighbor(proc,jp,ip,  idir+1);
    }
  }

  #if (NDIM==3 && GEOM==SPHERICAL)
  /* For 3-d spherical, if ip and jp overlap in phi ACROSS the north/south 
   * polar axis, then they are additionally neighbors in the -/+ theta-direction, 
   * respectively (idir=1, so use -/+ 2 for last argument)
   */
  if (idir!=0 && (proc[ip].istop[2] - (proc[jp].istart[2]+n3/2)%n3)*(proc[ip].istart[2] - (proc[jp].istop[2]+n3/2-1)%n3+1) < 0) {
    if (proc[ip].istart[0]==0 && proc[jp].istart[0]==0) {
      link_neighbor(proc,ip,jp,-1);
      link_neighbor(proc,jp,ip,-1);
      fnx_pout(0,"[is_neighbor]:  idir=%d, ip=%d [%d,%d]x[%d,%d]x[%d,%d] and jp=%d [%d,%d]x[%d,%d]x[%d,%d] linked in -1 direction\n",idir,
        ip,proc[ip].istart[0],proc[ip].istop[0],proc[ip].istart[1],proc[ip].istop[1],proc[ip].istart[2],proc[ip].istop[2],
        jp,proc[jp].istart[0],proc[jp].istop[0],proc[jp].istart[1],proc[jp].istop[1],proc[jp].istart[2],proc[jp].istop[2]);
    }
    if (proc[ip].istart[1]==0 && proc[jp].istart[1]==0) {
      link_neighbor(proc,ip,jp,-2);
      link_neighbor(proc,jp,ip,-2);
      fnx_pout(0,"[is_neighbor]:  idir=%d, ip=%d [%d,%d]x[%d,%d]x[%d,%d] and jp=%d [%d,%d]x[%d,%d]x[%d,%d] linked in -2 direction\n",idir,
        ip,proc[ip].istart[0],proc[ip].istop[0],proc[ip].istart[1],proc[ip].istop[1],proc[ip].istart[2],proc[ip].istop[2],
        jp,proc[jp].istart[0],proc[jp].istop[0],proc[jp].istart[1],proc[jp].istop[1],proc[jp].istart[2],proc[jp].istop[2]);
    }
    if (proc[ip].istop[1]==n2 && proc[jp].istop[1]==n2) {
      link_neighbor(proc,ip,jp, 2);
      link_neighbor(proc,jp,ip, 2);
      fnx_pout(0,"[is_neighbor]:  idir=%d, ip=%d [%d,%d]x[%d,%d]x[%d,%d] and jp=%d [%d,%d]x[%d,%d]x[%d,%d] linked in +2 direction\n",idir,
        ip,proc[ip].istart[0],proc[ip].istop[0],proc[ip].istart[1],proc[ip].istop[1],proc[ip].istart[2],proc[ip].istop[2],
        jp,proc[jp].istart[0],proc[jp].istop[0],proc[jp].istart[1],proc[jp].istop[1],proc[jp].istart[2],proc[jp].istop[2]);
    }
  }
  #endif

  return;
}

#if (NDIM>1)
int decompose(proc_info *proc)
{
  int i,j,k,ii,dd,tmp,ncart,np1,np2,np3,closest[3],nproc_dendritic,istart_cart;
  int istop_jrefine_local,istop_krefine_local,jstop_krefine_local,jstart_kcoarsen_local;
  uint64_t w[n1],ncells_cart_total;
  uint64_t ncells_active_below,ncells_active_above,ncells_active_total;
  uint64_t ncells_dendritic_below,ncells_dendritic_above,ncells_dendritic_total;
  int* factors=NULL;
  int nproc_dendritic_target,nproc_dendritic_max,nproc_cart_target,nproc_cart_min;
  int nproc_dendritic_total,nproc_dendritic_above,nproc_dendritic_below;
  int nproc_cart_best,nf,check;
  double obj_min,obj,nblock_target;

  int best_iskip,final,iskip_max,iskip,nproc_dendritic_rad,nproc_dendritic_rad_max,nassigned,ilast;
  uint64_t lowest_max,ncnt,next,ncheck;
  uint64_t ncells_dendritic_target_per_block,ncells_dendritic_per_iblock;
  int split_failed,iter,ip,jp,npmin=0,npmax;

  uint64_t mymax;
  int imymax;
  float dii,djj,dkk;
  int cproc,i0,i1,j0,j1,k00,k10,k01,k11;
  int imax_cnt,imax_neigh;
  uint64_t max_cnt,total_cnt,dendritic_cnt,cart_cnt,max_neigh;

  // If there is k-refinement at the center, must have at least MIN_IBLOCK_SIZE i-zones in dendritic region
  istop_krefine_local   = (istop_krefine         > 0 ) ? MAX(MAX(istop_jrefine,istop_krefine),MIN_IBLOCK_SIZE) : 0;
  istop_jrefine_local   = (istop_jrefine         > 0 ) ? MAX(MAX(istop_jrefine,istop_krefine_local),MIN_IBLOCK_SIZE) : 0;
  // If there is k-refinement at the poles, must have at least MIN_JBLOCK_SIZE j-zones in dendritic region
  jstop_krefine_local   = (jstop_krefine  [n1-1] > 0 ) ? MAX(jstop_krefine[n1-1],MIN_JBLOCK_SIZE)              : 0;
  jstart_kcoarsen_local = (jstart_kcoarsen[n1-1] < n2) ? MIN(jstart_kcoarsen[n1-1],n2-MIN_JBLOCK_SIZE)         : n2;
  fnx_pout(1,"[decompose]:  istop_jrefine_local   = %d\n",istop_jrefine_local);
  fnx_pout(1,"[decompose]:  istop_krefine_local   = %d\n",istop_krefine_local);
  fnx_pout(1,"[decompose]:  jstop_krefine_local   = %d\n",jstop_krefine_local);
  fnx_pout(1,"[decompose]:  jstart_kcoarsen_local = %d\n",jstart_kcoarsen_local);

  // Calculate the total number of active cells, Cartesian cells, and dendritic cells
  ncells_dendritic_total = 0;
  ncells_active_below = 0;
  for (i=0; i<istop_jrefine_local; i++) {
    #if (NDIM==3)
    for (j=0; j<n2/dj[i]; j++) ncells_active_below += n3/dk[i][j];
    #else
    ncells_active_below += n2/dj[i];
    #endif
  }
  ncells_active_above = 0;
  for (i=istop_jrefine_local; i<n1; i++) {
    #if (NDIM==3)
    for (j=0; j<n2/dj[i]; j++) ncells_active_above += n3/dk[i][j];
    #else
    ncells_active_above += n2/dj[i];
    #endif
  }
  ncells_active_total    = ncells_active_below + ncells_active_above;
  ncells_cart_total      = (n1-istop_jrefine_local)*(jstart_kcoarsen_local-jstop_krefine_local)*n3;
  ncells_dendritic_below = ncells_active_below;
  ncells_dendritic_above = ncells_active_above - ncells_cart_total;
  ncells_dendritic_total = ncells_dendritic_below + ncells_dendritic_above;
  fnx_pout(2,"[decompose]:  ncells_active_below                   = %llu\n",ncells_active_below);
  fnx_pout(2,"[decompose]:  ncells_active_above                   = %llu\n",ncells_active_above);
  fnx_pout(2,"[decompose]:  ncells_active_total                   = %llu\n",ncells_active_total);
  fnx_pout(2,"[decompose]:  ncells_cart_total                     = %llu\n",ncells_cart_total);
  fnx_pout(2,"[decompose]:  ncells_dendritic_below                = %llu\n",ncells_dendritic_below);
  fnx_pout(2,"[decompose]:  ncells_dendritic_above                = %llu\n",ncells_dendritic_above);
  fnx_pout(2,"[decompose]:  ncells_dendritic_total                = %llu\n",ncells_dendritic_total);

  // Target number of procs for dendritic region equals proportion of dendritic cells, but at least 1
  nproc_dendritic_target = (1.0*ncells_dendritic_total)/ncells_active_total * numprocs;
  if (istop_jrefine_local > 0) nproc_dendritic_target = MAX(nproc_dendritic_target,1);
  if (istop_krefine_local > 0) nproc_dendritic_target = MAX(nproc_dendritic_target,3);
  #endif
  // Target number of Cartesian cells is the compliment (but an even number or 1 at the min)
  nproc_cart_target  = numprocs - nproc_dendritic_target;
  nproc_cart_target -= nproc_cart_target % 2;  // floor to an even number
  nproc_cart_target  = MAX(nproc_cart_target,1);
  fnx_pout(1,"[decompose]:  numprocs                              = %d\n",numprocs);
  fnx_pout(1,"[decompose]:  nproc_dendritic_target                = %d\n",nproc_dendritic_target);
  fnx_pout(1,"[decompose]:  nproc_cart_target                     = %d\n",nproc_cart_target);

  // // Compute the minimum number of Cartesian procs
  // // Allow this number to be up to 10% below the target, but at least 2 below
  // // (and always at least 1 proc)
  // nproc_cart_min = MAX(MIN(0.9*nproc_cart_target,nproc_cart_target-2),1);
  nproc_cart_min = 1;
  // If no dendritic grid, then just use the target as the min
  if (ncells_dendritic_total == 0) nproc_cart_min = nproc_cart_target;

  /* Now, using blocks of minimum size (MIN_IBLOCK_SIZE,MIN_JBLOCK_SIZE,MIN_KBLOCK_SIZE), break
   * up the dendritic region into sub-blocks based on the refinement level at the lower corner
   * First, compute how many min-size blocks fit into the dendritic grid
   */
  nproc_dendritic_max = 0;
  for (i=0; i<n1; i+=MIN_IBLOCK_SIZE) {
    if (i<istop_jrefine_local) {
      #if (NDIM==3)
      for (j=0; j<n2/dj[i]/2; j+=MIN_JBLOCK_SIZE)
        nproc_dendritic_max += 2*n3/(dk[i][j]*MIN_KBLOCK_SIZE);
      #else
      nproc_dendritic_max += n2/(dj[i]*MIN_JBLOCK_SIZE);
      #endif
    }
    #if (NDIM==3)
    else {
      for (j=0; j<jstop_krefine_local/dj[i]; j+=MIN_JBLOCK_SIZE)
        nproc_dendritic_max += 2*n3/(dk[i][j]*MIN_KBLOCK_SIZE);
    }
    #endif
  }
  /* If the number of dendritic blocks is less than the compliment of the number of "Cartesian" blocks,
   * then increase the number of "Cartesian" blocks
   */
  if (numprocs - nproc_cart_min > nproc_dendritic_max) nproc_cart_min = numprocs - nproc_dendritic_max;
  fnx_pout(1,"[decompose]:  nproc_dendritic_max                   = %d\n",nproc_dendritic_max);
  fnx_pout(1,"[decompose]:  nproc_cart_min                        = %d\n",nproc_cart_min);


  /* ----------------------------------------------------------------------------------------------
   * I.  Split the Cartesian part of the domain:
   * Try different numbers of Cartesian processors from target down to min allowed
   * Check factor combinations for feasibility and find the minimizer of the objective function
   * This section does not actually assign any procs yet; that is done last
   */
  obj_min = 1.0e10;
  nproc_cart_best = -1;
  for (ncart=nproc_cart_target; ncart>=nproc_cart_min; ncart-=2) {
    // Find the integer factors of ncart
    nf = find_factors(ncart, &factors);
    // Check that the factors meet certain conditions (give min block size, etc.)
    check = check_factors(ncart, factors, nf, istop_krefine_local, jstop_krefine_local, jstart_kcoarsen_local);    
    
    if (check) {
      nblock_target = pow((1.0*ncells_cart_total)/ncart,1.0/NDIM);
      np1 = (n1-istop_jrefine_local)/nblock_target;
      np2 = (jstart_kcoarsen_local-jstop_krefine_local)/nblock_target;
      np3 = n3/nblock_target;
      np1 = MAX(np1,1);
      np2 = MAX(np2,1);
      np3 = MAX(np3,1);
      obj = split_cart(ncart,np1,np2,np3,closest,istop_jrefine_local,jstop_krefine_local,jstart_kcoarsen_local);
      fnx_pout(3,"[decompose]:  obj=%e, ncart=%d, np1=%d, np2=%d, np3=%d\n",obj,ncart,np1,np2,np3);
      
      if (obj < obj_min) {
        // New minimizer found
        if (nproc_cart_best == -1) {
          if (ncells_dendritic_total > 0) nproc_cart_min = MAX(MIN(0.98*ncart,ncart-2),1);
        }
        obj_min = obj;
        nproc_cart_best = ncart;
      }
    }
  }
  
  if (nproc_cart_best < 0) {
    #if (NDIM==3)
    // If no feasible combination found, default to the best trivial decomposition known
    // nproc_cart_best = (numprocs - 1)/3;
    nproc_cart_best = nproc_cart_target;
    nproc_cart_best = MIN(nproc_cart_best,(n1-istop_jrefine_local)/MIN_IBLOCK_SIZE);
    #else
    // If no feasible combination found, default to the target number of procs
    nproc_cart_best = nproc_cart_target;
    #endif
    np1 = nproc_cart_best;
    np2 = 1;
    np3 = 1;
    closest[0] = nproc_cart_best;
    closest[1] = 1;
    closest[2] = 1;
  }
  else {
    // Recompute the Cartesian splitting for the best processor count found
    nblock_target = pow((1.0*ncells_cart_total)/nproc_cart_best,1.0/NDIM);
    np1 = (n1-istop_jrefine_local)/nblock_target;
    np2 = (jstart_kcoarsen_local-jstop_krefine_local)/nblock_target;
    np3 = n3/nblock_target;
    np1 = MAX(np1,1);
    np2 = MAX(np2,1);
    np3 = MAX(np3,1);
    obj = split_cart(nproc_cart_best,np1,np2,np3,closest,istop_jrefine_local,jstop_krefine_local,jstart_kcoarsen_local);
  }
  if (closest[0]>(n1-istop_jrefine_local)) {
    fnx_perr(0,"[decompose]:  Splitting Cartesian grid failed!  closest[0]=%d, (n1-istop_jrefine_local)=%d\n",
      closest[0],(n1-istop_jrefine_local));
    return 0;
  }
  fnx_pout(1,"[decompose]:  nproc_cart_best                       = %d\n",nproc_cart_best);
  
  nproc_dendritic_total = numprocs - nproc_cart_best;
  if (nproc_dendritic_total > 0) {
    nproc_dendritic_above = (1.0*ncells_dendritic_above)/ncells_dendritic_total * nproc_dendritic_total;
    #if (NDIM==3)
    nproc_dendritic_above -= nproc_dendritic_above % 2;  // force this to be even
    nproc_dendritic_above = MAX(nproc_dendritic_above,2*closest[0]);  // force this to be at least 2*closest[0]
    #endif
    nproc_dendritic_below = nproc_dendritic_total - nproc_dendritic_above;
    nproc_dendritic_below = MAX(nproc_dendritic_below,1);
  } else {
    nproc_dendritic_above = 0;
    nproc_dendritic_below = 0;
  }
  #if (NDIM==3)
  fnx_pout(2,"[decompose]:  target Cartesian block size           = (%d, %d, %d)\n",np1,np2,np3);
  fnx_pout(1,"[decompose]:  closest Cartesian block size          = (%d, %d, %d)\n",closest[0],closest[1],closest[2]);
  #else
  fnx_pout(2,"[decompose]:  target Cartesian block size           = (%d, %d)\n",np1,np2);
  fnx_pout(1,"[decompose]:  closest Cartesian block size          = (%d, %d)\n",closest[0],closest[1]);
  #endif
  fnx_pout(1,"[decompose]:  nproc_dendritic_below                 = %d\n",nproc_dendritic_below);
  fnx_pout(1,"[decompose]:  nproc_dendritic_above                 = %d\n",nproc_dendritic_above);
  fnx_pout(1,"[decompose]:  nproc_dendritic_total                 = %d\n",nproc_dendritic_total);


  /* ----------------------------------------------------------------------------------------------
   * II.  Split the dendritic part of the grid above istop_krefine_local
   */
  if (istop_krefine_local > 0) {
    lowest_max = 2000000000;
    
    // First, do the left part (from 0 to jstop_krefine_local)
    npmin = 0;
    nassigned = 0;
    for (ip=0; ip<closest[0]; ip++) {
      i0 = istop_krefine_local +  ip   *(n1-istop_krefine_local)/closest[0];
      i1 = istop_krefine_local + (ip+1)*(n1-istop_krefine_local)/closest[0];
      i1 = MIN(i1,n1);
      fnx_pout(4,"[decompose]:  left above, ip=%d, i0=%d, i1=%d, dii=%d\n",ip,i0,i1,i1-i0);
      
      ncnt = 0;
      for (i=i0; i<i1; i++) {
        j0 = 0;
        j1 = jstop_krefine_local/dj[i];
        #if (NDIM==3)
        for (j=j0; j<j1; j++) ncnt += n3/dk[i][j];
        #else
        ncnt += (j1-j0);
        #endif
      }

      cproc = npmin + ip;
      proc[cproc].istart[0] = i0;
      proc[cproc].istop[0]  = i1;
      proc[cproc].istart[1] = 0;
      proc[cproc].istop[1]  = jstop_krefine_local;
      #if (NDIM==3)
      proc[cproc].istart[2] = 0;
      proc[cproc].istop[2]  = n3;
      #endif
      proc[cproc].ncells = ncnt;
      nassigned++;
    }

    split_failed = 0;
    iter = 0;
    while (nassigned < nproc_dendritic_above/2) {
      // Try to split the biggest proc transversely into two parts
      if (!split_proc(proc,npmin,npmin+nassigned)) {
        split_failed = 1;
        break;
      }
      fnx_pout(5,"[decompose]:  did a left above split, iter=%d\n",iter);        
      iter++;
      nassigned++;
    }
    fnx_pout(2,"\n[decompose]:  Done splitting left above!  nassigned=%d, should be=%d\n\n",nassigned,nproc_dendritic_above/2);

    for (ip=npmin; ip<npmin+nassigned; ip++) {
      #if (NDIM==3)
      fnx_pout(4,"[decompose]:  left above, ip=%d, ncells=%d, (%d,%d) x (%d,%d) x (%d,%d)\n",
        ip,proc[ip].ncells,
        proc[ip].istart[0],proc[ip].istop[0],
        proc[ip].istart[1],proc[ip].istop[1],
        proc[ip].istart[2],proc[ip].istop[2]);
      #else
      fnx_pout(4,"[decompose]:  left above, ip=%d, ncells=%d, (%d,%d) x (%d,%d)\n",
        ip,proc[ip].ncells,
        proc[ip].istart[0],proc[ip].istop[0],
        proc[ip].istart[1],proc[ip].istop[1]);
      #endif
    }

    if (split_failed || (nassigned != nproc_dendritic_above/2)) {
      fnx_perr(0,"[decompose]:  Splitting on left above failed!\n");
      nproc_dendritic_above = 2*nassigned;
      nproc_dendritic_below = nproc_dendritic_total - nproc_dendritic_above;
      // return 0;
      if (!decomp_only) exit(1);
    }

    // Second, do the right part
    npmin = nassigned;
    nassigned = 0;
    for (ip=0; ip<closest[0]; ip++) {
      i0 = istop_krefine_local +  ip   *(n1-istop_krefine_local)/closest[0];
      i1 = istop_krefine_local + (ip+1)*(n1-istop_krefine_local)/closest[0];
      i1 = MIN(i1,n1);
      fnx_pout(4,"[decompose]:  right above, ip=%d, i0=%d, i1=%d, dii=%d\n",ip,i0,i1,i1-i0);
      
      ncnt = 0;
      for (i=i0; i<i1; i++) {
        j0 = jstart_kcoarsen_local/dj[i];
        j1 = n2/dj[i];
        #if (NDIM==3)
        for (j=j0; j<j1; j++) ncnt += n3/dk[i][j];
        #else
        ncnt += (j1-j0);
        #endif
      }

      cproc = npmin + ip;
      proc[cproc].istart[0] = i0;
      proc[cproc].istop[0]  = i1;
      proc[cproc].istart[1] = jstart_kcoarsen_local;
      proc[cproc].istop[1]  = n2;
      #if (NDIM==3)
      proc[cproc].istart[2] = 0;
      proc[cproc].istop[2]  = n3;
      #endif
      proc[cproc].ncells = ncnt;
      nassigned++;
    }

    split_failed = 0;
    iter = 0;
    while (nassigned < nproc_dendritic_above/2) {
      // Try to split the biggest proc transversely into two parts
      if (!split_proc(proc,npmin,npmin+nassigned)) {
        split_failed = 1;
        break;
      }
      fnx_pout(5,"[decompose]:  did a right above split, iter=%d\n",iter);        
      iter++;
      nassigned++;
    }
    fnx_pout(2,"\n[decompose]:  Done splitting right above!  nassigned=%d, should be=%d\n\n",nassigned,nproc_dendritic_above/2);

    for (ip=npmin; ip<npmin+nassigned; ip++) {
      #if (NDIM==3)
      fnx_pout(4,"[decompose]:  right above, ip=%d, ncells=%d, (%d,%d) x (%d,%d) x (%d,%d)\n",
        ip,proc[ip].ncells,
        proc[ip].istart[0],proc[ip].istop[0],
        proc[ip].istart[1],proc[ip].istop[1],
        proc[ip].istart[2],proc[ip].istop[2]);
      #else
      fnx_pout(4,"[decompose]:  right above, ip=%d, ncells=%d, (%d,%d) x (%d,%d)\n",
        ip,proc[ip].ncells,
        proc[ip].istart[0],proc[ip].istop[0],
        proc[ip].istart[1],proc[ip].istop[1]);
      #endif
    }

    if (split_failed || (nassigned != nproc_dendritic_above/2)) {
      fnx_perr(0,"[decompose]:  Splitting on right above failed!\n");
      return 0;
      if (!decomp_only) exit(1);
    }

  } else {
    // No dendritic grid
    istart_cart = 0;
  }  // end if (istop_krefine_local > 0)
  fnx_pout(1,"[decompose]:  Done with dendritic grid above procs!\n");

  /* ----------------------------------------------------------------------------------------------
   * III.  Split the dendritic part of the grid below istop_jrefine_local
   */
  if (istop_jrefine_local > 0) {
    best_iskip = MIN_IBLOCK_SIZE;
    lowest_max = 2000000000;
    final = 0;
    iskip_max = istop_jrefine_local;
    // Try different iblock sizes
    for (iskip = MIN_IBLOCK_SIZE; iskip <= iskip_max; iskip++) {
      // For this iblock size, compute the number of radial divisions
      nproc_dendritic_rad_max = istop_jrefine_local/iskip;
      nproc_dendritic_rad_max = MAX(nproc_dendritic_rad_max,1);
      nproc_dendritic_rad     = MIN(nproc_dendritic_below,nproc_dendritic_rad_max);      
      fnx_pout(3,"[decompose]:  iskip                                 = %d\n",iskip);
      fnx_pout(3,"[decompose]:  nproc_dendritic_rad_max               = %d\n",nproc_dendritic_rad_max);
      fnx_pout(3,"[decompose]:  nproc_dendritic_rad                   = %d\n",nproc_dendritic_rad);

      ncells_dendritic_target_per_block = (1.0*ncells_dendritic_below)/nproc_dendritic_below;
      ncells_dendritic_per_iblock       = (1.0*ncells_dendritic_below)/nproc_dendritic_rad;
      fnx_pout(3,"[decompose]:  ncells_dendritic_below                = %llu\n",ncells_dendritic_below);
      fnx_pout(3,"[decompose]:  nproc_dendritic_below                 = %d\n",nproc_dendritic_below);
      fnx_pout(3,"[decompose]:  ncells_dendritic_target_per_block     = %llu\n",ncells_dendritic_target_per_block);
      fnx_pout(3,"[decompose]:  ncells_dendritic_per_iblock           = %llu\n",ncells_dendritic_per_iblock);
      npmin = nproc_dendritic_above;
      nassigned = 0;
      ilast = 0;
      ncheck = 0;
      while (nassigned < nproc_dendritic_rad) {
        ii = ilast + iskip;
        ii = MIN(ii,n1);
        fnx_pout(4,"[decompose]:  ii = %d, nassigned = %d\n",ii,nassigned);
        // ncnt is the number of zones in this i-block from ilast to ilast+iskip
        ncnt = 0;
        #if (NDIM==3)
        for (i=ilast; i<ii; i++) for (j=0; j<n2/dj[i]; j++) ncnt += n3/dk[i][j];
        #else
        for (i=ilast; i<ii; i++) ncnt += n2/dj[i];
        #endif
        // ncheck is the max number of i-block cells allowed to be added
        if (nassigned < nproc_dendritic_rad) {
          ncheck = ncells_dendritic_target_per_block;
        } else {
          ncheck = ncells_dendritic_per_iblock;
        }
        fnx_pout(4,"[decompose]:  ncheck = %llu, ncnt = %llu\n",ncheck,ncnt);
        // If another shell of cells can be added without going over, do it until limit reached
        while (ncnt < ncheck && ii < n1) {
          #if (NDIM==3)
          next = 0;
          for (j=0; j<n2/dj[ii]; j++)  next += n3/dk[ii][j];
          #else
          next = n2/dj[ii];
          #endif
          if (abs(ncnt-ncheck) < abs(next-ncheck)) {
            break;
          }
          ncnt += next;
          ii++;
          fnx_pout(5,"[decompose]:  added a shell, next = %llu, ncnt = %llu, ii = %d\n",next,ncnt,ii);
        }
        cproc = npmin + nassigned;
        fnx_pout(4,"[decompose]:  cproc=%d, ilast = %d, ii = %d, ncnt = %llu\n",cproc,ilast,ii,ncnt);
        proc[cproc].istart[0] = ilast;
        proc[cproc].istop[0]  = ii;
        proc[cproc].istart[1] = 0;
        proc[cproc].istop[1]  = n2;
        #if (NDIM==3)
        proc[cproc].istart[2] = 0;
        proc[cproc].istop[2]  = n3;
        #endif
        proc[cproc].ncells = ncnt;
        nassigned++;
        ilast = ii;
      } /* end while */

      // Shrink/grow if over/under the limit
      while (ilast > istop_jrefine_local && ilast >= MIN_IBLOCK_SIZE) {
        split_failed = shrink_biggest(proc,npmin,npmin+nassigned);
        ilast--;
        fnx_pout(5,"[decompose]:  Did shrink, ilast = %d\n",ilast);
      }
      fnx_pout(3,"[decompose]:  Done shrinking!\n");
      while (ilast < istop_jrefine_local) {
        split_failed = grow_smallest(proc,npmin,npmin+nassigned);
        ilast++;
        fnx_pout(5,"[decompose]:  Did grow, ilast = %d\n",ilast);
      }
      fnx_pout(3,"[decompose]:  Done growing!\n");
      split_failed = 0;
      iter = 0;
      while (nproc_dendritic_rad < nproc_dendritic_below) {
        // Try to split the biggest proc transversely into two parts
        if (!split_proc(proc,npmin,npmin+nproc_dendritic_rad)) {
          split_failed = 1;
          fnx_pout(3,"[decompose]:  Unable to split dendritic below procs further!\n");
          break;
        }
        fnx_pout(5,"[decompose]:  Did split iter=%d\n",iter);
        iter++;
        nproc_dendritic_rad++;
      }
      fnx_pout(3,"[decompose]:  Done splitting dendritic below procs!  nassigned=%d\n",nproc_dendritic_rad);

      if (final) break;

      // If the split failed, continue to the next trial
      if (split_failed || (nproc_dendritic_rad != nproc_dendritic_below)) {
        // If this was the last trial, go back to the best trial so far
        if (iskip==iskip_max-1) {
          iskip = best_iskip-1;
          final = 1;
        }
        continue;
      }

      mymax = 0;
      imymax = -1;
      for (ip = 0; ip < nproc_dendritic_below; ip++) {
        cproc = npmin + ip;
        if (proc[cproc].ncells > mymax) {
          imymax = cproc;
          mymax = proc[cproc].ncells;
        }
      }
      if (mymax < lowest_max) {
        lowest_max = mymax;
        best_iskip = iskip;
      }
      if (iskip==iskip_max-1) {
        iskip = best_iskip-1;
        final = 1;
      }
      fnx_pout(3,"[decompose]:  mymax = %llu, imymax = %d, lowest_max = %llu, best_iskip = %d\n",
        mymax,imymax,lowest_max,best_iskip);

    } // end for (iskip = MIN_IBLOCK_SIZE; iskip < iskip_max; iskip++)
    cproc = npmin + nproc_dendritic_below - 1;
    istart_cart = proc[cproc].istop[0];

  } else {
    // No dendritic grid
    istart_cart = 0;
  }  // end if (istop_krefine_local > 0)
  fnx_pout(1,"[decompose]:  Done with dendritic grid below procs!\n");


  // Now add procs on the Cartesian grid
  #if (NDIM==3)
  // Test that n3 is a multiple of closest[2] and that their quotient is a multiple of 4
  if (n3%closest[2] != 0 || (n3/closest[2])%4 != 0) {
    fnx_perr(0,"[decompose]:  problem in n3 decomposition for %d procs\n", numprocs);
    exit(2);
  }
  #endif
  cproc = nproc_dendritic_total;  // Cartesian procs are last
  for (i=0; i<closest[0]; i++) {
    i0 = istart_cart +  i   *(n1-istop_jrefine_local)/closest[0];  // why istart_cart instead of istop_krefine_local???
    i1 = istart_cart + (i+1)*(n1-istop_jrefine_local)/closest[0];  // why istart_cart instead of istop_krefine_local???
    for (j=0; j<closest[1]; j++) {
      j0 = jstop_krefine_local +  j   *(jstart_kcoarsen_local-jstop_krefine_local)/closest[1];
      j1 = jstop_krefine_local + (j+1)*(jstart_kcoarsen_local-jstop_krefine_local)/closest[1];
      #if (NDIM==3)
      for (k=0; k<closest[2]; k++) {
      #endif
        fnx_pout(3,"[decompose]:  Adding proc %d\n",cproc);
        proc[cproc].istart[0] = i0;
        proc[cproc].istop[0]  = i1;
        proc[cproc].istart[1] = j0;
        proc[cproc].istop[1]  = j1;
        proc[cproc].ncells    = (i1-i0)*(j1-j0);
        #if (NDIM==3)
        k00 =  k   *n3/closest[2];
        k10 = (k+1)*n3/closest[2];
        proc[cproc].istart[2] = k00;
        proc[cproc].istop[2]  = k10;
        proc[cproc].ncells   *= (k10-k00);
        #endif
        cproc++;
      #if (NDIM==3)
      }
      #endif
    }
  }
  fnx_pout(1,"[decompose]:  Done adding Cartesian procs!\n");
  
  // Link up neighboring procs
  for (ip=0; ip<numprocs; ip++) {
    if (proc[ip].nneigh > 0) {
      // HOW CAN WE EVER GET HERE????
      fnx_pout(5,"[decompose]:  Freeing stuff %d\n",ip);
      free(proc[ip].neigh_id);
      free(proc[ip].neigh_dir);
      proc[ip].nneigh = 0;
    }
  }
  for (ip=0; ip<numprocs; ip++) {
    for (jp=ip+1; jp<numprocs; jp++) {
      fnx_pout(5,"[decompose]:  Linking %d to %d\n",ip,jp);
      is_neighbor(proc,ip,jp);
    }
  }
  fnx_pout(1,"[decompose]:  Done linking!\n");

  // Check for correctness (i.e., that every proc has at least MIN_BLOCK_SIZE cells in each direction and total count is right)
  imax_cnt      = 0;
  imax_neigh    = 0;
  max_cnt       = 0;
  total_cnt     = 0;
  dendritic_cnt = 0;
  cart_cnt      = 0;
  max_neigh     = 0;
  for (ip=0; ip<numprocs; ip++) {
    i0  = proc[ip].istart[0];
    i1  = proc[ip].istop[0];
    j0  = proc[ip].istart[1]/dj[i0];
    j1  = proc[ip].istop[1]/dj[i0];
    #if (NDIM==3)
    k00 = proc[ip].istart[2]/dk[i0][j0];
    k10 = proc[ip].istop[2]/dk[i0][j0];
    k01 = proc[ip].istart[2]/dk[i0][j1-1];
    k11 = proc[ip].istop[2]/dk[i0][j1-1];
    #else
    k00 = k10 = k01 = k11 = 0;
    #endif
    total_cnt += proc[ip].ncells;
    if (ip<nproc_dendritic_total)
      dendritic_cnt += proc[ip].ncells;
    else
      cart_cnt += proc[ip].ncells;
    if (i1 - i0 < MIN_IBLOCK_SIZE) {
      proc[ip].ncells = n1*n2*n3;
      res.actual      = n1*n2*n3;
      fnx_perr(0,"[decompose]:  Error 1, ip=%d, i0=%d, i1=%d, j0=%d, j1=%d, k00=%d, k10=%d, k01=%d, k11=%d, closest[0]=%d\n", 
        ip,i0,i1,j0,j1,k00,k10,k01,k11,closest[0]);
      return 0;
    }
    if (j1-j0 < MIN_JBLOCK_SIZE && ncart > 0 && global_grid_dims[1]/dj[i0] >= MIN_JBLOCK_SIZE) {
      proc[ip].ncells = n1*n2*n3;
      res.actual      = n1*n2*n3;
      fnx_perr(0,"[decompose]:  Error 2A, ip=%d, i0=%d, i1=%d, j0=%d, j1=%d, k00=%d, k10=%d, k01=%d, k11=%d\n",
        ip,i0,i1,j0,j1,k00,k10,k01,k11);
      return 0;
    }
    #if (NDIM==3)
    if ((k10-k00 < MIN_KBLOCK_SIZE && global_grid_dims[2]/dk[i0][j0  ] >= MIN_KBLOCK_SIZE) 
     && (k11-k01 < MIN_KBLOCK_SIZE && global_grid_dims[2]/dk[i0][j1-1] >= MIN_KBLOCK_SIZE)) {
      proc[ip].ncells = n1*n2*n3;
      res.actual      = n1*n2*n3;
      fnx_perr(0,"[decompose]:  Error 2B, ip=%d, i0=%d, i1=%d, j0=%d, j1=%d, k00=%d, k10=%d, k01=%d, k11=%d\n",
        ip,i0,i1,j0,j1,k00,k10,k01,k11);
      return 0;
    }
    #endif
    if (proc[ip].ncells > max_cnt) {
      max_cnt = proc[ip].ncells;
      imax_cnt = ip;
    }
    if (proc[ip].nneigh > max_neigh) {
      max_neigh  = proc[ip].nneigh;
      imax_neigh = ip;
    }
  }
  if (max_cnt != n1*n2*n3 && total_cnt != ncells_active_total) {
    fnx_perr(0,"[decompose]:  Error 3, %llu %llu %llu %llu %llu %llu,  %d %d %d %d %d\n",
      dendritic_cnt, ncells_dendritic_total, cart_cnt, ncells_cart_total, total_cnt, ncells_active_total, 
      closest[0], closest[1], closest[2], numprocs, nproc_dendritic_total);
    max_cnt = n1*n2*n3;
    res.actual = max_cnt;
    return 0;
  }
  fnx_pout(1,"[decompose]:  Done checking correctness!\n");

  res.nproc_dendritic = nproc_dendritic_total;
  res.nproc_cart      = nproc_cart_best;
  res.ideal           = (1.0*ncells_active_total)/numprocs;
  res.actual          = max_cnt;
  res.worst           = imax_cnt;
  res.max_neigh       = max_neigh;
  DLOOP { res.cart_dims[dd] = closest[dd]; }
  
  fnx_pout(1,"[decompose]:  res.nproc_dendritic = %d\n",res.nproc_dendritic);
  fnx_pout(1,"[decompose]:  res.nproc_cart      = %d\n",res.nproc_cart);
  fnx_pout(1,"[decompose]:  res.ideal           = %f\n",res.ideal);
  fnx_pout(1,"[decompose]:  res.actual          = %llu\n",res.actual);
  fnx_pout(1,"[decompose]:  res.worst           = %d\n",res.worst);
  fnx_pout(1,"[decompose]:  res.max_neigh       = %lu\n",res.max_neigh);
  
  // Uncomment the following lines to view the decomposition in VisIt
#if 0
  if (mpi_io_proc()) dump_grid_VTK(proc);
  MPI_Barrier(MPI_COMM_WORLD);
  exit(1);
#endif

  return 1;
}
#endif /* NDIM>1 */

void decompose1d(proc_info *proc)
{
  double di = (1.0*n1)/numprocs;
  int i;

  for (i=0; i < numprocs; i++) {
    proc[i].istart[0] = i*di;
    proc[i].istop[0] = (i+1)*di;
    proc[i].nneigh = 2;
    proc[i].neigh_id = malloc_rank1(2, sizeof(uint32_t));
    proc[i].neigh_id[0] = i-1;
    proc[i].neigh_id[1] = i+1;
    proc[i].neigh_dir = malloc_rank1(2, sizeof(int));
    proc[i].neigh_dir[0] = -1;
    proc[i].neigh_dir[1] = 1;
  }
  if (!periodic[0]) {
    proc[0].nneigh = 1;
    proc[0].neigh_id[0] = 1;
    proc[numprocs-1].nneigh = 1;
    proc[numprocs-1].neigh_id[0] = numprocs-2;
  } else {
    proc[0].neigh_id[0] = numprocs-1;
    proc[numprocs-1].neigh_id[1] = 0;
  }

  return;
}

#if (NDIM>1)
int test_proc_info(int nproc, proc_info* proc)
{
  int ip,i,j,k,i0,i1,j0,j1,k0,k1;
  int j00,j01,j10,j11,k00,k01,k10,k11;
  uint64_t local_cnt,total_cnt;
  uint64_t ncells_active_below,ncells_active_above,ncells_active_total;
  int istop_krefine_local,istop_jrefine_local;
  int n,nid,status;
  int ok=1;
  uint64_t shell_cnt1,shell_cnt2;
  
  // Check total count
  fnx_pout(0,"[test_proc_info]:  checking total cell count...");
  total_cnt = 0;
  status = 1;
  for (ip=0; ip<nproc; ip++) {
    local_cnt = 0;
    i0 = proc[ip].istart[0];
    i1 = proc[ip].istop[0];
    for (i=i0; i<i1; i++) {
      j0 = proc[ip].istart[1]/dj[i];
      j1 = proc[ip].istop[1]/dj[i];
      #if (NDIM==2)
      local_cnt += (j1-j0);
      #elif (NDIM==3)      
      for (j=j0; j<j1; j++) local_cnt += (proc[ip].istop[2] - proc[ip].istart[2])/dk[i][j];
      #endif
    }
    if (local_cnt != proc[ip].ncells) {
      fnx_pout(0,"failed! ip=%d, local_cnt=%llu, ncells=%llu\n",ip,local_cnt,proc[ip].ncells);
      status = 0;
    }
    total_cnt += local_cnt;
  }
  if (status) fnx_pout(0,"passed! total_cnt = %llu\n",total_cnt);
  ok *= status;
    
  // Compute some stuff
  istop_krefine_local = (istop_krefine > 0) ? MAX(MAX(istop_jrefine,istop_krefine),MIN_IBLOCK_SIZE) : 0;
  istop_jrefine_local = (istop_jrefine > 0) ? MAX(MAX(istop_jrefine,istop_krefine_local),MIN_IBLOCK_SIZE) : 0;    
  ncells_active_below = 0;
  for (i=0; i<istop_jrefine_local; i++) {
    #if (NDIM==2)
    ncells_active_below += n2/dj[i];
    #elif (NDIM==3)
    for (j=0; j<n2/dj[i]; j++) ncells_active_below += n3/dk[i][j];
    #endif
  }
  
  ncells_active_above = 0;
  for (i=istop_jrefine_local; i<n1; i++) {
    #if (NDIM==2)
    ncells_active_above += n2/dj[i];
    #elif (NDIM==3)
    for (j=0; j<n2/dj[i]; j++) ncells_active_above += n3/dk[i][j];
    #endif
  }
  ncells_active_total = ncells_active_below + ncells_active_above;
      
  fnx_pout(1,"[test_proc_info]:  istop_jrefine_local = %d, istop_krefine_local = %d\n",istop_jrefine_local,istop_krefine_local);
  fnx_pout(1,"[test_proc_info]:  ncells_active_above = %llu, ncells_active_below = %llu, ncells_active_total = %llu\n",
    ncells_active_above,ncells_active_below,ncells_active_total);
  
  // Check for overlapping neighbors
  fnx_pout(0,"[test_proc_info]:  checking for overlapping neighbors...");
  status = 1;
  for (ip=0; ip<nproc; ip++) {
    for (n=0; n<proc[ip].nneigh; n++) {
      nid = proc[ip].neigh_id[n];
      if (proc[ip].neigh_dir[n]==-1) {
        #if (NDIM==3 && GEOM==SPHERICAL)
        if ((proc[ip].istop[2] - (proc[nid].istart[2]+n3/2)%n3)*(proc[ip].istart[2] - (proc[nid].istop[2]+n3/2-1)%n3+1) < 0
          && (proc[ip].istart[0]==0 && proc[nid].istart[0]==0))  continue;
        #else
        if (proc[ip].istart[0] != proc[nid].istop[0]%n1) {
          fnx_pout(0,"failed! ip=%d, nid=%d, dir=-1, istart=%d, istop=%d\n",ip,nid,proc[ip].istart[0],proc[nid].istop[0]);
          status = 0;
        }
        #endif
      } else if (proc[ip].neigh_dir[n]==+1) {
        #if (NDIM==3 && GEOM==SPHERICAL)
        if ((proc[ip].istop[2] - (proc[nid].istart[2]+n3/2)%n3)*(proc[ip].istart[2] - (proc[nid].istop[2]+n3/2-1)%n3+1) < 0
          && (proc[ip].istart[0]==0 && proc[nid].istart[0]==0))  continue;
        #else
        if (proc[ip].istop[0]%n1 != proc[nid].istart[0]) {
          fnx_pout(0,"failed! ip=%d, nid=%d, dir=+1, istart=%d, istop=%d\n",ip,nid,proc[ip].istop[0],proc[nid].istart[0]);
          status = 0;
        }
        #endif
      } else if (proc[ip].neigh_dir[n]==-2) {
        #if (NDIM==3 && GEOM==SPHERICAL)
        if ((proc[ip].istop[2] - (proc[nid].istart[2]+n3/2)%n3)*(proc[ip].istart[2] - (proc[nid].istop[2]+n3/2-1)%n3+1) < 0
          && (proc[ip].istart[1]==0 && proc[nid].istart[1]==0))  continue;
        #else
        if (proc[ip].istart[1] != proc[nid].istop[1]%n2) {
          fnx_pout(0,"failed! ip=%d, nid=%d, dir=-2, istart=%d, istop=%d\n",ip,nid,proc[ip].istart[1],proc[nid].istop[1]);
          status = 0;
        }
        #endif
      } else if (proc[ip].neigh_dir[n]==+2) {
        #if (NDIM==3 && GEOM==SPHERICAL)
        if ((proc[ip].istop[2] - (proc[nid].istart[2]+n3/2)%n3)*(proc[ip].istart[2] - (proc[nid].istop[2]+n3/2-1)%n3+1) < 0
          && (proc[ip].istop[1]==n2 && proc[nid].istop[1]==n2))  continue;
        #else
        if (proc[ip].istop[1]%n2 != proc[nid].istart[1]) {
          fnx_pout(0,"failed! ip=%d, nid=%d, dir=+2, istart=%d, istop=%d\n",ip,nid,proc[ip].istop[1],proc[nid].istart[1]);
          status = 0;
        }
        #endif
      }
      #if (NDIM==3)
      else if (proc[ip].neigh_dir[n]==-3) {
        if (proc[ip].istart[2] != proc[nid].istop[2]%n3) {
          fnx_pout(0,"failed! ip=%d, nid=%d, dir=-3, istart=%d, istop=%d\n",ip,nid,proc[ip].istart[2],proc[nid].istop[2]);
          status = 0;
        }
      } else if (proc[ip].neigh_dir[n]==+3) {
        if (proc[ip].istop[2]%n3 != proc[nid].istart[2]) {
          fnx_pout(0,"failed! ip=%d, nid=%d, dir=+3, istart=%d, istop=%d\n",ip,nid,proc[ip].istop[2],proc[nid].istart[2]);
          status = 0;
        }
      }
      #endif
    }
  }
  if (status) fnx_pout(0,"passed!\n");
  ok *= status;
  
  // check radial shells
  fnx_pout(0,"[test_proc_info]:  checking radial shells...");
  status = 1;
  for (i=0; i<n1; i++) {
    shell_cnt1 = 0;
    for (ip=0; ip<nproc; ip++) {
      if (i >= proc[ip].istart[0] && i < proc[ip].istop[0]) {
        j0 = proc[ip].istart[1]/dj[i];
        j1 = proc[ip].istop[1]/dj[i];
        #if (NDIM==2)
        shell_cnt1 += (j1-j0);
        #elif (NDIM==3)      
        for (j=j0; j<j1; j++) shell_cnt1 += (proc[ip].istop[2] - proc[ip].istart[2])/dk[i][j];
        #endif
      }
    }
    
    shell_cnt2 = 0;
    #if (NDIM==2)
    shell_cnt2 += n2/dj[i];
    #elif (NDIM==3)      
    for (j=0; j<n2/dj[i]; j++) shell_cnt2 += n3/dk[i][j];
    #endif
    
    if (shell_cnt1 != shell_cnt2) {
      fnx_pout(0,"failed, i=%d, shell_cnt1=%llu, shell_cnt2=%llu\n",i,shell_cnt1,shell_cnt2);
      for (ip=0; ip<nproc; ip++) {
        if (i >= proc[ip].istart[0] && i < proc[ip].istop[0]) {
        #if (NDIM==2)
          fnx_pout(0,"contributing procs:  ip=%d, istart[0]=%d, istop[0]=%d, istart[1]=%d, istop[1]=%d\n",
            ip,proc[ip].istart[0],proc[ip].istop[0],proc[ip].istart[1],proc[ip].istop[1]);
        #endif
        #if (NDIM==3)
          fnx_pout(0,"contributing procs:  ip=%d, istart[0]=%d, istop[0]=%d, istart[1]=%d, istop[1]=%d, istart[2]=%d, istop[2]=%d\n",
            ip,proc[ip].istart[0],proc[ip].istop[0],proc[ip].istart[1],proc[ip].istop[1],proc[ip].istart[2],proc[ip].istop[2]);
        #endif
        }
      }
      status = 0;
    }
  }
  if (status) fnx_pout(0,"passed!\n");
  ok *= status;
  
  // check refinement/coarsening domain criteria
  fnx_pout(0,"[test_proc_info]:  checking refinement/coarsening domain criteria...");
  status = 1;  
  for (ip=0; ip<nproc; ip++) {
    i0 = proc[ip].istart[0];
    i1 = proc[ip].istop[0];
    j00 = JS(i0,proc[ip].istart[1]);
    j01 = JS(i0,proc[ip].istop[1]);
    j10 = JS(i1,proc[ip].istart[1]);
    j11 = JS(i1,proc[ip].istop[1]);
    if (j11-j10 != (dj[i0]/dj[i1])*(j01-j00)) {
      fnx_pout(0,"failed, for proc %d, j11-j10 = %d, but should be (dj[i0]/dj[i1])*(j01-j00) = %d\n",
        ip,j11-j10,(dj[i0]/dj[i1])*(j01-j00));
      status = 0;
    }
    #if NDIM==3
    for (i=i0; i<i1; i++) {
      j0 = JS(i,proc[ip].istart[1]);
      j1 = JS(i,proc[ip].istop[1]);
      k00 = KS(i,j0,proc[ip].istart[2]);
      k01 = KS(i,j0,proc[ip].istop[2]);
      k10 = KS(i,j1,proc[ip].istart[2]);
      k11 = KS(i,j1,proc[ip].istop[2]);
      if (dk[i][j0] > dk[i][j1] && k11-k10 != (dk[i][j0]/dk[i][j1])*(k01-k00)) {
        fnx_pout(0,"failed, for proc %d at i=%d, k11-k10 = %d, but should be (dk[i][j0]/dk[i][j1])*(k01-k00) = %d\n",
          ip,i,k11-k10,(dk[i][j0]/dk[i][j1])*(k01-k00));
        fnx_pout(0,"dk[i][j0]=%d, dk[i][j1]=%d, istart[2]=%d, istop[2]=%d\n",dk[i][j0],dk[i][j1],proc[ip].istart[2],proc[ip].istop[2]);
        status = 0;
      } else if (dk[i][j0] <= dk[i][j1] && k01-k00 != (dk[i][j1]/dk[i][j0])*(k11-k10)) {
        fnx_pout(0,"failed, for proc %d at i=%d, k01-k00 = %d, but should be (dk[i][j1]/dk[i][j0])*(k11-k10) = %d\n",
          ip,i,k01-k00,(dk[i][j1]/dk[i][j0])*(k11-k10));
        fnx_pout(0,"dk[i][j0]=%d, dk[i][j1]=%d, istart[2]=%d, istop[2]=%d\n",dk[i][j0],dk[i][j1],proc[ip].istart[2],proc[ip].istop[2]);
        status = 0;
      }
    }
    #endif
  }
  if (status) fnx_pout(0,"passed!\n");  
  ok *= status;
  
  if (ok) fnx_pout(0,"[test_proc_info]:  ALL TESTS PASSED!\n");

  return ok;
}

void read_proc_info(int nproc, proc_info* proc, result* res)
{
  int ip,dd,n;
  int _nproc,_ndim,_n1,_n2,_n3,_ip;
  FILE* fp=NULL;
  double efficiency;
  char filename[512];
  
  if (mpi_io_proc()) {
    snprintf(filename,512,"%s/proc_info_%d_%d_%d_%d_%d.txt",decomp_path,nproc,NDIM,n1,n2,n3);
    fp = fopen(filename,"r");
    if (fp==NULL) {
      fnx_pout(0,"\n[read_proc_info]:  File %s not found!\n",filename);
      exit(1);
    }
    fscanf(fp,"%d %d %d %d %d %d %d %"SCNu32" %d %f %"SCNu64" %lf\n",
      &_nproc,&_ndim,&_n1,&_n2,&_n3,&res->nproc_dendritic,&res->nproc_cart,
      &res->max_neigh,&res->worst,&res->ideal,&res->actual,&efficiency);

    if ((_nproc != nproc) || (_ndim != NDIM) || (_n1 != n1) || (_n2 != n2) || (_n3 != n3)) {
      fnx_pout(0,"\n[read_proc_info]: ERROR!  Need %d %d %d %d %d, but have %d %d %d %d %d...\n",
        nproc,NDIM,n1,n2,n3,_nproc,_ndim,_n1,_n2,_n3);
      exit(1);
    }
  }
  #if (USE_MPI==TRUE)
  MPI_Bcast(&res->nproc_dendritic,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&res->nproc_cart,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&res->max_neigh,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
  MPI_Bcast(&res->worst,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(&res->ideal,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&res->actual,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
  MPI_Bcast(&efficiency,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  #endif

  for (ip=0; ip<nproc; ip++) {
    if (mpi_io_proc()) fscanf(fp,"%d %"SCNu64" ",&_ip,&proc[ip].ncells);
    #if (USE_MPI==TRUE)    
    MPI_Bcast(&proc[ip].ncells,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
    #endif
    DLOOP { 
      if (mpi_io_proc()) fscanf(fp,"%d %d ",&proc[ip].istart[dd],&proc[ip].istop[dd]);
      #if (USE_MPI==TRUE)
      MPI_Bcast(&proc[ip].istart[dd],1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&proc[ip].istop[dd],1,MPI_INT,0,MPI_COMM_WORLD);
      #endif
    }
    if (mpi_io_proc()) fscanf(fp,"%"SCNu32" ",&proc[ip].nneigh);
    #if (USE_MPI==TRUE)
    MPI_Bcast(&proc[ip].nneigh,1,MPI_UINT32_T,0,MPI_COMM_WORLD);
    #endif
    proc[ip].neigh_dir = malloc_rank1(proc[ip].nneigh,sizeof *proc[ip].neigh_dir);
    proc[ip].neigh_id  = malloc_rank1(proc[ip].nneigh,sizeof *proc[ip].neigh_id);
        
    for (n=0; n<proc[ip].nneigh; n++) {
      if (mpi_io_proc()) fscanf(fp,"%d %"SCNu32" ",&proc[ip].neigh_dir[n],&proc[ip].neigh_id[n]);
      #if (USE_MPI==TRUE)
      MPI_Bcast(&proc[ip].neigh_dir[n],1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&proc[ip].neigh_id[n],1,MPI_UINT32_T,0,MPI_COMM_WORLD);
      #endif
    }
    if (mpi_io_proc()) fscanf(fp,"\n");
  }
  
  if (mpi_io_proc()) {
    fclose(fp);
    fp = NULL;
  }
}

void write_proc_info(int nproc, proc_info* proc, result* res)
{
  int ip,dd,n;
  FILE* fp=NULL;
  double efficiency;
  char filename[512];

  if (mpi_io_proc() || decomp_only) {
    efficiency = 1.0 - (res->actual-res->ideal)/(1.0*res->ideal);
  
    snprintf(filename,512,"%s/proc_info_%d_%d_%d_%d_%d.txt",decomp_path,nproc,NDIM,n1,n2,n3);
    fprintf(stderr,"[write_proc_info]:  writing decomposition to file %s ...",filename);
    // before writing, check if file exists...
    fp = fopen(filename,"r");
    if (fp) {
      fnx_pout(0,"[write_proc_info]:  File %s exists!  Delete file to overwrite...\n",filename);
      fclose(fp);
    } else {
      fp = fopen(filename,"w");
      fprintf(fp,"%d %d %d %d %d %d %d %"SCNu32" %d %f %"SCNu64" %f\n",
        nproc,NDIM,n1,n2,n3,res->nproc_dendritic,res->nproc_cart,res->max_neigh,res->worst,res->ideal,res->actual,efficiency);
  
      for (ip=0; ip<nproc; ip++) {
        // After header line, each line is:  ncells istart[0] istop[0] ... nneigh neigh_dir[0] neigh_id[0] ...
        fprintf(fp,"%d %"SCNu64" ",ip,proc[ip].ncells);
        DLOOP { fprintf(fp,"%d %d ",proc[ip].istart[dd],proc[ip].istop[dd]); }
        fprintf(fp,"%"SCNu32" ",proc[ip].nneigh);
        for (n=0; n<proc[ip].nneigh; n++) {
          fprintf(fp,"%d %"SCNu32" ",proc[ip].neigh_dir[n],proc[ip].neigh_id[n]);
        }
        fprintf(fp,"\n");
      }
  
      fclose(fp);
      fp = NULL;
    }
    fprintf(stderr,"done!\n");
  }
}
#endif /* NDIM>1 */

proc_info* split_domain(int *ideal, int *actual, int *worst)
{
  proc_info* proc=NULL;
  int ip,dd,pstop[NDIM],ibest,neigh,ok,np,numprocs_actual;
  double efficiency,best;
  uint32_t nmin;
  FILE* fp=NULL;

  #if (NDIM>1)
  if (decomp_only) {
    // Just compute the decompositions and return
    best = -999.0;
    ibest = 0;
    // Dump the efficiency data to a file
    fp = fopen("efficiency.txt","a");
    numprocs_actual = numprocs;
    for (np=decomp_npmin; np<=decomp_npmax; np+=decomp_npskip*numprocs_actual) {
      numprocs = np + myrank*decomp_npskip;
      proc = malloc_rank1(numprocs+myrank, sizeof(proc_info));
      for (ip=0; ip<numprocs; ip++) proc[ip].nneigh = 0;
      
      ok = decompose(proc);
      if (!ok) {
        fnx_pout(0,"\n[split_domain]:  np=%d, FAILED!!\n",numprocs);
        exit(1);
      } else {
        efficiency = 1.0 - (res.actual-res.ideal)/(1.0*res.ideal);
        // It's not wise to have all procs write to same file, but maybe we can get away with it...
        fprintf(fp,"%d %d %d %"SCNu32" %d %f %"SCNu64" %f\n",
          numprocs,res.nproc_dendritic,res.nproc_cart,res.max_neigh,res.worst,res.ideal,res.actual,efficiency);
        fflush(fp);
        fprintf(stderr,"\n[split_domain]:  myrank=%d, np=%d, np_dend=%d, np_cart=%d, max_neigh=%"PRIu32", worst=%d, ideal=%f, actual=%"SCNu64", eff=%f\n",
          myrank,numprocs,res.nproc_dendritic,res.nproc_cart,res.max_neigh,res.worst,res.ideal,res.actual,efficiency);

        #if (NDIM==3)
        fprintf(stderr,"[split_domain]:  worst proc has %d neighbors with %"SCNu64" cells in (%d, %d) x (%d, %d) x (%d, %d)\n",
          proc[res.worst].nneigh,proc[res.worst].ncells,
          proc[res.worst].istart[0],proc[res.worst].istop[0],
          proc[res.worst].istart[1],proc[res.worst].istop[1],
          proc[res.worst].istart[2],proc[res.worst].istop[2]);
        #else
        fprintf(stderr,"[split_domain]:  worst proc has %d neighbors with %"SCNu64" cells in (%d, %d) x (%d, %d)\n",
          proc[res.worst].nneigh,proc[res.worst].ncells,
          proc[res.worst].istart[0],proc[res.worst].istop[0],
          proc[res.worst].istart[1],proc[res.worst].istop[1]);
        #endif
        
        if (efficiency>best) {
          best = efficiency;
          ibest = numprocs;
        }
      }
      
      ok = test_proc_info(numprocs,proc);      
      if (!ok) {
        fnx_pout(0,"\n[split_domain]:  np=%d, FAILED!!\n",numprocs);
        exit(1);
      }

      write_proc_info(numprocs,proc,&res);
      
      free(proc);
      #if (USE_MPI==TRUE)
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
    } /* end for */
    #if (USE_MPI==TRUE)
    double global_best;
    MPI_Allreduce(&best, &global_best, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    ibest = (best==global_best) ? ibest : -1;
    best  = global_best;
    MPI_Allreduce(MPI_IN_PLACE, &ibest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);    
    #endif
    fnx_pout(0,"[split_domain]:  best numprocs=%d, best efficiency=%f\n",ibest,best);
    fclose(fp);
    exit(1);
  } /* end if (decomp_only) */
  #endif /* NDIM>1 */
  
  proc = malloc_rank1(numprocs, sizeof *proc);

  #if (NDIM==1)
  decompose1d(proc);
  if (numprocs == 1) proc[0].nneigh = 0;
  #else
  if (decomp_from_file) {
    read_proc_info(numprocs,proc,&res);
    ok = test_proc_info(numprocs,proc);      
    if (!ok) {
      fnx_pout(0,"\n[split_domain]:  np=%d, FAILED!!\n",numprocs);
      exit(1);
    }
  } else {

    for (ip=0; ip<numprocs; ip++) proc[ip].nneigh = 0;

    if (numprocs==1) {
      ND_ARRAY_SET(pstop,n1,n2,n3);
      DLOOP {
        proc[0].istart[dd] = 0;
        proc[0].istop[dd]  = pstop[dd];
      }
      proc[0].nneigh = 0;
      return proc;
    }

    ok = decompose(proc);
    if (!ok) {
      fnx_pout(0,"\n[split_domain]:  np=%d, FAILED!!\n",numprocs);
      exit(1);
    }
  
    if (res.actual == n1*n2*n3) {
      fprintf(stderr,"[split_domain]:  Domain decomposition failed!\n");
      exit(3);
    }
    
    ok = test_proc_info(numprocs,proc);      
    if (!ok) {
      fnx_pout(0,"\n[split_domain]:  np=%d, FAILED!!\n",numprocs);
      exit(1);
    }
    
    write_proc_info(numprocs,proc,&res);    
  }
  
  *ideal  = res.ideal;
  *actual = res.actual;
  *worst  = res.worst;
  #endif

  return proc;
}
