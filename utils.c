#include "decs.h"
#include "constants.h"

// some routines for numerical integration
// just call romb(xmin,xmax,integrand)
// can't handle singularities in the integrand
#define ROMBERG_MAX	  15	// should not go above 15
#define ROMBERG_TOL	  1.e-15
#define ROMBERG_ATOL	1.e-15

// #define MAX_LOG2_ASPECT 1.0  // 0 <= log2(aspect) < 1
#define MAX_LOG2_ASPECT 0.5  // -0.5 <= log2(aspect) < 0.5, recommended

double rtrap(double x0, double x1, int rstage, double (*f)(double x))
{
	unsigned long long int nseg = 1<<rstage;
	unsigned long long int i;
	double dx = (x1-x0)/nseg;
  double sum = 0.0;

	if (rstage == 0) {
		return 0.5*dx*((*f)(x0) + (*f)(x1));
	}

	for (i=0; i<(1<<(rstage-1)); i++) {
		sum += (*f)(x0 + (2*i+1)*dx);
	}

	return dx*sum;
}

double romb(double x0, double x1, double (*f)(double x))
{
	double rvec[ROMBERG_MAX+1];
  double R0,R1;
	unsigned long long int fm;
  int rstage,m;

	rvec[0] = rtrap(x0,x1,0,f);
	if (isnan(rvec[0])) {
		fprintf(stderr,"caught isnan in romb: %g %g, singular integrand?\n", x0, x1);
		exit(1);
	}

	for (rstage=1; rstage<=ROMBERG_MAX; rstage++) {
		R0 = rvec[0];
		rvec[0] = 0.5*R0 + rtrap(x0,x1,rstage,f);
		for (m=1; m<=rstage; m++) {
			fm = 1<<(2*m);
			R1 = rvec[m];
			rvec[m] = 1.0/(fm-1)*(fm*rvec[m-1] - R0);
			R0 = R1;
		}
		if (fabs((rvec[rstage]-rvec[rstage-1])/rvec[rstage]) < ROMBERG_TOL || (rstage > 2 && fabs(rvec[rstage]) < ROMBERG_ATOL)) {
			return rvec[rstage];
		}
	}

	fprintf(stderr,"WARNING: Romberg integration failed to converge on %g %g: %g %g\n", x0, x1, rvec[ROMBERG_MAX], fabs((rvec[ROMBERG_MAX]-rvec[ROMBERG_MAX-1])/rvec[ROMBERG_MAX]));
	fflush(stderr);
	return rvec[ROMBERG_MAX];
}

static uint64_t xorx;
static uint64_t xors[16];
static int xor_p;

uint64_t xorshift_init(void)
{
	xorx ^= xorx >> 12; // a
	xorx ^= xorx << 25; // b
	xorx ^= xorx >> 27; // c
	return xorx * UINT64_C(2685821657736338717);
}


double fornax_rand(void)
{
	uint64_t s0 = xors[ xor_p ];
	uint64_t s1 = xors[ xor_p = ( xor_p + 1 ) & 15 ];
	s1 ^= s1 << 31;
	s1 ^= s1 >> 11;
	s0 ^= s0 >> 30;
	return (( xors[ xor_p ] = s0 ^ s1 ) * UINT64_C(1181783497276652981))/((double)UINT64_MAX);
}

void init_fornax_rand(void)
{
    xor_p = 0;
    xorx = time(NULL);
    for (int i=0;i<16;i++) xors[i] = xorshift_init();
}

double my_mod(double x, double a)
{
    return fmod(fmod(x,a)+a,a);
}

char *trim_whitespace(char *str)
{
	char *end;

	// Trim leading space
	while(isspace(*str)) str++;

	if (*str == 0) {  // All spaces?
		return str;
	}

	end = str;
	while(*end != 0 && !isspace(*end) && *end!='\n') end++;

	*end=0;	// null terminate

	return str;
}


int big_endian(void)
{
  short int n = 1;
  char *ep = (char *)&n;

  return (*ep == 0); // Returns 1 on a big endian machine
}


// swap bytes, code stolen from athena, who stole it from NEMO
void bswap(void *vdat, int len, int cnt)
{
  char tmp, *dat = (char *) vdat;
  int k;

  if (len==1)
    return;
  else if (len==2)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
      dat += 2;
    }
  else if (len==4)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
      tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
      dat += 4;
    }
  else if (len==8)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
      tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
      tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
      tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
      dat += 8;
    }
  else {  // the general SLOOOOOOOOOW case
    for (k=0; k<len/2; k++) {
      tmp = dat[k];
      dat[k] = dat[len-1-k];
      dat[len-1-k] = tmp;
    }
  }
}

void* malloc_rank1(uint64_t n1, uint64_t size)
{
  void* A=NULL;

  if ((A = malloc(n1*size)) == NULL) {
    fprintf(stderr,"[malloc_rank1]:  malloc failure on proc %d, n1=%" PRIu64 ", size=%" PRIu64 "\n",myrank,n1,size);
  }
  memset(A, 0, n1*size);

  return A;
}

/* [compute_dj]:  Compute the number of j-zones for a given i-zone
 *   Note that if you start with nj_min < 2, then nj can't fast enough
 *   to keep the ratio within a factor of 2, since we won't allow more
 *   than 2-to-1 refinement.  Finally, compute the dendritic factor dj.
 */
void compute_dj()
{
  const int nj_min=2;
  const int nj_max=n2;
  int i,s0;
  double r0,r1,rc,dr,nj_ideal,dth;

  /* Allocate arrays */
  s0  = n1 + 2*NG;
  nj  = malloc_rank1(s0, sizeof *nj);
  nj += NG;
  dj  = malloc_rank1(s0, sizeof *dj);
  dj += NG;

  #if (DENDRITIC_GRID==TRUE && NDIM>1 && GEOM==SPHERICAL)
  nj[0] = nj_min;
  dth = th_of_x(thx_info,startx[1] + n2*dx[1]) - th_of_x(thx_info,startx[1]);
  for (i=1; i<n1; i++) {
    /* Assume a uniform angular discretization */
    r0 = r_of_x(rx_info,startx[0] + (i+0.0)*dx[0]);
    rc = r_of_x(rx_info,startx[0] + (i+0.5)*dx[0]);
    r1 = r_of_x(rx_info,startx[0] + (i+1.0)*dx[0]);
    dr = r1 - r0;
    nj_ideal = rc*dth/dr;

    nj[i] = nj[i-1];
    /* Increase if aspect ratio exceeds 1 too greatly */
    if (log2(nj_ideal/nj[i]) > MAX_LOG2_ASPECT) nj[i] *= 2;
    nj[i] = MIN(MAX(nj[i], nj[i-1]), nj_max);
  }
  #else /* Not multi-D SPHERICAL, but compute anyway */
  for (i=0; i<n1; i++) {
    nj[i] = n2;
  }
  #endif

  /* Extend to ghost zones */
  for (i=-NG; i<0; i++) {
    #if (GEOM==SPHERICAL)
    /* Reflect across i=0 */
    nj[i] = nj[-i-1];
    #else
    /* Otherwise, just copy from i=0 */
    nj[i] = nj[0];
    #endif
  }
  for (i=n1; i<n1+NG; i++) nj[i] = nj[n1-1];

  /* Compute dendritic factor */
  for (i=-NG; i<n1+NG; i++) {
    dj[i] = n2/nj[i];
  }
}

/* [compute_dk]:  Compute the number of k-zones for a given i,j-zone
 *   Note that if you start with nk_min < 4, then nk can't fast enough
 *   to keep the ratio within a factor of 2, since we won't allow more
 *   than 2-to-1 refinement.  Finally, compute the dendritic factor, dk.
 */
void compute_dk()
{
  const int nk_min = 4;
  const int nk_max = n3;
  int i,j,s0,s1,cnt;
  int* tmp1=NULL;
  int* tmp2=NULL;
  double r0,r1,rc,dr,thc,dph,nk_ideal;

  /* Allocate arrays */
  s0 = n1 + 2*NG;
  s1 = n2 + 2*NG;
  tmp1 = malloc_rank1(s0*s1, sizeof *tmp1);
  tmp2 = malloc_rank1(s0*s1, sizeof *tmp1);
  nk   = malloc_rank1(s0, sizeof *nk);
  dk   = malloc_rank1(s0, sizeof *dk);
  nk  += NG;
  dk  += NG;
  cnt  = 0;
  for (i=-NG; i<s0-NG; i++) {
    nk[i] = tmp1 + cnt + NG;
    dk[i] = tmp2 + cnt + NG;
    cnt += s1;
  }

  #if (DENDRITIC_GRID==TRUE && NDIM==3 && GEOM==SPHERICAL)
  /* Compute nk for all j-zones at i=0 */
  nk[0][0] = nk_min;
  dph = n3*dx[2];
  for (j=1; j<nj[0]; j++) {
    r0  = r_of_x(rx_info,startx[0] + (0.0)*dx[0]);
    rc  = r_of_x(rx_info,startx[0] + (0.5)*dx[0]);
    r1  = r_of_x(rx_info,startx[0] + (1.0)*dx[0]);
    dr  = r1 - r0;
    thc = th_of_x(thx_info,startx[1] + (j+0.5)*dj[0]*dx[1]);
    nk_ideal = rc*sin(thc)*dph/dr;

    /* Compute nk for the first nj/2 j-zones, otherwise reflect */
    if (j < nj[0]/2) {
      nk[0][j] = nk[0][j-1];
      if (log2(nk_ideal/nk[0][j]) > MAX_LOG2_ASPECT) nk[0][j] *= 2;
      nk[0][j] = MIN(nk[0][j],nk_max);
    }
    else {
      nk[0][j] = nk[0][nj[0]-1-j];
    }
  }

  /* Compute nk for all other i,j-zones */
  for (i=1; i<n1; i++) {
    r0  = r_of_x(rx_info,startx[0] + (i+0.0)*dx[0]);
    rc  = r_of_x(rx_info,startx[0] + (i+0.5)*dx[0]);
    r1  = r_of_x(rx_info,startx[0] + (i+1.0)*dx[0]);
    dr  = r1 - r0;

    nk[i][0] = nk_min;
    for (j=1; j<nj[i]; j++) {
      thc = th_of_x(thx_info,startx[1] + (j+0.5)*dj[i]*dx[1]);
      nk_ideal = rc*sin(thc)*dph/dr;

      if (j < nj[i]/2) {
        nk[i][j] = nk[i][j-1];
        if (log2(nk_ideal/nk[i][j]) > MAX_LOG2_ASPECT) nk[i][j] *= 2;
        if (nj[i-1] == nj[i]) {
          nk[i][j] = MIN(MAX(nk[i][j], nk[i-1][j]),nk_max);
        } else {
          nk[i][j] = MIN(MAX(nk[i][j], nk[i-1][j/2]),nk_max);
        }
      }
      else {
        nk[i][j] = nk[i][nj[i]-1-j];
      }
    }
  }
  #else /* Not 3D SPHERICAL, but compute anyway */
  for (i=0; i<n1; i++) {
    for (j=0; j<nj[i]; j++) {
      nk[i][j] = n3;
    }
  }
  #endif

  /* Extend to ghost zones */
  for (i=-NG; i<0; i++) {
    for (j=0; j<nj[i]; j++) {
      #if (GEOM==SPHERICAL)
      /* Reflect across i=0 */
      nk[i][j] = nk[-i-1][j];
      #else
      /* Otherwise, just copy from i=0 */
      nk[i][j] = nk[0][j];
      #endif
    }
  }
  for (i=n1; i<n1+NG; i++) {
    for (j=0; j<nj[i]; j++) {
      nk[i][j] = nk[n1-1][j];
    }
  }
  for (i=-NG; i<n1+NG; i++) {
    for (j=-NG; j<0; j++) {
      #if (GEOM==SPHERICAL && NDIM>1)
      /* Reflect across j=0 */
      nk[i][j] = nk[i][-j-1];
      #else
      /* Otherwise, just copy from j=0 */
      nk[i][j] = nk[i][0];
      #endif
    }
  }
  for (i=-NG; i<n1+NG; i++) {
    for (j=nj[i]; j<nj[i]+NG; j++) {
      #if (GEOM==SPHERICAL && NDIM>1)
      /* Reflect across j=nj[i]-1 */
      nk[i][j] = nk[i][2*nj[i]-1-j];
      #else
      /* Otherwise, just copy from j=nj[i]-1 */
      nk[i][j] = nk[i][nj[i]-1];
      #endif
    }
  }

  /* Compute dendritic factor */
  for (i=-NG; i<n1+NG; i++) {
    for (j=-NG; j<nj[i]+NG; j++) {
      dk[i][j] = n3/nk[i][j];
    }
  }
}

void init_dendritic()
{
  int i,j,k,ii,jj,kk,i0,i1,j00,j01,j10,j11;
  int tmp,dkmin,s0,refine_error;
  double r0,r1,rc,dr,th0,th1,thc,dth,dph;
  double aspect,min_aspect,max_aspect;
  FILE *fp=NULL;

  refine_error = 0;

  #if (NDIM>1)
  if (mpi_io_proc()) fp = fopen("dendritic.txt","w");
  #else
  if (mpi_io_proc()) fp = stderr;
  #endif

  compute_dj();


  /* Determine where to stop j-refinement */
  for (i=0; i<n1; i++) {
    if (dj[i] == 1) {
      istop_jrefine = i;
      break;
    }
  }
  if (istop_jrefine==0)  istop_jrefine = -2*NG;
  else istop_jrefine += NG;
	if (mpi_io_proc()) fprintf(fp,"istop_jrefine = %d\n", istop_jrefine);

  /* Determine max j-refinement level */
  tmp = dj[0];
  max_jrefine_level = 0;
  while (tmp >>= 1) max_jrefine_level++;
	if (mpi_io_proc()) fprintf(fp,"max_jrefine_level = %d\n", max_jrefine_level);

  /* Print dj to file */
  if (mpi_io_proc()) {
    for (i=-NG; i<=istop_jrefine; i++) {
      fprintf(fp,"nj[%d]=%d, dj[%d]=%d\n", i,nj[i],i,dj[i]);
    }
    fflush(fp);
  }

  #if (GEOM==SPHERICAL || GEOM==CYLINDRICAL)
  /* Test aspect ratios in theta-direction */
  min_aspect = 10.0;
  max_aspect = 0.0;
  for (i=0; i<n1; i++) {
    r0  = r_of_x(rx_info,startx[0] + (i+0.0)*dx[0]);
    rc  = r_of_x(rx_info,startx[0] + (i+0.5)*dx[0]);
    r1  = r_of_x(rx_info,startx[0] + (i+1.0)*dx[0]);
    dr  = r1 - r0;
    for (j=0; j<nj[i]; j++) {
      th0 = th_of_x(thx_info,startx[1] + (j+0.0)*dj[i]*dx[1]);
      th1 = th_of_x(thx_info,startx[1] + (j+1.0)*dj[i]*dx[1]);
      dth = th1 - th0;
      aspect = rc*dth/dr;
      min_aspect = MIN(min_aspect,aspect);
      max_aspect = MAX(max_aspect,aspect);
      // printf("i=%d, j=%d, nj[i]=%d, aspect=%f, min_aspect=%f, max_aspect=%f\n",i,j,nj[i],aspect,min_aspect,max_aspect);
    }
  }
  if (mpi_io_proc()) fprintf(fp,"Aspect ratios in the theta-direction:  MIN=%f, MAX=%f\n",min_aspect,max_aspect);
  #endif

  compute_dk();

  s0 = n1 + 2*NG;
  jstop_krefine    = malloc_rank1(s0,sizeof *jstop_krefine);
  jstart_kcoarsen  = malloc_rank1(s0,sizeof *jstart_kcoarsen);
  jstop_krefine   += NG;
  jstart_kcoarsen += NG;

  /* Determine where to start/stop k-refinement */
  istop_krefine = 0;
  for (j=-NG; j<n2+NG; j++) {
    for (i=-NG; i<n1+NG; i++) {
      jj = j/dj[i];
      if (dk[i][jj] == 1) {
        istop_krefine = MAX(istop_krefine,i);
        break;
      }
    }
  }
  if (istop_krefine==0)  istop_krefine = -2*NG;
  else istop_krefine += NG;
  if (mpi_io_proc()) fprintf(fp,"istop_krefine = %d\n",istop_krefine);

  for (i=-NG; i<n1+NG; i++) {
    jstop_krefine[i] = 0;
    for (j=0; j<nj[i]/2; j++) {
      if (dk[i][j] == dk[i][nj[i]/2]) {
        jstop_krefine[i] = j;
        break;
      }
    }
    if (jstop_krefine[i]==0) jstop_krefine[i] = -2*NG;
    else jstop_krefine[i] += NG;
    jstop_krefine[i] = MIN(jstop_krefine[i],nj[i]/2);

    jstart_kcoarsen[i] = nj[i];
    for (j=nj[i]; j>nj[i]/2; j--) {
      if (dk[i][j-1] == dk[i][nj[i]/2]) {
        jstart_kcoarsen[i] = j;
        break;
      }
    }
    if (jstart_kcoarsen[i]==nj[i])  jstart_kcoarsen[i] = nj[i]+2*NG;
    else jstart_kcoarsen[i] -= NG;
    jstart_kcoarsen[i] = MAX(jstart_kcoarsen[i],nj[i]/2);

    /* Determine max k-refine level for each i */
    tmp = dk[i][0];
    max_krefine_level = 0;
    while (tmp >>= 1) max_krefine_level++;
  }

  if (mpi_io_proc()) fprintf(fp,"max_krefine_level = %d\n",max_krefine_level);
  /* Print info on where dendritic refinement/coarsening in k stops/starts for each i */
  for (i=0; i<=istop_krefine; i++) {
    if (mpi_io_proc()) fprintf(fp,"jstop_krefine[%3d] = %3d, jstart_kcoarsen[%3d] = %3d, %6.2f%% of j-zones are unrefined in k\n",
    i,jstop_krefine[i],i,jstart_kcoarsen[i],100.0*(jstart_kcoarsen[i]-jstop_krefine[i])/nj[i]);
  }

  /* Print dk to file */
  if (mpi_io_proc()) {
    for (i=-NG; i<=istop_krefine; i++) {
      for (j=-NG; j<nj[i]+NG; j++) {
        if (mpi_io_proc()) fprintf(fp,"nk[%d][%d]=%d, dk[%d][%d]=%d\n", i,j,nk[i][j],i,j,dk[i][j]);
      }
    }
    fflush(fp);
  }

  #if (GEOM==SPHERICAL || GEOM==CYLINDRICAL)
  /* Test aspect ratios in phi-direction */
  min_aspect = 10.0;
  max_aspect = 0.0;
  for (i=0; i<n1; i++) {
    r0  = r_of_x(rx_info,startx[0] + (i+0.0)*dx[0]);
    rc  = r_of_x(rx_info,startx[0] + (i+0.5)*dx[0]);
    r1  = r_of_x(rx_info,startx[0] + (i+1.0)*dx[0]);
    dr  = r1 - r0;
    for (j=0; j<nj[i]; j++) {
      thc = th_of_x(thx_info,startx[1] + (j+0.5)*dj[i]*dx[1]);
      dph = (n3*dx[2])/nk[i][j];
      aspect = rc*sin(thc)*dph/dr;
      min_aspect = MIN(min_aspect,aspect);
      max_aspect = MAX(max_aspect,aspect);
      // printf("i=%d, i=%d, nk[i][j]=%d, aspect=%f, min_aspect=%f, max_aspect=%f\n",i,j,nk[i][j],aspect,min_aspect,max_aspect);
    }
  }
  if (mpi_io_proc()) fprintf(fp,"Aspect ratios in the   phi-direction:  MIN=%f, MAX=%f\n",min_aspect,max_aspect);
  #endif

#if 0
  if (mpi_io_proc()) {
    /* "Graphically" represent results */
    fprintf(stderr,"\n");
    for (i=0; i<n1; i++) {
      fprintf(stderr,"%03d ",i);
      for (j=0; j<nj[i]; j++) {
        #if 1
        /* Print multiple copies */
        for (jj=0; jj<n2/nj[i]; jj++) {
          fprintf(stderr,"%d",(int)log2(dk[i][j]));
        }
        #else
        /* Print just one copy */
        fprintf(stderr,"%d",(int)log2(dk[i][j]));
        #endif
      }
      fprintf(stderr,"\n");
    }
  }
#endif

  #if (GEOM==SPHERICAL || GEOM==CYLINDRICAL)
  /* Test for a few forbidden cases (they should never happen!) */
  refine_error = 0;
  for (ii=0; ii<n1; ii++) {
    for (jj=0; jj<n2/dj[ii]; jj++) {
      for (kk=0; kk<n3/dk[ii][jj]; kk++) {
        i0   = ii;
        i1   = ii + 1;

        j00  = jj;
        j01  = jj + 1;
        j10  = jj*dj[i0]/dj[i1];
        j11  = j10 + 1;

        if (dj[i0]/dj[i1]==2) {
          if (dk[i0][j00]/dk[i1][j10]==4 || dk[i0][j00]/dk[i1][j11]==4) {
            fprintf(stdout,"[init_dendritic]:  Refinement is 4-to-1 in cell (%d,%d,%d)...\n",ii,jj,kk);
            fprintf(stdout,"[init_dendritic]:  dj[%d]=%d, dj[%d]=%d, dk[%d][%d]=%d, dk[%d][%d]=%d, dk[%d][%d]=%d\n",
              i0,dj[i0],i1,dj[i1],i0,j00,dk[i0][j00],i1,j10,dk[i1][j10],i1,j11,dk[i1][j11]);
            refine_error = 1;
          }
          if (dk[i0][j00]/dk[i1][j10]==0 || dk[i0][j00]/dk[i1][j11]==0) {
            // fprintf(stdout,"[init_dendritic]:  Refinement is 1/2-to-1 in cell (%d,%d,%d)...\n",ii,jj,kk);
            // fprintf(stdout,"[init_dendritic]:  dj[%d]=%d, dj[%d]=%d, dk[%d][%d]=%d, dk[%d][%d]=%d, dk[%d][%d]=%d\n",
            //   i0,dj[i0],i1,dj[i1],i0,j00,dk[i0][j00],i1,j10,dk[i1][j10],i1,j11,dk[i1][j11]);
            // refine_error = 1;
            // This type of refinement error seems simple enough to correct...
            dk[i0][j00] = MAX(dk[i1][j10],dk[i1][j11]);
          }
        }
      }
    }
  }
  #endif

  #if (NDIM>1)
  if (mpi_io_proc()) fclose(fp);
  #endif

  if (refine_error) exit(1);
  else if (mpi_io_proc()) fprintf(stderr,"No refinement errors detected!\n");

}

/* Count total active + ghost cells on this proc */
uint64_t my_gridcell_count()
{
  int s0,s1,s2;
  int i,j;
  int offset[NDIM];
  static uint64_t cnt = 0;
  //uint64_t cnt = 0;

  if (cnt==0) {
    s0 = my_grid_dims[0] + 2*NG;
    offset[0] = NG - istart[0];
    #if (NDIM==1)
    cnt += s0;
    #else
    for (i=-offset[0]; i<s0-offset[0]; i++) {
      s1 = my_grid_dims[1]/DJS(i) + 2*NG;
      offset[1] = NG - istart[1]/DJS(i);
      #if (NDIM==2)
      cnt += s1;
      // if (myrank==0) printf("i=%d, s1=%d, cnt=%llu\n",i,s1,cnt);
      #else
      for (j=-offset[1]; j<s1-offset[1]; j++) {
        s2 = my_grid_dims[2]/DKS(i,j) + 2*NG;
        offset[2] = NG - istart[2]/DKS(i,j);
        cnt += s2;
      }
      #endif
    }
    #endif
    // printf("[my_gridcell_count]:  myrank=%d, cnt=%llu\n",myrank,cnt);
  }

  return cnt;
}

uint64_t my_gridface_count(int dir)
{
  int s0,s1,s2;
  int i,j;
  int offset[NDIM];
  static uint64_t cnt[3] = {0,0,0};

  if (cnt[dir]==0) {
    s0 = my_grid_dims[0] + 1;
    offset[0] = -istart[0];
    #if (NDIM==1)
    cnt[dir] += s0;
    #else
    for (i=-offset[0]; i<s0-offset[0]; i++) {
      s1 = my_grid_dims[1]/DJS(i) + 1;
      offset[1] = -istart[1]/DJS(i);
      #if (NDIM==2)
      cnt[dir] += s1;
      #else
      for (j=-offset[1]; j<s1-offset[1]; j++) {
        // If using spherical coordinates in 3-d, need to be careful here.
        // There will be coarsening in phi-direction as theta increases near
        // the southern pole, so must check for this case.
        if (dir==1) {
          s2 = my_grid_dims[2]/MIN(DKS(i,j),DKS(i,j-1)) + 1;
          offset[2] = -istart[2]/MIN(DKS(i,j),DKS(i,j-1));
        } else {
          s2 = my_grid_dims[2]/DKS(i,j) + 1;
          offset[2] = -istart[2]/DKS(i,j);
        }
        cnt[dir] += s2;
      }
      #endif
    }
    #endif
    // printf("[my_gridface_count]:  myrank=%d, dir=%d, cnt=%llu\n",myrank,dir,cnt[dir]);
  }

  return cnt[dir];
}

int ND_PTR dendritic_malloc_int()
{
  int s0,s1,s2;
  uint64_t cnt;
  int i,j;
  int offset[NDIM];
  int ND_PTR A = NULL;
  int *tmp = NULL;

  tmp = malloc_rank1(my_gridcell_count(), sizeof *tmp);

  offset[0] = NG - istart[0];
  #if (NDIM==1)
  A = tmp + offset[0];
  #else
  s0 = my_grid_dims[0] + 2*NG;
  A  = malloc_rank1(s0, sizeof *A);
  A += offset[0];
  cnt = 0;
  for (i=-offset[0]; i<s0-offset[0]; i++) {
    s1 = my_grid_dims[1]/DJS(i) + 2*NG;
    offset[1] = NG - istart[1]/DJS(i);
    #if (NDIM==2)
    A[i] = tmp + cnt + offset[1];
    cnt += s1;
    #else
    A[i]  = malloc_rank1(s1, sizeof *A[i]);
    A[i] += offset[1];
    for (j=-offset[1]; j<s1-offset[1]; j++) {
      s2 = my_grid_dims[2]/DKS(i,j) + 2*NG;
      offset[2] = NG - istart[2]/DKS(i,j);
      A[i][j] = tmp + cnt + offset[2];
      cnt += s2;
    }
    #endif
  }
  #endif
  
  return A;
}

int ND_PTR dendritic_malloc_int_square()
{
  int s0=n1+2*NG,s1=n2+2*NG,s2=n3+2*NG;
  uint64_t cnt;
  int i,j,k;
  int offset[NDIM];
  int ND_PTR A = NULL;
  int *tmp = NULL;
 
  #if(NDIM==1)
  tmp = malloc_rank1(s0, sizeof *tmp);
  A = tmp + NG;
  for(i=-NG;i<s0-NG;i++){
    A[i] = -100;
  }
  #elif(NDIM==2)
  tmp = malloc_rank1(s0*s1, sizeof *tmp);
  A = malloc_rank1(s0,sizeof *A);
  A += NG;
  cnt = 0;
  for(i=-NG;i<s0-NG;i++){
    A[i] = tmp+cnt+NG;
    cnt += s1;
  }
  for(i=-NG;i<s0-NG;i++){
    for(j=-NG;j<s1-NG;j++){
      A[i][j] = -100;
    }
  }
  #else
  tmp = malloc_rank1(s0*s1*s2, sizeof *tmp);
  A = malloc_rank1(s0,sizeof *A);
  A += NG;
  cnt = 0;
  for(i=-NG;i<s0-NG;i++){
    A[i] = malloc_rank1(s1, sizeof *A[i]);
    A[i] += NG;
    for(j=-NG;j<s1-NG;j++){
      A[i][j] = tmp + cnt + NG;
      cnt += s2;
    }
  }
  for(i=-NG;i<s0-NG;i++){
    for(j=-NG;j<s1-NG;j++){
      for(k=-NG;k<s2-NG;k++){
        A[i][j][k] = -100;
      }
    }
  }
  #endif

  return A;
}

double ND_PTR dendritic_malloc_double()
{
  int s0,s1,s2;
  uint64_t cnt;
  int i,j;
  int offset[NDIM];
  double* tmp = NULL;
  double ND_PTR A = NULL;

  tmp = malloc_rank1(my_gridcell_count(), sizeof *tmp);
  offset[0] = NG - istart[0];
  #if (NDIM==1)
  A = tmp + offset[0];
  #else
  s0 = my_grid_dims[0] + 2*NG;
  A  = malloc_rank1(s0, sizeof *A);
  A += offset[0];
  cnt = 0;
  for (i=-offset[0]; i<s0-offset[0]; i++) {
    s1 = my_grid_dims[1]/DJS(i) + 2*NG;
    offset[1] = NG - istart[1]/DJS(i);
    #if (NDIM==2)
    A[i] = tmp + cnt + offset[1];
    cnt += s1;
    #else
    A[i]  = malloc_rank1(s1, sizeof *A[i]);
    A[i] += offset[1];
    for (j=-offset[1]; j<s1-offset[1]; j++) {
      s2 = my_grid_dims[2]/DKS(i,j) + 2*NG;
      offset[2] = NG - istart[2]/DKS(i,j);
      A[i][j] = tmp + cnt + offset[2];
      cnt += s2;
    }
    #endif
  }
  #endif

  return A;
}

zone_geom ND_PTR dendritic_malloc_geom()
{
  int s0,s1,s2;
  uint64_t cnt;
  int i,j;
  int offset[NDIM];
  zone_geom* tmp = NULL;
  zone_geom ND_PTR A = NULL;

  tmp = malloc_rank1(my_gridcell_count(), sizeof *tmp);
  offset[0] = NG - istart[0];
  #if (NDIM==1)
  A = tmp + offset[0];
  #else
  s0 = my_grid_dims[0] + 2*NG;
  A  = malloc_rank1(s0, sizeof *A);
  A += offset[0];
  cnt = 0;
  for (i=-offset[0]; i<s0-offset[0]; i++) {
    s1 = my_grid_dims[1]/DJS(i) + 2*NG;
    offset[1] = NG - istart[1]/DJS(i);
    #if (NDIM==2)
    A[i] = tmp + cnt + offset[1];
    cnt += s1;
    #else
    A[i]  = malloc_rank1(s1, sizeof *A[i]);
    A[i] += offset[1];
    for (j=-offset[1]; j<s1-offset[1]; j++) {
      s2 = my_grid_dims[2]/DKS(i,j) + 2*NG;
      offset[2] = NG - istart[2]/DKS(i,j);
      A[i][j] = tmp + cnt + offset[2];
      cnt += s2;
    }
    #endif
  }
  #endif

  return A;
}

double NDP_PTR dendritic_malloc_vec(int vlen)
{
  int s0,s1,s2;
  uint64_t cnt;
  int i,j,k;
  int offset[NDIM];
  double* tmp = NULL;
  double NDP_PTR A = NULL;

  tmp = malloc_rank1(my_gridcell_count()*vlen, sizeof *tmp);
  offset[0] = NG - istart[0];
  s0 = my_grid_dims[0] + 2*NG;
  A  = malloc_rank1(s0, sizeof *A);
  A += offset[0];
  cnt = 0;
  for (i=-offset[0]; i<s0-offset[0]; i++) {
    #if (NDIM==1)
    A[i] = tmp + cnt;
    cnt += vlen;
    #else
    s1 = my_grid_dims[1]/DJS(i) + 2*NG;
    offset[1] = NG - istart[1]/DJS(i);
    A[i]  = malloc_rank1(s1, sizeof *A[i]);
    A[i] += offset[1];
    for (j=-offset[1]; j<s1-offset[1]; j++) {
      #if (NDIM==2)
      A[i][j] = tmp + cnt;
      cnt += vlen;
      #else
      s2 = my_grid_dims[2]/DKS(i,j) + 2*NG;
      offset[2] = NG - istart[2]/DKS(i,j);
      A[i][j]  = malloc_rank1(s2, sizeof *A[i][j]);
      A[i][j] += offset[2];
      for (k=-offset[2]; k<s2-offset[2]; k++) {
        A[i][j][k] = tmp + cnt;
        cnt += vlen;
      }
      #endif
    }
    #endif
  }

  return A;
}

double NDP_PTR dendritic_malloc_vec_face(int vlen, int dir)
{
  int s0,s1,s2;
  uint64_t cnt;
  int i,j,k;
  int offset[NDIM];
  double* tmp = NULL;
  double NDP_PTR A = NULL;

  tmp = malloc_rank1(my_gridface_count(dir)*vlen, sizeof *tmp);
  s0 = my_grid_dims[0] + 1;
  offset[0] = -istart[0];
  A  = malloc_rank1(s0, sizeof *A);
  A += offset[0];
  cnt = 0;
  for (i=-offset[0]; i<s0-offset[0]; i++) {
    #if (NDIM==1)
    A[i] = tmp + cnt;
    cnt += vlen;
    #else
    s1 = my_grid_dims[1]/DJS(i) + 1;
    offset[1] = -istart[1]/DJS(i);
    A[i]  = malloc_rank1(s1, sizeof *A[i]);
    A[i] += offset[1];
    for (j=-offset[1]; j<s1-offset[1]; j++) {
      #if (NDIM==2)
      A[i][j] = tmp + cnt;
      cnt += vlen;
      #else
      // If using spherical coordinates in 3-d, need to be careful here.
      // There will be coarsening in phi-direction as theta increases near
      // the southern pole, so must check for this case.
      if (dir==1) {
        s2 = my_grid_dims[2]/MIN(DKS(i,j),DKS(i,j-1)) + 1;
        offset[2] = -istart[2]/MIN(DKS(i,j),DKS(i,j-1));
      } else {
        s2 = my_grid_dims[2]/DKS(i,j) + 1;
        offset[2] = -istart[2]/DKS(i,j);
      }
      A[i][j]  = malloc_rank1(s2, sizeof *A[i][j]);
      A[i][j] += offset[2];
      for (k=-offset[2]; k<s2-offset[2]; k++) {
        A[i][j][k] = tmp + cnt;
        cnt += vlen;
      }
      #endif
    }
    #endif
  }

  return A;
}

double ** dendritic_malloc_vec_collapsed(int vlen)
{
  double* tmp = NULL;
  double** A = NULL;

  int len = my_gridcell_count();
  tmp = malloc_rank1(my_gridcell_count()*vlen, sizeof *tmp);
  A  = malloc_rank1(my_gridcell_count(), sizeof *A);
  for(int vv=0;vv<my_gridcell_count();vv++){
    A[vv] = tmp+vv*vlen;
  }
  return A;
}

int * dendritic_malloc_int_collapsed()
{
  int* A = NULL;

  A  = malloc_rank1(my_gridcell_count(), sizeof *A);
  return A;
}

double * dendritic_malloc_double_collapsed()
{
  double* A = NULL;

  A  = malloc_rank1(my_gridcell_count(), sizeof *A);
  return A;
}

zone_geom * dendritic_malloc_geom_collapsed()
{
  zone_geom* A = NULL;

  A  = malloc_rank1(my_gridcell_count(), sizeof *A);
  return A;
}

double ** dendritic_malloc_vec_face_collapsed(int vlen, int dir)
{
  int s0,s1,s2;
  int i,j,k;
  double ** A = NULL;
  double *tmp;

  tmp = malloc_rank1(my_gridface_count(dir)*vlen, sizeof *tmp);
  A  = malloc_rank1(my_gridface_count(dir), sizeof *A);
  for(int vv=0;vv<vlen;vv++){
    A[vv] = tmp+vv*my_gridface_count(dir);
  }
  //for(int vv=0;vv<my_gridface_count(dir);vv++){
  //  A[vv] = tmp+vv*vlen;
  //}
  return A;
}

int ND_PTR dendritic_malloc_int_face_square(int dir)
{
  int s0=n1+1,s1=n2+1,s2=n3+1;
  uint64_t cnt;
  int i,j,k;
  int* tmp = NULL;
  int ND_PTR A = NULL;

  #if(NDIM==1)
  tmp = malloc_rank1(s0, sizeof *tmp);
  A = tmp;
  for(i=0;i<s0;i++){
    A[i] = -100;
  }
  #elif(NDIM==2)
  tmp = malloc_rank1(s0*s1, sizeof *tmp);
  A = malloc_rank1(s0,sizeof *A);
  cnt = 0;
  for(i=0;i<s0;i++){
    A[i] = tmp+cnt;
    cnt += s1;
  }
  for(i=0;i<s0;i++){
    for(j=0;j<s1;j++){
      A[i][j] = -100;
    }
  }
  #else
  tmp = malloc_rank1(s0*s1*s2, sizeof *tmp);
  A = malloc_rank1(s0,sizeof *A);
  cnt = 0;
  for(i=0;i<s0;i++){
    A[i] = malloc_rank1(s1, sizeof *A[i]);
    for(j=0;j<s1;j++){
      A[i][j] = tmp + cnt;
      cnt += s2;
    }
  }
  for(i=0;i<s0;i++){
    for(j=0;j<s1;j++){
      for(k=0;k<s2;k++){
        A[i][j][k] = -100;
      }
    }
  }
  #endif
  return A;
}


double** pencil_malloc()
{
  int i;
  double* tmp = NULL;
  double** A = NULL;

  tmp = malloc_rank1(ninterp*max_grid_dim, sizeof *tmp);
  A   = malloc_rank1(ninterp, sizeof *A);
  for (i=0; i<ninterp; i++) {
    A[i] = tmp + i*max_grid_dim;// + NG;
  }

  return A;
}

int get_prims_count()
{
  int ii = 0, jj = 0, kk = 0;
  int vv = 0;
  int count = 0;
  VLOOP {
    ZGLOOP {
      count += 1;
    }
  }
  return count;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void pack_prims(double ** p, double * out)
#else
void pack_prims(double NDP_PTR p, double * out)
#endif
{
  int ii = 0, jj = 0, kk = 0;
  int vv = 0;
  ptrdiff_t i4d = 0;
  ZGLOOP {
    VLOOP {
      out[i4d++] = NDP_ELEM_LINEAR(p, ii, jj, kk, vv);
    }
  }
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void unpack_prims(double const * in, double ** p)
#else
void unpack_prims(double const * in, double NDP_PTR p)
#endif
{
  int ii = 0, jj = 0, kk = 0;
  int vv = 0;
  ptrdiff_t i4d = 0;
  ZGLOOP {
    VLOOP {
      NDP_ELEM_LINEAR(p, ii, jj, kk, vv) = in[i4d++];
    }
  }

}

void calc_cellcenter_maps(){
  int ii,jj,kk,II;
  ZGLOOP{
    ND_ELEM(ijk_to_I,ii,jj,kk)=-1;
  }
  ii=0;
  jj=0;
  kk=0;
  II=0;
  //Bulk
  //Send: bulk cell next to ghost cells

  #if (NDIM==1)
    #if (USE_MPI==TRUE)
    int i,count,blocks;
    proc_info* myproc = &proc[myrank];
    proc_info* nproc = NULL;
  
    for (i=0; i < myproc->nneigh; i++) {
      nproc = &proc[myproc->neigh_id[i]];
      count    = 1;
      blocks   = NG*nvars;
      if (myproc->istart[0] < nproc->istart[0]) {
        for(int isend=0;isend<NG;isend++){
          ND_ELEM(ijk_to_I,istart[0]+my_grid_dims[0]-NG+isend,jj,kk)=II;
          I_to_ijk[3*II]=istart[0]+my_grid_dims[0]-NG+isend;
          I_to_ijk[3*II+1]=0;
          I_to_ijk[3*II+2]=0;
          II++;
        }
      } else {
        for(int isend=0;isend<NG;isend++){
          ND_ELEM(ijk_to_I,istart[0]+isend,jj,kk)=II;
          I_to_ijk[3*II]=istart[0]+isend;
          I_to_ijk[3*II+1]=0;
          I_to_ijk[3*II+2]=0;
          II++;
        }
      }
    }
    #endif
  #endif
  #if (NDIM==2)
    #if (USE_MPI==TRUE)
    int i,j,idir,ip,c,ix,iy;
    int i0,j0,k0;
    int overlap_start,overlap_stop,count,*blocks,*blocks_shock;
    proc_info* myproc = &proc[myrank];
    proc_info* nproc = NULL;
  
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
          i0 = istart[0];
          j0 = istart[1]/DJS(i0);
          k0 = istart[2]/DKS(i0,j0);
          // Send types
          for (c=0; c<count; c++) {
            ix = (idir > 0) ? myproc->istop[0] - NG + c : myproc->istart[0] + c;
            iy = JS(ix,overlap_start);
            for(int isend=0;isend<JS(ix,overlap_stop) - JS(ix,overlap_start);isend++){
              ND_ELEM(ijk_to_I,ix,iy+isend,0)=II;
              I_to_ijk[3*II]=ix;
              I_to_ijk[3*II+1]=iy+isend;
              I_to_ijk[3*II+2]=0;
              II++;
            }
          }
          break;
        case 2: // neighbor is in x1-direction
          // build 1-type, this should be easy
          // figure out overlap in 0-direction
          overlap_start = MAX(myproc->istart[0],nproc->istart[0]);
          overlap_stop  = MIN(myproc->istop[0] ,nproc->istop[0] );
  
          count         = overlap_stop - overlap_start;
          i0 = istart[0];
          j0 = istart[1]/DJS(i0);
          k0 = istart[2]/DKS(i0,j0);
          // Send types
          for (c=0; c<count; c++) {
            // y is the fastest-changing index, hence the need for small blocks...
            ix = overlap_start + c;
            iy = (idir > 0) ? JS(ix,istop[1]) - NG : JS(ix,istart[1]);
            blocks[c] = NG * nvars;
            for(int isend=0;isend<NG;isend++){
              ND_ELEM(ijk_to_I,ix,iy+isend,0)=II;
              I_to_ijk[3*II]=ix;
              I_to_ijk[3*II+1]=iy+isend;
              I_to_ijk[3*II+2]=0;
              II++;
            }
          }
          break;
        default:
          // should never be here
          fprintf(stderr,"Invalid idir in build_types2\n");
          exit(3);
      }  /* end switch */
    }  /* end for */
    #endif  /* USE_MPI==TRUE */
  #endif
  #if (NDIM==3)
    #if (USE_MPI==TRUE)
    int i,idir,ip,ix,iy,iz;
    int i0,j0,k0;
    int overlap_start[NDIM],overlap_stop[NDIM],count,*blocks,c,*blocks_shock;
    proc_info* myproc = &proc[myrank];
    proc_info* nproc = NULL;
  
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
  
          // Send types
          count = 0;
          for (ii=0; ii<NG; ii++) {
            ix = (idir > 0) ? myproc->istop[0] + ii - NG : myproc->istart[0] + ii;
            count += JS(ix,overlap_stop[1]) - JS(ix,overlap_start[1]);
          }
  
          c = 0;
          for (ii=0; ii<NG; ii++) {
            ix = (idir > 0) ? myproc->istop[0] + ii - NG : myproc->istart[0] + ii;
            for (iy = JS(ix,overlap_start[1]); iy < JS(ix,overlap_stop[1]); iy++) {
              iz = KS(ix,iy,overlap_start[2]);
              for(int isend=0;isend<KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2]);isend++){
                ND_ELEM(ijk_to_I,ix,iy,iz+isend)=II;
                I_to_ijk[3*II]=ix;
                I_to_ijk[3*II+1]=iy;
                I_to_ijk[3*II+2]=iz+isend;
                II++;
              }
              c++;
            }
          }
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
          i0 = istart[0];
          j0 = JS(i0,istart[1]);
          k0 = KS(i0,j0,istart[2]);
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
              for(int isend=0;isend<KS(ix,iy,overlap_stop[2]) - KS(ix,iy,overlap_start[2]);isend++){
                ND_ELEM(ijk_to_I,ix,iy,iz+isend)=II;
                I_to_ijk[3*II]=ix;
                I_to_ijk[3*II+1]=iy;
                I_to_ijk[3*II+2]=iz+isend;
                II++;
              }
              c++;
            }
          }
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
          count = 0;
          for (ix = overlap_start[0]; ix < overlap_stop[0]; ix++) {
            count += JS(ix,overlap_stop[1]) - JS(ix,overlap_start[1]);
          }
          // Send types
          c = 0;
          for (ix = overlap_start[0]; ix < overlap_stop[0]; ix++) {
            for (iy = JS(ix,overlap_start[1]); iy < JS(ix,overlap_stop[1]); iy++) {
              iz = (idir > 0) ? KS(ix,iy,istop[2]) - NG : KS(ix,iy,istart[2]);
              blocks[c] = NG * nvars;
              for(int isend=0;isend<NG;isend++){
                ND_ELEM(ijk_to_I,ix,iy,iz+isend)=II;
                I_to_ijk[3*II]=ix;
                I_to_ijk[3*II+1]=iy;
                I_to_ijk[3*II+2]=iz+isend;
                II++;
              }
              c++;
            }
          }
          break;
  
        default:
          // should never be here
          fprintf(stderr,"Invalid idir in build_types2\n");
          exit(3);
      }
    }
    #endif /* USE_MPI==TRUE */
  #endif /* NDIM==3 */

  cell_count_send = II;
  ii=0;
  jj=0;
  kk=0;
  //Other bulk cells
  ZLOOP {
    if(ND_ELEM(ijk_to_I,ii,jj,kk)<0){
      ND_ELEM(ijk_to_I,ii,jj,kk)=II;
      I_to_ijk[3*II]=ii;
      I_to_ijk[3*II+1]=jj;
      I_to_ijk[3*II+2]=kk;
      II++;
    }
  }
  cell_count = II;
  //Recv: ghost cells
  ii=0;
  jj=0;
  kk=0;
  ZGLOOP {
    if(ND_ELEM(ijk_to_I,ii,jj,kk)<0){
      ND_ELEM(ijk_to_I,ii,jj,kk)=II;
      I_to_ijk[3*II]=ii;
      I_to_ijk[3*II+1]=jj;
      I_to_ijk[3*II+2]=kk;
      II++;
    }
  }
  if(II!=my_gridcell_count()){printf("Error: Number of cells wrong!\n");exit(3);}
  cell_count_all = II;
  cell_count_recv = cell_count_all-cell_count;


  #if(NDIM==1)
  GPU_PRAGMA(omp target enter data map(to:ijk_to_I[-NG:n1+2*NG]))
  #elif(NDIM==2)
  GPU_PRAGMA(omp target enter data map(to:ijk_to_I[-NG:n1+2*NG]))
  for(int itemp=-NG;itemp<n1+NG;itemp++){
    GPU_PRAGMA(omp target enter data map(to:ijk_to_I[itemp][-NG:n2+2*NG]))
  }
  #else
  GPU_PRAGMA(omp target enter data map(to:ijk_to_I[-NG:n1+2*NG]))
  for(int itemp=-NG;itemp<n1+NG;itemp++){
    GPU_PRAGMA(omp target enter data map(to:ijk_to_I[itemp][-NG:n2+2*NG]))
    for(int jtemp=-NG;jtemp<n2+NG;jtemp++){
      GPU_PRAGMA(omp target enter data map(to:ijk_to_I[itemp][jtemp][-NG:n3+2*NG]))
    }
  }
  //GPU_PRAGMA(omp target enter data map(to:ijk_to_I[-NG:n1+2*NG][-NG:n2+2*NG][-NG:n3+2*NG]))
  #endif
  GPU_PRAGMA(omp target enter data map(to:I_to_ijk[:3*cell_count_all]))

  int itest,jtest,ktest;
  int error=0;
  ZGLOOP{
    II = ND_ELEM(ijk_to_I,ii,jj,kk);
    GET_IJK_FROM_I(II,itest,jtest,ktest);
    if(ii!=itest||jj!=jtest||kk!=ktest){error=1;}
  }
  if(error==1){printf("Collapsed Indices Wrong!\n");exit(1);}
}

void calc_facecenter_maps(){
  int s0,s1,s2;
  int offset[NDIM];
  int II,dd;
  s0 = my_grid_dims[0] + 1;
  offset[0] = -istart[0];
  DLOOP {
    II = 0;
    for (int i=-offset[0]; i<s0-offset[0]; i++) {
      #if (NDIM==1)
      ND_ELEM(face_ijk_to_I[dd],i,j,k)=II;
      II += 1;
      #else
      s1 = my_grid_dims[1]/DJS(i) + 1;
      offset[1] = -istart[1]/DJS(i);
      for (int j=-offset[1]; j<s1-offset[1]; j++) {
        #if (NDIM==2)
        ND_ELEM(face_ijk_to_I[dd],i,j,k)=II;
        II += 1;
        #else
        // If using spherical coordinates in 3-d, need to be careful here.
        // There will be coarsening in phi-direction as theta increases near
        // the southern pole, so must check for this case.
        if (dd==1) {
          s2 = my_grid_dims[2]/MIN(DKS(i,j),DKS(i,j-1)) + 1;
          offset[2] = -istart[2]/MIN(DKS(i,j),DKS(i,j-1));
        } else {
          s2 = my_grid_dims[2]/DKS(i,j) + 1;
          offset[2] = -istart[2]/DKS(i,j);
        }
        for (int k=-offset[2]; k<s2-offset[2]; k++) {
          ND_ELEM(face_ijk_to_I[dd],i,j,k)=II;
          II += 1;
        }
        #endif
      }
      #endif
    }
  }

	if(myrank==1){
  printf("%d %d %d\n",face_ijk_to_I[0][516],-offset[0],s0-offset[0]);
	}

  s0=n1+1;
  s1=n2+1;
  s2=n3+1;
  #if(NDIM==1)
  GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[:3]))
  for(int dirtemp=0;dirtemp<3;dirtemp++){
    GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[dirtemp][:s0]))
  }
  #elif(NDIM==2)
  GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[:3]))
  for(int dirtemp=0;dirtemp<3;dirtemp++){
    GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[dirtemp][:s0]))
    for(int itemp=0;itemp<s0;itemp++){
      GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[dirtemp][itemp][:s1]))
    }
  }
  #else
  GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[:3]))
  for(int dirtemp=0;dirtemp<3;dirtemp++){
    GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[dirtemp][:s0]))
    for(int itemp=0;itemp<s0;itemp++){
      GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[dirtemp][itemp][:s1]))
      for(int jtemp=0;jtemp<s1;jtemp++){
        GPU_PRAGMA(omp target enter data map(to:face_ijk_to_I[dirtemp][itemp][jtemp][:s2]))
      }
    }
  }
  #endif

}


void allocate_memory()
{
	int dd,vv;
	int pmemsize[NDIM+1],fmemsize[NDIM+1];
	int ememsize[NDIM+1],gmemsize[NDIM];
	int pencil_memsize[2];
	int pbuf[NDIM], fbuf[NDIM], gbuf[NDIM];
	int pencil_buf[2];
	int *tmp=NULL;

	if (mpi_io_proc()) {
		fprintf(stderr,"Entering allocate_memory...");
		fflush(stderr);
	}

	bc = malloc_rank1(nvars, sizeof *bc);
  I_to_ijk = malloc(my_gridcell_count()*3*(sizeof *I_to_ijk));
  ijk_to_I = dendritic_malloc_int_square();

  calc_cellcenter_maps();
		fprintf(stderr,"cell centered maps done...");
  //test
//  GPU_PRAGMA(omp target teams distribute parallel for schedule(static,1))
//  for(II=0;II<cell_count_all;II++){
//    GET_IJK_FROM_I(II,ii,jj,kk);
//    printf("%d %d %d %d %d\n",II,ii,jj,kk,omp_get_thread_num());
//  }

  face_ijk_to_I = malloc_rank1(3,sizeof *face_ijk_to_I);
  face_ijk_to_I[0] = dendritic_malloc_int_face_square(0);
  face_ijk_to_I[1] = dendritic_malloc_int_face_square(1);
  face_ijk_to_I[2] = dendritic_malloc_int_face_square(2);
  face_count_0 = my_gridface_count(0);
  face_count_1 = my_gridface_count(1);
  face_count_2 = my_gridface_count(2);
  calc_facecenter_maps();
		fprintf(stderr,"face centered maps done...");


  #if (USE_LINEAR_ALLOCATION==TRUE)
  sim_p   = dendritic_malloc_vec_collapsed(nvars);
  sim_ph  = dendritic_malloc_vec_collapsed(nvars);
  sim_eos = dendritic_malloc_vec_collapsed(NEOS);
  sim_src = dendritic_malloc_vec_collapsed(nvars);
  // shock_flag
  sim_shock_flag = dendritic_malloc_int_collapsed();
  // gravitational potential
  sim_Phi = dendritic_malloc_double_collapsed();
  // geometry on the grid
  geom = dendritic_malloc_geom_collapsed();
#if (USE_LARGER_STEP==TRUE)
  sim_dudt   = dendritic_malloc_vec_collapsed(nvars);
  sim_deosdt = dendritic_malloc_vec_collapsed(NEOS);
#endif
  // data on faces
  DLOOP {
    if(dd==0) sim_fdir0 = dendritic_malloc_vec_face_collapsed(nvars,dd);
    if(dd==1) sim_fdir1 = dendritic_malloc_vec_face_collapsed(nvars,dd);
    if(dd==2) sim_fdir2 = dendritic_malloc_vec_face_collapsed(nvars,dd);
    //sim_f[dd]     = dendritic_malloc_vec_face(nvars,dd);
    if(dd==0) sim_vedgedir0 = dendritic_malloc_vec_face_collapsed(SPACEDIM,dd);
    if(dd==1) sim_vedgedir1 = dendritic_malloc_vec_face_collapsed(SPACEDIM,dd);
    if(dd==2) sim_vedgedir2 = dendritic_malloc_vec_face_collapsed(SPACEDIM,dd);
  }
  #else
  sim_p   = dendritic_malloc_vec(nvars);
  sim_ph  = dendritic_malloc_vec(nvars);
  sim_eos = dendritic_malloc_vec(NEOS);
  sim_src = dendritic_malloc_vec(nvars);
  // shock_flag
  sim_shock_flag = dendritic_malloc_int();
  // gravitational potential
  sim_Phi = dendritic_malloc_double();
  // geometry on the grid
  geom = dendritic_malloc_geom();
#if (USE_LARGER_STEP==TRUE)
  sim_dudt   = dendritic_malloc_vec(nvars);
  sim_deosdt = dendritic_malloc_vec(NEOS);
#endif
  // data on faces
  DLOOP {
    if(dd==0) sim_fdir0 = dendritic_malloc_vec_face(nvars,dd);
    if(dd==1) sim_fdir1 = dendritic_malloc_vec_face(nvars,dd);
    if(dd==2) sim_fdir2 = dendritic_malloc_vec_face(nvars,dd);
    //sim_f[dd]     = dendritic_malloc_vec_face(nvars,dd);
    if(dd==0) sim_vedgedir0 = dendritic_malloc_vec_face(SPACEDIM,dd);
    if(dd==1) sim_vedgedir1 = dendritic_malloc_vec_face(SPACEDIM,dd);
    if(dd==2) sim_vedgedir2 = dendritic_malloc_vec_face(SPACEDIM,dd);
  }
  #endif

  #if (GR_MONOPOLE==TRUE)
  // lapse and lapse_edge
  gr_lapse_edge = malloc_rank1(n1+1, sizeof *gr_lapse_edge);
  gr_lapse      = malloc_rank1(n1,   sizeof *gr_lapse     );
  #endif

  // rad_fluid.c
  Trad = malloc_rank1(ngroups*SQR(SPACEDIM), sizeof *Trad);

  max_grid_dim=0;
  DLOOP {
    pmemsize[dd] = my_grid_dims[dd] + 2*NG;
    if (pmemsize[dd] > max_grid_dim) max_grid_dim = pmemsize[dd];
  }
  
  {
    //pencil = pencil_malloc();
    //pleft  = pencil_malloc();
    //pright = pencil_malloc();
    //flatten = malloc_rank1(max_grid_dim, sizeof *flatten);
    ////flatten += NG;
  }

  alpha0  = malloc_rank1(n1+2*NG, sizeof *alpha0);
  alpha0 += NG;
  beta0   = malloc_rank1(n1+2*NG, sizeof *beta0 );
  beta0  += NG;
  Gamma0  = malloc_rank1(n1+2*NG, sizeof *Gamma0);
  Gamma0 += NG;

  interp_order = malloc_rank1(ninterp, sizeof *interp_order);
  
  if (mpi_io_proc()) {
  	fprintf(stderr,"done!\n");
  	fflush(stderr);
  }

  
  return;
}

#if (NDIM==3 && GEOM==SPHERICAL)
static double * u    = NULL;
static double * uavg = NULL;

// EXPERIMENTAL!!!
/* Average the solution over all polar zones.  Cartesian components of
 * vectors are averaged.  This code assumes a single processor owns all
 * the zones around the pole.
 */
#if (USE_LINEAR_ALLOCATION==TRUE)
void avg_poles_3d(double ** p)
#else
void avg_poles_3d(double NDP_PTR p)
#endif
{
  static int firstc = 1;
  int l,m,g;
  double lam[SPACEDIM][SPACEDIM];
  double r,th,ph,sth,cth,sph,cph;
  int s,lev;
  int ii,jj,kk,dd,vv;
  zone_geom *gm;
  double vol_sum,vol_frac;
  int pole;
  short do_north_pole, do_south_pole;

  TIMER_START("avg_poles_3d");

  if(firstc) {
    {
      u    = malloc_rank1(nvars, sizeof(*u));
      uavg = malloc_rank1(nvars, sizeof(*uavg));
    }
    firstc = 0;
  }

  for(int pole=0; pole<2; pole++) {
    // North pole
    if(pole == 0) {
      if(istart[1]==0) {
        do_north_pole = 1;
        do_south_pole = 0;
      }
      else {
        continue;
      }
    }
    if(pole == 1) {
      if(istop[1]==global_grid_dims[1]) {
        do_north_pole = 0;
        do_south_pole = 1;
      }
      else {
        continue;
      }
    }
    ISLOOP(ii) {
      // North pole
      if(do_north_pole) {
        jj = JS(ii,istart[1]);
      }
      // South pole
      else if(do_south_pole) {
        jj = JS(ii,istop[1]) - 1;
      }
      else {
        fprintf(stdout, "Fatal bug in avg_poles_3d\n");
        fflush(stdout);
        abort(); // this is a bug
      }

      lev = 0;
      s = DJS(ii);
      while (s >>= 1) lev++;
      r = r_of_x(rx_info,beta0[ii]/Gamma0[ii]);

      th = th_of_x(thx_info,beta1s[lev][jj]/Gamma1s[lev][jj]);
      sth = sin(th);
      cth = cos(th);

      vol_sum = 0.0;
      KSLOOP(ii,jj,kk) vol_sum += ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
      vol_sum = 1.0/vol_sum;

      memset(uavg, 0, nvars*sizeof(*uavg));

      KSLOOP(ii,jj,kk) {
        ph = startx[2] + (kk+0.5)*DKS(ii,jj)*dx[2];
        sph = sin(ph);
        cph = cos(ph);
        vol_frac = ND_ELEM_LINEAR(geom,ii,jj,kk).volume*vol_sum;

        // Compute conservative variables
        gm = &ND_ELEM_LINEAR(geom,ii,jj,kk);
        p_to_u(ND_ELEM_LINEAR(p,ii,jj,kk),u,gm,nhydro);

        // Compute transformation matrix at zone volume centers
        lam[0][0] = sth*cph;
        lam[0][1] = cth*cph;
        lam[0][2] = -sph;

        lam[1][0] = sth*sph;
        lam[1][1] = cth*sph;
        lam[1][2] = cph;

        lam[2][0] = cth;
        lam[2][1] = -sth;
        lam[2][2] = 0.0;

        // Compute volume-weighted mean of density and energy around the pole
        uavg[URHO] += u[URHO]*vol_frac;
        uavg[ETOT] += u[ETOT]*vol_frac;
        for(dd=U1+SPACEDIM; dd<nhydro; dd++) {
          uavg[dd] += u[dd]*vol_frac;
        }

        #if (DO_RADIATION)
        GLOOP {
          uavg[irad1+g] += u[irad1+g]*vol_frac;
        }
        #endif

        for (l=0; l<3; l++) {
          for (m=0; m<3; m++) {
            // Compute volume-weighted mean Cartesian components of momenta
            uavg[U1+l] += lam[l][m] * u[U1+m] * gm->scale[0][m]*gm->gcon[m] * vol_frac;
            #if (DO_RADIATION)
            GLOOP {
              uavg[ifrad1+g*NDIM+l] += lam[l][m] * u[ifrad1+g*NDIM+m] *
                gm->scale[0][m]*gm->gcon[m] * vol_frac;
            }
            #endif
          }
        }
      }

      KSLOOP(ii,jj,kk) {
        ph = startx[2] + (kk+0.5)*DKS(ii,jj)*dx[2];
        sph = sin(ph);
        cph = cos(ph);

        gm = &ND_ELEM_LINEAR(geom,ii,jj,kk);

        lam[0][0] = sth*cph;
        lam[0][1] = cth*cph;
        lam[0][2] = -sph;

        lam[1][0] = sth*sph;
        lam[1][1] = cth*sph;
        lam[1][2] = cph;

        lam[2][0] = cth;
        lam[2][1] = -sth;
        lam[2][2] = 0.0;

        u[URHO] = uavg[URHO];
        u[ETOT] = uavg[ETOT];
        for(dd=U1+SPACEDIM; dd<nhydro; dd++) {
          u[dd] = uavg[dd];
        }

        #if (DO_RADIATION)
        GLOOP {
          u[irad1+g] = uavg[irad1+g];
        }
        #endif

        for (l=0; l<3; l++) {
          u[U1+l] = 0;
          for (m=0; m<3; m++) {
            // The transformation matrix is orthonormal, so its inverse is its transpose
            u[U1+l] += lam[m][l] * uavg[U1+m] * gm->gcov[0][l]/gm->scale[0][l];
          }
          #if (DO_RADIATION)
          GLOOP {
            u[ifrad1+g*NDIM+l] = 0.0;
            for(m=0; m<3; m++) {
              u[irad1+g*NDIM+l] += lam[m][l] * uavg[irad1+g*NDIM+m] *
                gm->gcov[0][l]/gm->scale[0][l];
            }
          }
          #endif
        }

        u_to_p(u, ND_ELEM_LINEAR(p,ii,jj,kk), gm,nhydro,rho_floor,e_floor);
      }
    }
  }

  TIMER_STOP;
  return;
}
#endif /* NDIM==3 */

/*------------------------------------------------------------------------------
 *  Function ran2
 *
 *  The routine ran2() is extracted from the Numerical Recipes in C
 *  (version 2) code.  I've modified it to use doubles instead of
 *  floats. -- T. A. Gardiner -- Aug. 12, 2003
 *
 *  Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 *  with Bays-Durham shuffle and added safeguards.  Returns a uniform
 *  random deviate between 0.0 and 1.0 (exclusive of the endpoint
 *  values).  Call with idum = a negative integer to initialize;
 *  thereafter, do not alter idum between successive deviates in a
 *  sequence.  RNMX should appriximate the largest floating point
 *  value that is less than 1.
 */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)
#define NTAB 32

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

#undef OFST
#undef KCOMP
#undef KWVM

#if (USE_LINEAR_ALLOCATION==TRUE)
int check_prim(double ** p, int ii, int jj, int kk, char *msg)
#else
int check_prim(double NDP_PTR p, int ii, int jj, int kk, char *msg)
#endif
{
  int vv;
  int status = 0;

  VLOOP {
    if (isnan(NDP_ELEM_LINEAR(p,ii,jj,kk,vv))) {
      printf("[check_prim]:  %s, proc %d detected NaN in variable %d at zone (%d,%d,%d)\n",msg,myrank,vv,ii,jj,kk);
      status = 1;
    }
    if (isinf(NDP_ELEM_LINEAR(p,ii,jj,kk,vv))) {
      printf("[check_prim]:  %s, proc %d detected Inf in variable %d at zone (%d,%d,%d)\n",msg,myrank,vv,ii,jj,kk);
      status = 1;
    }
    if (vv==RHO && NDP_ELEM_LINEAR(p,ii,jj,kk,vv) <= 0.0) {
      printf("[check_prim]:  %s, proc %d detected RHO=%e at zone (%d,%d,%d)\n",msg,myrank,NDP_ELEM_LINEAR(p,ii,jj,kk,vv),ii,jj,kk);
      status = 1;
    }
    if (vv==UU && NDP_ELEM_LINEAR(p,ii,jj,kk,vv) <= 0.0) {
      printf("[check_prim]:  %s, proc %d detected UU=%e at zone (%d,%d,%d)\n",msg,myrank,NDP_ELEM_LINEAR(p,ii,jj,kk,vv),ii,jj,kk);
      status = 1;
    }
    #if (0 && DO_RADIATION)
    if (vv>=irad1 && vv<irad1+ngroups && NDP_ELEM_LINEAR(p,ii,jj,kk,vv) < 0.0) {
      printf("[check_prim]:  %s, proc %d detected Er[%d]=%e at zone (%d,%d,%d)\n",msg,myrank,vv-irad1,NDP_ELEM_LINEAR(p,ii,jj,kk,vv),ii,jj,kk);
      status = 0;
    }
    #endif
  }

  return status;
}

int check_cons(double *u, char *msg, int ii, int jj, int kk)
{
  int vv;
  int status = 0;

  VLOOP {
    if (isnan(u[vv])) {
      printf("[check_cons]:  %s, proc %d detected NaN in variable %d at zone (%d,%d,%d)\n",msg,myrank,vv,ii,jj,kk);
      status = 1;
    }
    if (isinf(u[vv])) {
      printf("[check_cons]:  %s, proc %d detected Inf in variable %d at zone (%d,%d,%d)\n",msg,myrank,vv,ii,jj,kk);
      status = 1;
    }
    if (vv==URHO && u[vv] <= 0.0) {
      printf("[check_cons]:  %s, proc %d detected URHO=%e at zone (%d,%d,%d)\n",msg,myrank,u[vv],ii,jj,kk);
      status = 1;
    }
    if (vv==ETOT && u[vv] <= 0.0) {
      printf("[check_cons]:  %s, proc %d detected ETOT=%e at zone (%d,%d,%d)\n",msg,myrank,u[vv],ii,jj,kk);
      status = 1;
    }
    #if (0 && DO_RADIATION)
    if (vv>=irad1 && vv<irad1+ngroups && u[vv] < 0.0) {
      printf("[check_cons]:  %s, proc %d detected Er[%d]=%e at zone (%d,%d,%d)\n",msg,myrank,vv-irad1,u[vv],ii,jj,kk);
      status = 1;
    }
    #endif
  }

  return status;
}

// check for Inf and NaN, negative or zero densities or energies
#if (USE_LINEAR_ALLOCATION==TRUE)
int find_prim_errors(double ** p, int check_ghost, char *msg)
#else
int find_prim_errors(double NDP_PTR p, int check_ghost, char *msg)
#endif
{
  int ii,jj=0,kk=0,jstart,jstop,kstart,kstop;
  int status = 0;
  double var;

  ZLOOP {
    status |= check_prim(p,ii,jj,kk,msg);
  }

  if (check_ghost) {
    // Inner i-boundary
    for (ii=istart[0]-NG; ii<istart[0]; ii++) {
      #if (NDIM>1)
      JSLOOP(ii,jj) {
      #endif
        #if (NDIM==3)
        KSLOOP(ii,jj,kk) {
        #endif
          status |= check_prim(p,ii,jj,kk,msg);
        #if (NDIM==3)
        }
        #endif
      #if (NDIM>1)
      }
      #endif
    }

    // Outer i-boundary
    for (ii=istop[0]; ii<istop[0]+NG; ii++) {
      #if (NDIM>1)
      JSLOOP(ii,jj) {
      #endif
        #if (NDIM==3)
        KSLOOP(ii,jj,kk) {
        #endif
          status |= check_prim(p,ii,jj,kk,msg);
        #if (NDIM==3)
        }
        #endif
      #if (NDIM>1)
      }
      #endif
    }

    #if (NDIM>1)
    // Inner j-boundary
    ISLOOP(ii) {
      jstart = JS(ii,istart[1]);
      for (jj=jstart-NG; jj<jstart; jj++) {
        #if (NDIM==3)
        KSLOOP(ii,jj,kk) {
        #endif
          status |= check_prim(p,ii,jj,kk,msg);
        #if (NDIM==3)
        }
        #endif
      }
    }

    // Outer j-boundary
    ISLOOP(ii) {
      jstop = JS(ii,istop[1]);
      for (jj=jstop; jj<jstop+NG; jj++) {
        #if (NDIM==3)
        KSLOOP(ii,jj,kk) {
        #endif
          status |= check_prim(p,ii,jj,kk,msg);
        #if (NDIM==3)
        }
        #endif
      }
    }
    #endif

    #if (NDIM==3)
    // Inner k-boundary
    ISLOOP(ii) {
      JSLOOP(ii,jj) {
        kstart = KS(ii,jj,istart[2]);
        for (kk=kstart-NG; kk<kstart; kk++) {
          status |= check_prim(p,ii,jj,kk,msg);
        }
      }
    }

    // Outer k-boundary
    ISLOOP(ii) {
      JSLOOP(ii,jj) {
        kstop = KS(ii,jj,istop[2]);
        for (kk=kstop; kk<kstop+NG; kk++) {
          status |= check_prim(p,ii,jj,kk,msg);
        }
      }
    }
    #endif
  }

  return status;
}

double check_flux_differencing(int verbose, int check_shells)
{
  int ii,jj=0,kk=0,ip,jp,kp,njp,nkp,vv,j,k;
  static double* fluxsum0a = NULL;
  static double* fluxsum0b = NULL;
  static double* fluxsum1a = NULL;
  static double* fluxsum1b = NULL;
  static double* fluxsum2a = NULL;
  static double* fluxsum2b = NULL;
  static double* fluxmax0 = NULL;
  static double* fluxmax1 = NULL;
  static double* fluxmax2 = NULL;
  double fact,fluxin,fluxout;
  double relflux0,relflux1,relflux2;
  double relfluxmax = 0.0;
  double shellrelfluxmax0 = 0.0;
  double shellrelfluxmax1 = 0.0;
  double shellrelfluxmax2 = 0.0;
  double shellrelfluxmax = 0.0;
  static FILE* fp = NULL;
  static FILE* fps = NULL;
  static int firstc = 1;
  static double* shellfluxsum0a = NULL;
  static double* shellfluxsum0b = NULL;
  static double* shellfluxsum1a = NULL;
  static double* shellfluxsum1b = NULL;
  static double* shellfluxsum2a = NULL;
  static double* shellfluxsum2b = NULL;
  static double* shellfluxmax0 = NULL;
  static double* shellfluxmax1 = NULL;
  static double* shellfluxmax2 = NULL;


  if (firstc) {
    if (mpi_io_proc()) {
      fp = fopen("flux_differencing.txt","w");
      if (check_shells) fps = fopen("shell_flux_differencing.txt","w");
    }
    firstc = 0;
  }

  fluxsum0a = malloc_rank1(nvars, sizeof *fluxsum0a);
  fluxsum0b = malloc_rank1(nvars, sizeof *fluxsum0b);
  fluxsum1a = malloc_rank1(nvars, sizeof *fluxsum1a);
  fluxsum1b = malloc_rank1(nvars, sizeof *fluxsum1b);
  fluxsum2a = malloc_rank1(nvars, sizeof *fluxsum2a);
  fluxsum2b = malloc_rank1(nvars, sizeof *fluxsum2b);
  fluxmax0 = malloc_rank1(nvars, sizeof *fluxmax0);
  fluxmax1 = malloc_rank1(nvars, sizeof *fluxmax1);
  fluxmax2 = malloc_rank1(nvars, sizeof *fluxmax2);
  VLOOP {
    fluxsum0a[vv] = 0.0;
    fluxsum0b[vv] = 0.0;
    fluxsum1a[vv] = 0.0;
    fluxsum1b[vv] = 0.0;
    fluxsum2a[vv] = 0.0;
    fluxsum2b[vv] = 0.0;
    fluxmax0[vv] = 0.0;
    fluxmax1[vv] = 0.0;
    fluxmax2[vv] = 0.0;
  }

  if (check_shells) {
    shellfluxsum0a = malloc_rank1(n1, sizeof *shellfluxsum0a);
    shellfluxsum0b = malloc_rank1(n1, sizeof *shellfluxsum0b);
    shellfluxsum1a = malloc_rank1(n1, sizeof *shellfluxsum1a);
    shellfluxsum1b = malloc_rank1(n1, sizeof *shellfluxsum1b);
    shellfluxsum2a = malloc_rank1(n1, sizeof *shellfluxsum2a);
    shellfluxsum2b = malloc_rank1(n1, sizeof *shellfluxsum2b);
    shellfluxmax0 = malloc_rank1(n1, sizeof *shellfluxmax0);
    shellfluxmax1 = malloc_rank1(n1, sizeof *shellfluxmax1);
    shellfluxmax2 = malloc_rank1(n1, sizeof *shellfluxmax2);
    for (ii=0; ii<n1; ii++) {
      shellfluxsum0a[ii] = 0.0;
      shellfluxsum0b[ii] = 0.0;
      shellfluxsum1a[ii] = 0.0;
      shellfluxsum1b[ii] = 0.0;
      shellfluxsum2a[ii] = 0.0;
      shellfluxsum2b[ii] = 0.0;
      shellfluxmax0[ii] = 0.0;
      shellfluxmax1[ii] = 0.0;
      shellfluxmax2[ii] = 0.0;
    }
  }

  VLOOP {
    ZLOOP {
      /* Flux differences in 0-direction */
      fact = (ii==global_grid_dims[0]-1) ? 0.0 : 1.0;
      jp = jj*DJS(ii)/DJS(ii+1);
      njp = (DJS(ii+1) < DJS(ii)) ? DJS(ii)/DJS(ii+1) : 1;
      for (j=jp; j<jp+njp; j++) {
        nkp = DKS(ii,jj)/DKS(ii+1,j);
        kp = kk*nkp;
        for (k=kp; k<kp+nkp; k++) {
          fluxout = fact*NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,j,k,vv);
          fluxsum0a[vv] += fluxout;
          fluxmax0[vv] = MAX(fluxmax0[vv],fabs(fluxout));
          if (check_shells && ii<n1-1 && vv==URHO) {
            shellfluxsum0a[ii+1] += fluxout;
            shellfluxmax0[ii+1] = MAX(shellfluxmax0[ii+1],fabs(fluxout));
          }
        }
      }
      fluxin = NDP_ELEM_LINEAR_F(sim_fdir0,0,ii,jj,kk,vv);
      fluxsum0b[vv] -= fluxin;
      fluxmax0[vv] = MAX(fluxmax0[vv],fabs(fluxin));
      if (check_shells && vv==URHO) {
        shellfluxsum0b[ii] -= fluxin;
        shellfluxmax0[ii] = MAX(shellfluxmax0[ii],fabs(fluxin));
      }

      #if (NDIM>1)
      /* Flux differences in 1-direction */
      /* Outer 1-face */
      if (DKS(ii,jj+1) > DKS(ii,jj)) {
        // Coarsening boundary, 1 outer neighbor, shifted area index, half area, original flux index
        fluxout = 0.5*NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,kk,vv);
      } else if (DKS(ii,jj+1) < DKS(ii,jj)) {
        // Refinement boundary, 2 outer neighbors, shifted area index, original area, shifted flux index
        fluxout = NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,2*kk,vv) + NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,2*kk+1,vv);
      } else {
        // Regular boundary, 1 outer neighbor, original area index, original area, original flux index
        fluxout = NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,kk,vv);
      }
      fluxsum1a[vv] += fluxout;
      fluxmax1[vv] = MAX(fluxmax1[vv],fabs(fluxout));
      if (check_shells && vv==URHO) {
        shellfluxsum1a[ii] += fluxout;
        shellfluxmax1[ii] = MAX(shellfluxmax1[ii],fabs(fluxout));
      }

      /* Inner 1-face */
      if (DKS(ii,jj-1) < DKS(ii,jj)) {
        // Coarsening boundary, 2 inner neighbors, original area index, half area, shifted flux index
        fluxin = 0.5*NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,2*kk,vv) + 0.5*NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,2*kk+1,vv);
      } else {
        // Refinement boundary, 1 inner neighbor, original area index, original area, original flux index
        // Regular boundary, 1 inner neighbor, original area index, original area, original flux index
        fluxin = NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,kk,vv);
      }
      fluxsum1b[vv] -= fluxin;
      fluxmax1[vv] = MAX(fluxmax1[vv],fabs(fluxin));
      if (check_shells && vv==URHO) {
        shellfluxsum1b[ii] -= fluxin;
        shellfluxmax1[ii] = MAX(shellfluxmax1[ii],fabs(fluxin));
      }
      #endif

      #if (NDIM>2)
      /* Flux differences in 2-direction */
      fluxout = NDP_ELEM_LINEAR_F(sim_fdir2,2,ii,jj,kk+1,vv);
      fluxin  = NDP_ELEM_LINEAR_F(sim_fdir2,2,ii,jj,kk  ,vv);
      fluxsum2a[vv] += fluxout;
      fluxsum2b[vv] -= fluxin;
      fluxmax2[vv] = MAX(fluxmax2[vv],MAX(fabs(fluxout),fabs(fluxin)));
      if (check_shells && vv==URHO) {
        shellfluxsum2a[ii] += fluxout;
        shellfluxsum2b[ii] -= fluxin;
        shellfluxmax2[ii] = MAX(shellfluxmax2[ii],MAX(fabs(fluxout),fabs(fluxin)));
      }
      #endif
    }
  }

  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE, fluxsum0a, nvars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxsum0b, nvars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxsum1a, nvars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxsum1b, nvars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxsum2a, nvars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxsum2b, nvars, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxmax0, nvars, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxmax1, nvars, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, fluxmax2, nvars, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (check_shells) {
    MPI_Allreduce(MPI_IN_PLACE, shellfluxsum0a, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxsum0b, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxsum1a, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxsum1b, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxsum2a, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxsum2b, n1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxmax0, n1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxmax1, n1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, shellfluxmax2, n1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  #endif

  if (mpi_io_proc()) {
    VLOOP {
      double fluxmax = MAX(fluxmax0[vv],MAX(fluxmax1[vv],fluxmax2[vv]));
      relflux0 = fabs(fluxsum0a[vv]+fluxsum0b[vv])/(fluxmax + 1.0e-50);
      relflux1 = fabs(fluxsum1a[vv]+fluxsum1b[vv])/(fluxmax + 1.0e-50);
      relflux2 = fabs(fluxsum2a[vv]+fluxsum2b[vv])/(fluxmax + 1.0e-50);
      //relflux0 = fabs(fluxsum0a[vv]+fluxsum0b[vv])/(fluxmax0[vv] + 1.0e-50);
      //relflux1 = fabs(fluxsum1a[vv]+fluxsum1b[vv])/(fluxmax1[vv] + 1.0e-50);
      //relflux2 = fabs(fluxsum2a[vv]+fluxsum2b[vv])/(fluxmax2[vv] + 1.0e-50);
      relfluxmax = MAX(relfluxmax,MAX(relflux0,MAX(relflux1,relflux2)));
    }
    if (check_shells) {
      for (ii=0; ii<n1; ii++) {
        double fluxmax = MAX(shellfluxmax0[ii],MAX(shellfluxmax1[ii],shellfluxmax2[ii]));
        relflux0 = fabs(shellfluxsum0a[ii]+shellfluxsum0b[ii])/(fluxmax + 1.0e-50);
        relflux1 = fabs(shellfluxsum1a[ii]+shellfluxsum1b[ii])/(fluxmax + 1.0e-50);
        relflux2 = fabs(shellfluxsum2a[ii]+shellfluxsum2b[ii])/(fluxmax + 1.0e-50);
        //relflux0 = fabs(shellfluxsum0a[ii]+shellfluxsum0b[ii])/(shellfluxmax0[ii] + 1.0e-50);
        //relflux1 = fabs(shellfluxsum1a[ii]+shellfluxsum1b[ii])/(shellfluxmax1[ii] + 1.0e-50);
        //relflux2 = fabs(shellfluxsum2a[ii]+shellfluxsum2b[ii])/(shellfluxmax2[ii] + 1.0e-50);
        shellrelfluxmax0 = MAX(shellrelfluxmax0,relflux0);
        shellrelfluxmax1 = MAX(shellrelfluxmax1,relflux1);
        shellrelfluxmax2 = MAX(shellrelfluxmax2,relflux2);
        shellrelfluxmax = MAX(shellrelfluxmax0,MAX(shellrelfluxmax1,shellrelfluxmax2));
        // BEWARE:  Since fluxes in the d-direction are proportional to velocities in the d-direction,
        //   which may be close to zero, these quantities can be relatively large even when absolutely small
        printf("[check_flux_differencing]:  ii = %d, shellrelfluxsum0 = %e, shellrelfluxsum1 = %e, shellrelfluxsum2 = %e\n",
          ii,relflux0,relflux1,relflux2);
      }
    }

    if (verbose) {
      fprintf(stderr,"[check_flux_differencing]:  max relative flux difference = %e\n",relfluxmax);
      if (check_shells) fprintf(stderr,"[check_flux_differencing]:  shell max relative flux difference = %e\n",shellrelfluxmax);
    }
    fprintf(fp,"%e %1.16e\n",t,relfluxmax);
    fflush(fp);
    if (check_shells) fprintf(fps,"%e %1.16e %1.16e %1.16e %1.16e\n",t,shellrelfluxmax0,shellrelfluxmax1,shellrelfluxmax2,shellrelfluxmax);
    fflush(fps);
  }

  free(fluxsum0a);
  free(fluxsum0b);
  free(fluxsum1a);
  free(fluxsum1b);
  free(fluxsum2a);
  free(fluxsum2b);
  free(fluxmax0);
  free(fluxmax1);
  free(fluxmax2);
  if (check_shells) {
    free(shellfluxsum0a);
    free(shellfluxsum0b);
    free(shellfluxsum1a);
    free(shellfluxsum1b);
    free(shellfluxsum2a);
    free(shellfluxsum2b);
    free(shellfluxmax0);
    free(shellfluxmax1);
    free(shellfluxmax2);
  }

  return relfluxmax;
}

double calculate_total_mass(int verbose)
{
  int ii,jj,kk;
  double total_mass = 0.0;
  static FILE* fp = NULL;
  static int firstc = 1;

  if (firstc) {
    if (mpi_io_proc()) fp = fopen("total_mass.txt","w");
    firstc = 0;
  }

  ZLOOP { total_mass += NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO)*ND_ELEM_LINEAR(geom,ii,jj,kk).volume; }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE,&total_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #endif

  if (mpi_io_proc()) {
    if (verbose) fprintf(stderr,"[calculate_total_mass]:  Total mass = %1.15e\n",total_mass);
    fprintf(fp,"%e %1.16e\n",t,total_mass);
    fflush(fp);
  }

  return total_mass;
}

#if (GEOM==SPHERICAL)
static double *ms,*mr,*grav,*dgrav;
#if (GR_MONOPOLE==TRUE)
static double *rhoavg,*pgasavg,*pradavg,*uavg,*vavg,*eradavg,*vdF,*gtov;
#endif
double calculate_grav_energy()
{
  static int firstc = 1;
  int ii,jj,kk,dd,g;
  int iter,j,k,jp,kp,njp,nkp;
  double momsrc,Etot,vdotF,Fmag,fred,chi,Prr,vol;
  double Ftot[SPACEDIM];
  double gtov_avg,mlast,mstar;
  double Vshell,xl,A,B,vedge,vr,vl,area;
  double pgasedge,pradedge,rhoedge,uedge;
  double C,xr,grav_energy,gravcon;
  if (firstc) {
    ms            = malloc_rank1(n1,   sizeof *ms           );
    mr            = malloc_rank1(n1+1, sizeof *mr           );
    grav          = malloc_rank1(n1+1, sizeof *grav         );
    dgrav         = malloc_rank1(n1,   sizeof *dgrav        );
    #if (GR_MONOPOLE==TRUE)
    rhoavg        = malloc_rank1(n1,   sizeof *rhoavg       );
    pgasavg       = malloc_rank1(n1,   sizeof *pgasavg      );
    pradavg       = malloc_rank1(n1,   sizeof *pradavg      );
    uavg          = malloc_rank1(n1,   sizeof *uavg         );
    vavg          = malloc_rank1(n1,   sizeof *vavg         );
    eradavg       = malloc_rank1(n1,   sizeof *eradavg      );
    vdF           = malloc_rank1(n1,   sizeof *vdF          );
    gtov          = malloc_rank1(n1+1, sizeof *gtov         );
    gr_grav       = malloc_rank1(n1,   sizeof *gr_grav      );
    #endif
    firstc=0;
  }

  // assumes that the first dimension is radial
  for (ii=0; ii<n1; ii++) {
    ms[ii]      = 0.0;
    #if (GR_MONOPOLE==TRUE)
    rhoavg[ii]  = 0.0;
    pgasavg[ii] = 0.0;
    pradavg[ii] = 0.0;
    uavg[ii]    = 0.0;
    vavg[ii]    = 0.0;
    eradavg[ii] = 0.0;
    vdF[ii]     = 0.0;
    #endif
  }

  #if (GR_MONOPOLE==TRUE)
  ZLOOP {
    #if (DO_RADIATION==TRUE)
    Etot = 0.0;
    SLOOP { Ftot[dd] = 0.0; }
    GLOOP {
      Etot += NDP_ELEM_LINEAR(sim_p,ii,jj,kk,irad1+g);
      SLOOP { Ftot[dd] += NDP_ELEM_LINEAR(sim_p,ii,jj,kk,ifrad1+g*NDIM+dd); }
    }
    Etot /= SQR(CLIGHT);
    vdotF = geom_dot(Ftot, &NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1), ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0])/pow(CLIGHT,4);
    // calculate Prr for the M1 closure
    Fmag  = sqrt(geom_dot(Ftot,Ftot,ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0]));
    fred  = Fmag/(CLIGHT*fabs(Etot) + 1.0e-16);
    fred  = MAX(0.0,MIN(1.0,fred));
    chi   = (3.0 + 4.0*SQR(fred))/(5.0 + 2.0*sqrt(4.0 - 3.0*SQR(fred)));
    Prr   = 0.5*Etot*((1.0-chi) + (3.0*chi-1.0)*SQR(Ftot[0])*ND_ELEM_LINEAR(geom,ii,jj,kk).gcov[0][0]/(SQR(Fmag) + 1.0e-16));
    #else
    Etot  = 0.0;
    vdotF = 0.0;
    Prr   = 0.0;
    #endif
    vol          = ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
    rhoavg[ii]  += vol*NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
    pgasavg[ii] += vol*NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,PRESS);
    pradavg[ii] += vol*Prr;
    uavg[ii]    += vol*NDP_ELEM_LINEAR(sim_p,ii,jj,kk,UU);
    vavg[ii]    += vol*NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1);
    eradavg[ii] += vol*Etot;
    vdF[ii]     += vol*vdotF;
  }
  #else  /* NOT GR_MONOPOLE */
  ZLOOP { ms[ii] += ND_ELEM_LINEAR(geom,ii,jj,kk).volume*NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO); }
  #endif  /* ifdef GR_MONOPOLE */

#if (USE_MPI==TRUE)
  #if (GR_MONOPOLE==TRUE)
  MPI_Allreduce(MPI_IN_PLACE,rhoavg,  n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,pgasavg, n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,pradavg, n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,uavg,    n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,vavg,    n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,eradavg, n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,vdF,     n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #else
  MPI_Allreduce(MPI_IN_PLACE,ms,      n1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #endif
#endif


  #if (GR_MONOPOLE==TRUE)
  // need to get volume averages
  for (ii=0; ii<n1; ii++) {
    Vshell       = 1.0/(4.0*M_PI/3.0 * (pow(redge[ii+1],3) - pow(redge[ii],3)));
    rhoavg[ii]  *= Vshell;
    pgasavg[ii] *= Vshell;
    pradavg[ii] *= Vshell;
    uavg[ii]    *= Vshell;
    vavg[ii]    *= Vshell;
    eradavg[ii] *= Vshell;
    vdF[ii]     *= Vshell;
  }

  // now iterate to converge on Mtov and Gamma
  mr[0]   = 0.0;
  grav[0] = 0.0;
  if (firstc) {
    // get a reasonable estimate for gtov
    gtov[0] = 1.0;
    mstar   = 0.0;
    for (ii=1; ii<=n1; ii++) {
      Vshell   = 4.0*M_PI/3.0 * (pow(redge[ii],3) - pow(redge[ii-1],3));
      mstar   += Vshell*rhoavg[ii-1];
      gtov[ii] = sqrt(1.0 - 2.0*GNEWT*mstar/(SQR(CLIGHT)*redge[ii]));
    }
    firstc = 0;
  }

  iter = 0;
  while (1) {
    mstar = 0.0;
    for (ii=1; ii<=n1; ii++) {
      xl = startx[0] + ii*dx[0];
      A  = (gtov[ii] - gtov[ii-1])/dx[0];
      B  = gtov[ii-1] - A*xl;
      gtov_avg = A*beta0[ii]/Gamma0[ii] + B;
      Vshell = 4.0*M_PI/3.0 * (pow(redge[ii],3) - pow(redge[ii-1],3));
      vedge = (ii < n1) ? 0.5*(vavg[ii-1]+vavg[ii]) : vavg[ii-1];
      gtov[ii] = sqrt(1.0 + pow(vedge/CLIGHT,2) - 2.0*GNEWT*(mr[ii-1] + Vshell*((rhoavg[ii-1]+uavg[ii-1]/SQR(CLIGHT) + eradavg[ii-1])*gtov_avg + vdF[ii-1]))/(redge[ii]*SQR(CLIGHT)));
      A = (gtov[ii] - gtov[ii-1])/dx[0];
      B = gtov[ii-1] - A*xl;
      gtov_avg = A*beta0[ii]/Gamma0[ii] + B;
      mr[ii] = mr[ii-1] + Vshell*((rhoavg[ii-1]+uavg[ii-1]/SQR(CLIGHT) + eradavg[ii-1])*gtov_avg + vdF[ii-1]);
      mstar += Vshell*rhoavg[ii-1];
    }
    if (iter > 0) {
      if (fabs(mlast - mr[n1])/mr[n1] < 1.0e-10) {
        break;
      }
    }
    mlast = mr[n1];
    iter++;
    if (iter > 20) {
      if (mpi_io_proc()) fprintf(stderr,"iter = %d\n", iter);
      fprintf(stderr,"%d %g\n", myrank, mlast);
      fflush(stderr);
      #if (USE_MPI==TRUE)
      MPI_Barrier(MPI_COMM_WORLD);
      #endif
      exit(3);
    }
  }
  if (mpi_io_proc() && iter > 3) fprintf(stderr,"gr monopole required %d iterations\n", iter);

  for (ii=1; ii<n1; ii++) {
    pgasedge = 0.5*(pgasavg[ii-1] + pgasavg[ii]);
    pradedge = 0.5*(pradavg[ii-1] + pradavg[ii]);
    rhoedge  = 0.5*( rhoavg[ii-1] +  rhoavg[ii]);
    uedge    = 0.5*(   uavg[ii-1] +    uavg[ii]);
    grav[ii] = -GNEWT*(mr[ii] + 4.0*M_PI*pow(redge[ii],3)*(pgasedge+pradedge)/SQR(CLIGHT))/pow(redge[ii]*gtov[ii],2) * (rhoedge + (uedge+pgasedge)/SQR(CLIGHT))/rhoedge * dr_dx(startx[0]+ii*dx[0]);
  }
  grav[n1] = -GNEWT*(mr[n1] + 4.0*M_PI*pow(redge[n1],3)*(pgasavg[n1-1]+pradavg[n1-1])/SQR(CLIGHT))/pow(redge[n1]*gtov[n1],2) * (rhoavg[n1-1] + (uavg[n1-1]+pgasavg[n1-1])/SQR(CLIGHT))/rhoavg[n1-1] * dr_dx(startx[0]+n1*dx[0]);

  // now calculate the "GR" potential
  gr_lapse_edge[n1] = -GNEWT*mstar/redge[n1];
  for (ii=n1-1; ii>=0; ii--) {
    gr_lapse_edge[ii] = gr_lapse_edge[ii+1] + 0.5*(grav[ii]+grav[ii+1])*dx[0];   // + because g=-dphi/dr & I'm integrating inwards
  }

  for (ii=0; ii<n1; ii++) {
    xl = startx[0] + ii*dx[0];
    A  = (grav[ii+1]-grav[ii])/dx[0];
    B  = grav[ii] - A*xl;
    dgrav[ii]   = A;
    grav[ii]    = A*beta0[ii]/Gamma0[ii] + B;
  }

  #else  /* NOT GR_MONOPOLE */

  mr[0]   = 0.0;
  grav[0] = 0.0;
  for (ii=1; ii<=n1; ii++) {
    mr[ii]   = mr[ii-1] + ms[ii-1];
    grav[ii] = -GNEWT*mr[ii]/SQR(redge[ii])*dr_dx(startx[0]+ii*dx[0]);
  }

  for (ii=0; ii<n1; ii++) {
    xl = startx[0] + ii*dx[0];
    A  = (grav[ii+1]-grav[ii])/dx[0];
    B  = grav[ii] - A*xl;
    dgrav[ii] = A;
    grav[ii]  = A*beta0[ii]/Gamma0[ii] + B;
  }
  #endif  /* ifdef GR_MONOPOLE */

  // Compute the binding energy as in eq. (2.17) of Binney & Tremaine (2nd ed.)
  grav_energy = 0.0;
  ZLOOP {
    gravcon = grav[ii]*ND_ELEM_LINEAR(geom,ii,jj,kk).gcon[0];  // raise indices
    grav_energy += grav[ii]*gravcon*ND_ELEM_LINEAR(geom,ii,jj,kk).volume;  // |g|^2 = g_1 g^1
  }
  grav_energy *= -1.0/(8.0*M_PI*GNEWT);
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE,&grav_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #endif

  return grav_energy;
}

double calculate_total_energy(int verbose)
{
  static double* u = NULL;
  static int firstc = 1;
  int ii,jj,kk;
  double Etot = 0.0;
  double grav_energy,total_energy;
  static FILE* fp = NULL;

  if (firstc) {
    u = malloc_rank1(nvars, sizeof *u);
    if (mpi_io_proc()) fp = fopen("total_energy.txt","w");
    firstc = 0;
  }

  ZLOOP {
    p_to_u(ND_ELEM_LINEAR(sim_p,ii,jj,kk),u,&ND_ELEM_LINEAR(geom,ii,jj,kk),nhydro);
    Etot += u[ETOT]*ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
  }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE,&Etot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  #endif
  grav_energy = calculate_grav_energy();
  total_energy = Etot + grav_energy;

  if (mpi_io_proc()) {
    if (verbose) fprintf(stderr,"[calculate_total_energy]:  Total energy = %1.15e\n",total_energy);
    fprintf(fp,"%e %1.16e %1.16e %1.16e\n",t,Etot,grav_energy,total_energy);
    fflush(fp);
  }

  return total_energy;
}
#endif /* GEOM==SPHERICAL */

int new_timer(char const * name) {
  int timer_index = -1;
  if(ntimers == MAX_NTIMERS) {
    fprintf(stdout, "Could not allocate timer \"%s\": too many timers!\n", name);
    fflush(stdout);
    abort();
  }
  else {
    timer_index = ntimers++;
    int ierr = timer_init(name, &timers[timer_index]);
    if(ierr != ERROR_NONE) {
      fprintf(stdout, "Could not allocate timer \"%s\": internal error.\n", name);
      fflush(stdout);
      abort();
    }
  }
  return timer_index;
}

Timer * get_timer(int const timer_index) {
  if(timer_index < 0 || timer_index >= ntimers) {
    fprintf(stdout, "Requesting invalid timer: \"%d\"\n", timer_index);
    fflush(stdout);
    abort();
  }
  return timers[timer_index];
}

void free_timers() {
  for(int tidx = 0; tidx < ntimers; ++tidx) {
    timer_free(timers[tidx]);
    timers[tidx] = NULL;
  }
  ntimers = 0;
}

double dynamic_root_find(double (*f)(double), double (*dfdx)(double),
                         double x0, double xlow, double xup, const double tol,
                         int max_iter) {
  double next, x1, y0, y1, ylow, yup;
  ylow = (*f)(xlow);
  yup = (*f)(xup);
  if (ylow*yup > 0.0) {
    fprintf(stderr,"dynamic_root_find [root not bracketed]:  ylow=%g, yup=%g\n",ylow,yup);
    exit(1);
  }
  if (fabs(ylow) < tol) return xlow;
  if (fabs(yup) < tol) return xup;
  y0 = (*f)(x0);
  x1 = 1.01 * x0;
  y1 = (*f)(x1);
  while (fabsl(y0) > tol && max_iter-- > 0) {
    // If x0 is outside the bracket use bisection
    if (x0 < xlow || x0 > xup) {
      x1 = x0;
      y1 = y0;
      if ((xlow > 0) && (xup > 0)) {
        x0 = sqrt(xlow*xup);  // Geometric bisection
      } else {
        x0 = 0.5*(xlow+xup);  // Arithmetic bisection
      }
    } else {
      // update the brackets
      if (ylow*y0 < 0.0) {
        xup = x0;
        yup = y0;
      } else {
        xlow = x0;
        ylow = y0;
      }
      if ((fabs(ylow) > 1e-2) || !(dfdx)) {
        // Secant method step
        next = xlow - y0*(x0-x1)/(y0-y1);
      } else {
        // Newton step
        next = x0 - y0 / (*dfdx)(x0);
      }
      x1 = x0;
      y1 = y0;
      x0 = next;
    }
    y0 = (*f)(x0);
  }
  if (max_iter <= 0) {
    fprintf(stderr,"dynamic_root_find [max_iter reached]:  y0=%g\n",y0);
    exit(1);
  }
  return x0;
}

