
#include "../decs.h"

void init_coords()
{
  return;
}

void ijk_to_x(int i, int j, int k, double x[NDIM])
{
  int dd;
  int index[] = {i, j, k};

  DLOOP {
    x[dd] = (index[dd]+0.5)*dx[dd] + startx[dd];
  }

  return;
}

void ijk_to_Cart(int i, int j, int k, double *x)
{
  x[0] = startx[0] + i*dx[0];
  #if (NDIM>1)
  x[1] = startx[1] + j*dx[1];
  #if (NDIM==3)
  x[2] = startx[2] + k*dx[2];
  #endif
  #endif
  return;
}

double ijk_to_r(int i, int j, int k, double *rhat)
{
  int dd;
  double x[NDIM],r;
  
  ijk_to_x(i,j,k,x);
  r = 0.0;
  DLOOP { r += SQR(x[dd]); }
  r = sqrt(r);
  DLOOP { rhat[dd] = x[dd]/r; }

  return r;
}

void init_interp_vol()
{
	int i,j,k,n2s,n3s;
	double x1,x0;

	for (i=istart[0]-NG; i<istop[0]+NG; i++) {
		x0 = startx[0] + i*dx[0];
		x1 = x0 + dx[0];
		alpha0[i] = (CUBE(x1) - CUBE(x0))/3.0;
		beta0[i]  = (SQR(x1) - SQR(x0))/2.0;
		Gamma0[i] = dx[0];
	}

	#if (NDIM>1)
	Gamma1s = malloc_rank1(1, sizeof *Gamma1s);
	beta1s  = malloc_rank1(1, sizeof *beta1s );
	alpha1s = malloc_rank1(1, sizeof *alpha1s);

	n2s = global_grid_dims[1] + 2*NG;
	Gamma1s[0]  = malloc_rank1(n2s, sizeof *Gamma1s[0]);
	Gamma1s[0] += NG - istart[1];
	beta1s [0]  = malloc_rank1(n2s, sizeof *beta1s [0]);
	beta1s [0] += NG - istart[1];
	alpha1s[0]  = malloc_rank1(n2s, sizeof *alpha1s[0]);
	alpha1s[0] += NG - istart[1];

	for (j=istart[1]-NG; j<istop[1]+NG; j++) {
		x0 = startx[1] + j*dx[1];
		x1 = x0 + dx[1];
		alpha1s[0][j] = (CUBE(x1) - CUBE(x0))/3.0;
		beta1s[0][j]  = (SQR(x1) - SQR(x0))/2.0;
		Gamma1s[0][j] = dx[1];
	}
	#endif

	#if (NDIM==3)
	Gamma2s = malloc_rank1(1, sizeof *Gamma2s);
	beta2s  = malloc_rank1(1, sizeof *beta2s );
	alpha2s = malloc_rank1(1, sizeof *alpha2s);

	n3s = global_grid_dims[2] + 2*NG;
	Gamma2s[0]  = malloc_rank1(n3s, sizeof *Gamma2s[0]);
	Gamma2s[0] += NG - istart[2];
	beta2s [0]  = malloc_rank1(n3s, sizeof *beta2s [0]);
	beta2s [0] += NG - istart[2];
	alpha2s[0]  = malloc_rank1(n3s, sizeof *alpha2s[0]);
	alpha2s[0] += NG - istart[2];

	for (k=istart[2]-NG; k<istop[2]+NG; k++) {
		x0 = startx[2] + k*dx[2];
		x1 = x0 + dx[2];
		alpha2s[0][k] = (CUBE(x1) - CUBE(x0))/3.0;
		beta2s[0][k]  = (SQR(x1) - SQR(x0))/2.0;
		Gamma2s[0][k] = dx[2];
	}
	#endif

	return;
}

void calc_interp_coeffs1(int i, double *alpha, double *beta, double *Gamma, double *x)
{
	for (int j=istart[1]-NG; j<istop[1]+NG; j++) {
		alpha[j-istart[1]] = alpha1s[0][j];
		beta[j-istart[1]] = beta1s[0][j];
		Gamma[j-istart[1]] = Gamma1s[0][j];
		x[j-istart[1]] = startx[1] + j*dx[1];
	}
	return;
}

void calc_interp_coeffs2(int i, int j, double *alpha, double *beta, double *Gamma, double *x)
{
	for (int k=istart[2]-NG; k<istop[2]+NG; k++) {
		alpha[k-istart[2]] = alpha2s[0][k];
		beta[k-istart[2]] = beta2s[0][k];
		Gamma[k-istart[2]] = Gamma2s[0][k];
		x[k-istart[2]] = startx[2] + k*dx[2];
	}
	return;
}

void vec_transform_to_xyz(double *vp, int i, int j, int k)
{
	return;
}

void init_volume()
{
	int ii,jj,kk,dd;

	for (ii=istart[0]-NG; ii<istop[0]+NG; ii++) {
		#if (NDIM>1)
		for (jj=istart[1]-NG; jj<istop[1]+NG; jj++) {
		#endif
  		#if (NDIM==3)
  		for (kk=istart[2]-NG; kk<istop[2]+NG; kk++) {
  		#endif
    		ND_ELEM(geom,ii,jj,kk).volume = 1.0;
    		DLOOP {
    			ND_ELEM(geom,ii,jj,kk).volume *= dx[dd];
    		}
  		#if (NDIM==3)
  		}
		#endif
		#if (NDIM>1)
		}
		#endif
	}
}

	

void init_area()
{
	int ii,jj,kk,loc,dd;

	for (ii=istart[0]; ii<=istop[0]; ii++) {
		#if (NDIM>1)
		for (jj=istart[1]; jj<=istop[1]; jj++) {
		#endif
  		#if (NDIM==3)
  		for (kk=istart[2]; kk<=istop[2]; kk++) {
  		#endif
    		for (loc=0; loc<NDIM; loc++) {		
    			ND_ELEM(geom,ii,jj,kk).area[loc] = 1.0;
    			DLOOP {
    				ND_ELEM(geom,ii,jj,kk).area[loc] *= dx[dd];
    			}
    			ND_ELEM(geom,ii,jj,kk).area[loc] /= dx[loc];
    		}
  		#if (NDIM==3)
  		}
  		#endif
		#if (NDIM>1)
		}
		#endif
	}

	return;
}

void init_conn()
{
	int ii,jj,kk,i,j,k;

	for (ii=istart[0]; ii<istop[0]; ii++) {
		#if (NDIM>1)
		for (jj=istart[1]; jj<istop[1]; jj++) {
		#endif
  		#if (NDIM==3)
  		for (kk=istart[2]; kk<istop[2]; kk++) {
  		#endif
    		for (i=0; i<SPACEDIM; i++) {
    			for (j=0; j<SPACEDIM; j++) {
    				for (k=0; k<SPACEDIM; k++) {
    					ND_ELEM(geom,ii,jj,kk).conn[i][j][k] = 0.0;
    				}
    			}
    		}
  		#if (NDIM==3)
  		}
  		#endif
		#if (NDIM>1)
		}
		#endif
	}

	return;
}

void init_gcov()
{
	int ii,jj,kk,loc,i,j;

	for (ii=istart[0]-NG; ii<istop[0]+NG; ii++) {
		#if (NDIM>1)
		for (jj=istart[1]-NG; jj<istop[1]+NG; jj++) {
		#endif
  		#if (NDIM==3)
  		for (kk=istart[2]-NG; kk<istop[2]+NG; kk++) {
  		#endif
    		for (loc=0; loc<=NDIM; loc++) {
    			for (i=0; i<SPACEDIM; i++) {
    				ND_ELEM(geom,ii,jj,kk).gcov[loc][i] = 1.0;
    			}
    		}
  		#if (NDIM==3)
  		}
  		#endif
		#if (NDIM>1)
		}
		#endif
	}

	return;
}

void init_gcon()
{
	int ii,jj,kk,loc,i,j;

	for (ii=istart[0]-NG; ii<istop[0]+NG; ii++) {
		#if (NDIM>1)
		for (jj=istart[1]-NG; jj<istop[1]+NG; jj++) {
		#endif
  		#if (NDIM==3)
  		for (kk=istart[2]-NG; kk<istop[2]+NG; kk++) {
  		#endif
    		for (i=0; i<SPACEDIM; i++) {
    			ND_ELEM(geom,ii,jj,kk).gcon[i] = 1.0;
    		}
  		#if (NDIM==3)
  		}
  		#endif
		#if (NDIM>1)
		}
		#endif
	}

	return;
}

void init_scale()
{
	int ii,jj,kk,dd,loc;

	for (ii=istart[0]-NG; ii<istop[0]+NG; ii++) {
		#if (NDIM>1)
		for (jj=istart[1]-NG; jj<istop[1]+NG; jj++) {
		#endif
  		#if (NDIM==3)
  		for (kk=istart[2]-NG; kk<istop[2]+NG; kk++) {
  		#endif
    		for (loc=0; loc<=NDIM; loc++) {
    			SLOOP {
    				ND_ELEM(geom,ii,jj,kk).scale[loc][dd] = 1.0;
    			}
    		}
  		#if (NDIM==3)
  		}
  		#endif
		#if (NDIM>1)
		}
		#endif
	}

	return;
}

// TODO: init_dendritic should not need to call these for a cartesian grid
double r_of_x(double x) {
    return x;
}
double th_of_x(double x) {
    return x;
}
