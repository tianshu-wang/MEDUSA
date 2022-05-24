
#include "../decs.h"

void ijk_to_x(int i, int j, int k, double x[SPACEDIM]){

	x[0] = exp((i+0.5)*dx[0] + startx[0]);
	#if(NDIM>1)
	x[1] = (j+0.5)*dx[1] + startx[1];
	#else
	x[1] = 0.5;
	#endif
	#if(NDIM==3)
	x[2] = (k+0.5)*dx[2] + startx[2];
	#else
	x[2] = 0.;
	#endif

	return;
}

double ijk_to_r(int i, int j, int k, double rhat[SPACEDIM]){

	double r;

	r = exp((i+0.5)*dx[0] + startx[0]);

	rhat[0] = 1.;
	rhat[1] = 0.;
	rhat[2] = 0.;

	return r;
}

void ijk_to_r_dr(int i, int j, int k, double r[]) {

	r[0] = exp(i*dx[0] + startx[0]);
	r[1] = exp((i+0.5)*dx[0] + startx[0]);
	r[2] = exp((i+1)*dx[0] + startx[0]);
}

void vec_transform_to_xyz(float *vp, int i, int j, int k) {

	double x[SPACEDIM];
	float v[SPACEDIM];
	double lam[SPACEDIM][SPACEDIM];
	double r,sth,cth,sph,cph;
	ijk_to_x(i, j, k, x);

	r = x[0];
	sth = sin(x[1]);
	cth = cos(x[1]);
	sph = sin(x[2]);
	cph = cos(x[2]);

	// x
	lam[0][0] = sth*cph;	// r
	lam[0][1] = r*cth*cph;	// theta
	lam[0][2] = -r*sth*sph;	// phi

	// y
	lam[1][0] = sth*sph;
	lam[1][1] = r*cth*sph;
	lam[1][2] = r*sth*cph;


	// z
	lam[2][0] = cth;
	lam[2][1] = -r*sth;
	lam[2][2] = 0.;

	v[0] = v[1] = v[2] = 0.;

	for(i=0;i<SPACEDIM;i++) {
		for(j=0;j<SPACEDIM;j++) {
			v[i] += lam[i][j]*vp[j];
		}
	}

	#if(NDIM==2)
	vp[0] = v[0];
	vp[1] = v[2];
	vp[2] = v[1];
	#else
	vp[0] = v[0];
	vp[1] = v[1];
	vp[2] = v[2];
	#endif

	return;
}

void init_coords() {

	for(ii=-NG;ii<my_grid_dims[0]+NG;ii++) {
		redge[ii] = exp(ii*dx[0]+startx[0]);
		rcenter[ii] = exp((ii+0.5)*dx[0]+startx[0]);
		dr[ii] = exp((ii+1)*dx[0]+startx[0]) - redge[ii];
		fprintf(stderr,"%d %g\n", ii, redge[ii]);
	}

}

void init_gdet() {

	int loc;
	double x[SPACEDIM];
	double r[3];

	for(ii=-1;ii<my_grid_dims[0];ii++) {
		#if(NDIM==1)
		ijk_to_x(ii,0,0,x);
		ijk_to_r_dr(ii,0,0,r);
		ND_ELEM(geom,ii,jj,kk).gdet[0] = x[0]*x[0]*x[0];
		x[0] = r[2];
		ND_ELEM(geom,ii,jj,kk).gdet[1] = x[0]*x[0]*x[0];
		#endif
		#if(NDIM>1)
		for(jj=-1;jj<my_grid_dims[1]/2;jj++) {
			ijk_to_x(ii,jj,0,x);
			ijk_to_r_dr(ii,jj,0,r);
			#if(NDIM==3)
			for(kk=-1;kk<my_grid_dims[2];kk++) {
			#endif
			ND_ELEM(geom,ii,jj,kk).gdet[0] = x[0]*x[0]*x[0]*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			x[0] = r[2];
			ND_ELEM(geom,ii,jj,kk).gdet[1] = x[0]*x[0]*x[0]*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			x[0] = r[1];
			x[1] += 0.5*dx[1];
			ND_ELEM(geom,ii,jj,kk).gdet[2] = x[0]*x[0]*x[0]*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			#if(NDIM==3)
			ND_ELEM(geom,ii,jj,kk).gdet[3] = ND_ELEM(geom,ii,jj,kk).gdet[0];
			}
			#endif
		}
		for(jj=my_grid_dims[1]/2;jj<my_grid_dims[1];jj++) {
			#if(NDIM==3)
			for(kk=-1;kk<my_grid_dims[2];kk++) {
			#endif
			ND_ELEM(geom,ii,jj,kk).gdet[0] = ND_ELEM(geom,ii,my_grid_dims[1]-jj-1,kk).gdet[0];
			ND_ELEM(geom,ii,jj,kk).gdet[1] = ND_ELEM(geom,ii,my_grid_dims[1]-jj-1,kk).gdet[1];
			ND_ELEM(geom,ii,jj,kk).gdet[2] = ND_ELEM(geom,ii,my_grid_dims[1]-jj-2,kk).gdet[2];
			#if(NDIM==3)
			ND_ELEM(geom,ii,jj,kk).gdet[3] = ND_ELEM(geom,ii,jj,kk).gdet[0];
			}
			#endif
		}
		#endif
		fprintf(stderr,"gdet: %d %g %g\n", ii, ND_ELEM(geom,ii,jj,kk).gdet[0], ND_ELEM(geom,ii,jj,kk).gdet[1]);
	}

	//ND_ELEM(geom,-1,jj,kk).gdet[0] = ND_ELEM(geom,0,jj,kk).gdet[0];
	//ND_ELEM(geom,-1,jj,kk).gdet[1] = 0;

	return;
}

void init_conn() {

	double x[SPACEDIM],sth,cth,discrete_cth;

	for(ii=-1;ii<my_grid_dims[0];ii++) {
		#if(NDIM>1)
		for(jj=-1;jj<my_grid_dims[1];jj++) {
		#endif
		#if(NDIM==3)
		for(kk=-1;kk<my_grid_dims[2];kk++) {
		#endif
		ijk_to_x(ii,jj,kk,x);
		memset(&ND_ELEM(geom,ii,jj,kk).conn[0][0][0], 0, 27*sizeof(double));
		sth = sin(x[1]);//0.5*(sin(x[1]) + sin(M_PI-x[1]));
		cth = cos(x[1]);//0.5*(cos(x[1]) - cos(x[1]+M_PI);
		discrete_cth = (sin(x[1]+0.5*dx[1])-sin(x[1]-0.5*dx[1]))/dx[1];
		ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = 2*sinh(1.5*dx[0])/(3.*dx[0]);
		ND_ELEM(geom,ii,jj,kk).conn[0][1][1] = -1.;
		ND_ELEM(geom,ii,jj,kk).conn[0][2][2] = -sth*sth;
		ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = 2.*sinh(1.5*dx[0])/(3.*dx[0]);
		ND_ELEM(geom,ii,jj,kk).conn[1][1][0] = 1.;
		ND_ELEM(geom,ii,jj,kk).conn[1][2][2] = -sth*cth;
		ND_ELEM(geom,ii,jj,kk).conn[2][0][2] = 2.*sinh(1.5*dx[0])/(3.*dx[0]);
		ND_ELEM(geom,ii,jj,kk).conn[2][1][2] = discrete_cth/sth;
		ND_ELEM(geom,ii,jj,kk).conn[2][2][0] = 1.;
		ND_ELEM(geom,ii,jj,kk).conn[2][2][1] = discrete_cth/sth;
		//fprintf(stderr,"conn: %g %g %g %g %g %g %g %g %g %g\n", ND_ELEM(geom,ii,jj,kk).conn[0][0][0], ND_ELEM(geom,ii,jj,kk).conn[0][1][1], ND_ELEM(geom,ii,jj,kk).conn[0][2][2], ND_ELEM(geom,ii,jj,kk).conn[1][0][1], ND_ELEM(geom,ii,jj,kk).conn[1][1][0], ND_ELEM(geom,ii,jj,kk).conn[1][2][2], ND_ELEM(geom,ii,jj,kk).conn[2][0][2], ND_ELEM(geom,ii,jj,kk).conn[2][1][2], ND_ELEM(geom,ii,jj,kk).conn[2][2][0], ND_ELEM(geom,ii,jj,kk).conn[2][2][1] = discrete_cth/sth);

		#if(NDIM==3)
		}
		#endif
		#if(NDIM>1)
		}
		#endif
	}

	fprintf(stderr,"conn2 -> %g\n", ND_ELEM(geom,0,0,0).conn[1][0][1]);

	return;
}

void init_gcov() {

	double x[SPACEDIM],r[3];
	int loc,i,j;

	for(ii=-1;ii<my_grid_dims[0];ii++) {
	#if(NDIM>1)
	for(jj=-1;jj<my_grid_dims[1]/2;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=-1;kk<my_grid_dims[2];kk++) {
	#endif
		for(loc=0;loc<=NDIM;loc++) {
			for(i=0;i<SPACEDIM;i++) {
				for(j=0;j<SPACEDIM;j++) {
					ND_ELEM(geom,ii,jj,kk).gcov[loc][i][j] = 0.;
				}
			}
		}
		ijk_to_x(ii,jj,kk,x);
		ijk_to_r_dr(ii,jj,kk,r);
		// centered
		ND_ELEM(geom,ii,jj,kk).gcov[0][0][0] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[0][1][1] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[0][2][2] = pow(x[0]*sin(x[1]),2.);
		// + 0.5*dx[0]
		x[0] = r[2];
		ND_ELEM(geom,ii,jj,kk).gcov[1][0][0] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[1][1][1] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[1][2][2] = pow(x[0]*sin(x[1]),2.);
		// + 0.5*dx[1]
		x[0] = r[1];
		x[1] += 0.5*dx[1];
		ND_ELEM(geom,ii,jj,kk).gcov[2][0][0] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[2][1][1] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[2][2][2] = pow(x[0]*sin(x[1]),2.);
		
	#if(NDIM==3)
		// + 0.5*dx[2]
		ND_ELEM(geom,ii,jj,kk).gcov[3][0][0] = ND_ELEM(geom,ii,jj,kk).gcov[0][0][0];
		ND_ELEM(geom,ii,jj,kk).gcov[3][1][1] = ND_ELEM(geom,ii,jj,kk).gcov[0][1][1];
		ND_ELEM(geom,ii,jj,kk).gcov[3][2][2] = ND_ELEM(geom,ii,jj,kk).gcov[0][2][2];

	}
	#endif
	#if(NDIM>1)
	}
	for(jj=my_grid_dims[1]/2;jj<my_grid_dims[1];jj++) {
	#endif
	#if(NDIM==3)
	for(kk=-1;kk<my_grid_dims[2];kk++) {
	#endif
		for(i=0;i<SPACEDIM;i++) {
			for(j=0;j<SPACEDIM;j++) {
				ND_ELEM(geom,ii,jj,kk).gcov[0][i][j] = ND_ELEM(geom,ii,my_grid_dims[1]-1-jj,kk).gcov[0][i][j];
				ND_ELEM(geom,ii,jj,kk).gcov[1][i][j] = ND_ELEM(geom,ii,my_grid_dims[1]-1-jj,kk).gcov[1][i][j];
				ND_ELEM(geom,ii,jj,kk).gcov[2][i][j] = ND_ELEM(geom,ii,my_grid_dims[1]-2-jj,kk).gcov[2][i][j];
				#if(NDIM==3)
				ND_ELEM(geom,ii,jj,kk).gcov[3][i][j] = ND_ELEM(geom,ii,jj,kk).gcov[0][i][j];
				#endif
			}
		}
	#if(NDIM==3)
	}
	#endif
	#if(NDIM>1)
	}
	#endif
	fprintf(stderr,"gcov: %d %g\n", ii, ND_ELEM(geom,ii,jj,kk).gcov[0][0][0]);
	}

	return;
}

void init_gcon() {

	int loc,i,j;

	for(ii=-1;ii<my_grid_dims[0];ii++) {
	#if(NDIM>1)
	for(jj=-1;jj<my_grid_dims[1];jj++) {
	#endif
	#if(NDIM==3)
	for(kk=-1;kk<my_grid_dims[2];kk++) {
	#endif
		for(i=0;i<SPACEDIM;i++) {
			for(j=0;j<SPACEDIM;j++) {
				if(fabs(ND_ELEM(geom,ii,jj,kk).gcov[0][i][j]) > 1.e-50) {
					ND_ELEM(geom,ii,jj,kk).gcon[i][j] = 1./ND_ELEM(geom,ii,jj,kk).gcov[0][i][j];
				} else {
					ND_ELEM(geom,ii,jj,kk).gcon[i][j] = 0.;
				}
			}
		}
	#if(NDIM==3)
	}
	#endif
	#if(NDIM>1)
	}
	#endif
	}

	return;
}

void init_scale() {

	int loc;

	for(ii=-1;ii<my_grid_dims[0];ii++) {
	#if(NDIM>1)
	for(jj=-1;jj<my_grid_dims[1];jj++) {
	#endif
	#if(NDIM==3)
	for(kk=-1;kk<my_grid_dims[2];kk++) {
	#endif
		for(loc=0;loc<=NDIM;loc++) {
			DLOOP {
				ND_ELEM(geom,ii,jj,kk).scale[loc][dd] = sqrt(ND_ELEM(geom,ii,jj,kk).gcov[loc][dd][dd]);
			}
		}
	#if(NDIM==3)
	}
	#endif
	#if(NDIM>1)
	}
	#endif
	}

	return;
}
