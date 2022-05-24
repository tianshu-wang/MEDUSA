
#include "../decs.h"

void ijk_to_x(int i, int j, int k, double x[SPACEDIM]){

	x[0] = (i+0.5)*dx[0] + startx[0];
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

	r = (i+0.5)*dx[0] + startx[0];
	//r /= (1. + dx[0]*dx[0]/(dx[0]*dx[0] + 12.*r*r));

	rhat[0] = 1.;
	rhat[1] = 0.;
	rhat[2] = 0.;

	return r;
}

double r_of_x(double x) {

	return x;
}

void x_to_rthphi(double x[SPACEDIM], double rvec[SPACEDIM]) {

	rvec[0] = r_of_x(x[0]);
	rvec[1] = x[1];
	rvec[2] = x[2];
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

void init_volume() {

	int ii,jj;
	double rin,rout,th0,th1;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		rin = ii*dx[0] + startx[0];
		rout = (ii+1)*dx[0] + startx[0];
		#if(NDIM==1)
		ND_ELEM(geom,ii,jj,kk).volume = 4./3.*M_PI*fabs(pow(rout,3.) - pow(rin,3.));
		#endif
		#if(NDIM==2)
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
			th0 = jj*dx[1] + startx[1];
			th1 = th0+dx[1];
			ND_ELEM(geom,ii,jj,kk).volume = 2.*M_PI/3.*fabs(pow(rout,3.)-pow(rin,3.))*fabs(cos(th1)-cos(th0));
		}
		#endif
	}

}

void init_coords() {

	int ii;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		redge[ii] = ii*dx[0]+startx[0];
		rcenter[ii] = (ii+0.5)*dx[0]+startx[0];
		dr[ii] = dx[0];//(ii+1)*dx[0]+startx[0] - redge[ii];
		//fprintf(stderr,"%d %g\n", ii, redge[ii]);
	}

}

void init_gdet() {

	int ii,jj,kk,loc;
	double x[SPACEDIM];
	double rtemp;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		#if(NDIM==1)
		ijk_to_x(ii,0,0,x);
		//ND_ELEM(geom,ii,jj,kk).gdet[0] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gdet[0] = (dx[0]*dx[0] + 3.*dx[0]*(ii*dx[0]+startx[0]) + 3.*pow(ii*dx[0]+startx[0],2.))/3.;
		x[0] -= 0.5*dx[0];
		ND_ELEM(geom,ii,jj,kk).gdet[1] = x[0]*x[0];
		#endif
		#if(NDIM>1)
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
			ijk_to_x(ii,jj,0,x);
			#if(NDIM==3)
			for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
			#endif
			ND_ELEM(geom,ii,jj,kk).gdet[0] = (x[0]*x[0] + dx[0]*dx[0]/12.)*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			x[0] -= 0.5*dx[0];
			ND_ELEM(geom,ii,jj,kk).gdet[1] = x[0]*x[0]*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			x[0] += 0.5*dx[0];
			x[1] -= 0.5*dx[1];
			ND_ELEM(geom,ii,jj,kk).gdet[2] = (x[0]*x[0] + dx[0]*dx[0]/12.)*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			#if(NDIM==3)
			ND_ELEM(geom,ii,jj,kk).gdet[3] = ND_ELEM(geom,ii,jj,kk).gdet[0];
			}
			#endif
		}
		#endif
	}

	//ND_ELEM(geom,-1,jj,kk).gdet[0] = ND_ELEM(geom,0,jj,kk).gdet[0];
	//ND_ELEM(geom,-1,jj,kk).gdet[1] = 0;

	return;
}

void init_conn() {

	int ii,jj,kk;
	double x[SPACEDIM],sth,cth,discrete_cth,r0;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		#if(NDIM>1)
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
		#endif
		#if(NDIM==3)
		for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
		#endif
		ijk_to_x(ii,jj,kk,x);
		r0 = ii*dx[0] + startx[0];
		memset(&ND_ELEM(geom,ii,jj,kk).conn[0][0][0], 0, 27*sizeof(double));
		sth = sin(x[1]);//0.5*(sin(x[1]) + sin(M_PI-x[1]));
		cth = cos(x[1]);//0.5*(cos(x[1]) - cos(x[1]+M_PI);
		discrete_cth = (sin(x[1]+0.5*dx[1])-sin(x[1]-0.5*dx[1]))/dx[1];
		ND_ELEM(geom,ii,jj,kk).conn[0][1][1] = -x[0];
		ND_ELEM(geom,ii,jj,kk).conn[0][2][2] = -x[0]*sth*sth;
		//ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = 1./x[0];
		ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = 3.*(dx[0]+2.*r0)/(2.*(dx[0]*dx[0] + 3.*dx[0]*r0+3.*r0*r0));
		ND_ELEM(geom,ii,jj,kk).conn[1][1][0] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];//1./x[0];
		ND_ELEM(geom,ii,jj,kk).conn[1][2][2] = -sth*cth;
		//ND_ELEM(geom,ii,jj,kk).conn[2][0][2] = 1./x[0];
		ND_ELEM(geom,ii,jj,kk).conn[2][0][2] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[2][1][2] = discrete_cth/sth;
		ND_ELEM(geom,ii,jj,kk).conn[2][2][0] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];//1./x[0];
		ND_ELEM(geom,ii,jj,kk).conn[2][2][1] = discrete_cth/sth;
		#if(NDIM==3)
		}
		#endif
		#if(NDIM>1)
		}
		#endif
	}

	//fprintf(stderr,"conn2 -> %g\n", ND_ELEM(geom,0,0,0).conn[1][0][1]);

	return;
}

void init_gcov() {

	double x[SPACEDIM];
	int ii,jj,kk,loc,i,j;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
	#if(NDIM>1)
	for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
	#endif
		for(loc=0;loc<=NDIM;loc++) {
			for(i=0;i<SPACEDIM;i++) {
				for(j=0;j<SPACEDIM;j++) {
					ND_ELEM(geom,ii,jj,kk).gcov[loc][i][j] = 0.;
				}
			}
		}
		ijk_to_x(ii,jj,kk,x);
		// centered
		ND_ELEM(geom,ii,jj,kk).gcov[0][0][0] = 1.;
		ND_ELEM(geom,ii,jj,kk).gcov[0][1][1] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[0][2][2] = pow(x[0]*sin(x[1]),2.);
		// - 0.5*dx[0]
		x[0] -= 0.5*dx[0];
		ND_ELEM(geom,ii,jj,kk).gcov[1][0][0] = 1.;
		ND_ELEM(geom,ii,jj,kk).gcov[1][1][1] = x[0]*x[0];
		ND_ELEM(geom,ii,jj,kk).gcov[1][2][2] = pow(x[0]*sin(x[1]),2.);
		// + 0.5*dx[1]
		x[0] += 0.5*dx[0];
		x[1] -= 0.5*dx[1];
		ND_ELEM(geom,ii,jj,kk).gcov[2][0][0] = 1.;
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
	#endif
	}

	return;
}

void init_gcon() {

	int ii,jj,kk,loc,i,j;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
	#if(NDIM>1)
	for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
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

	int ii,jj,kk,dd,loc;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
	#if(NDIM>1)
	for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
	#endif
		/*ND_ELEM(geom,ii,jj,kk).scale[0][0] = 1.;
		ND_ELEM(geom,ii,jj,kk).scale[0][1] = (ii+0.5)*dx[0]+startx[0];
		//ND_ELEM(geom,ii,jj,kk).scale[0][2]
		ND_ELEM(geom,ii,jj,kk).scale[1][0] = 1.;
		ND_ELEM(geom,ii,jj,kk).scale[2][1] = (ii+0.5)*dx[0]+startx[0];*/
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
