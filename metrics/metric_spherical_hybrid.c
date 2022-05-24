
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

double r_of_x(double x) {

	return rtrans*sinh(x/rtrans);
}

double x_of_r(double r) {

	double y = r/rtrans;

	return rtrans*log(y + sqrt(y*y + 1.));
}

double dr_dx(double x) {

	return cosh(x/rtrans);
}

void x_to_rthphi(double x[SPACEDIM], double rvec[SPACEDIM]) {

	rvec[0] = r_of_x(x[0]);
	rvec[1] = x[1];
	rvec[2] = x[2];
}	

void ijk_to_rthphi(int i, int j, int k, double rvec[SPACEDIM]) {

	double x[SPACEDIM];

	ijk_to_x(i, j, k, x);

	rvec[0] = r_of_x(x[0]);
	rvec[1] = x[1];
	rvec[2] = x[2];

	return;
}

void ijk_to_r_dr(int i, int j, int k, double r[]) {

	double x[SPACEDIM];

	ijk_to_x(i, j, k, x);

	r[0] = r_of_x(x[0]-0.5*dx[0]);
	r[1] = r_of_x(x[0]);
	r[2] = r_of_x(x[0]+0.5*dx[0]);
}


double ijk_to_r(int i, int j, int k, double rhat[SPACEDIM]){

	double r,x[SPACEDIM];

	ijk_to_x(i, j, k, x);
	/*r = redge[i]+0.5*dr[i];
	r *= (1. + 2.*dr[i]*dr[i]/(dr[i]*dr[i] + 12.*r*r));
	x[0] = x_of_r(r);*/
	r = rcenter[i];//r_of_x(x[0]);
	x[0] = x_of_r(r);
	rhat[0] = dr_dx(x[0]);
	rhat[1] = 0.;
	rhat[2] = 0.;

	return r;
}


void vec_transform_to_xyz(float *vp, int i, int j, int k) {

	int l,m;
	double x[SPACEDIM];
	float v[SPACEDIM];
	double lam[SPACEDIM][SPACEDIM];
	double r,sth,cth,sph,cph;

	ijk_to_x(i, j, k, x);

	r = r_of_x(x[0]);
	sth = sin(x[1]);
	cth = cos(x[1]);
	sph = sin(x[2]);
	cph = cos(x[2]);

	// x
	lam[0][0] = sth*cph;	// r
	lam[0][1] = cth*cph;	// theta
	lam[0][2] = -sth*sph;	// phi

	// y
	lam[1][0] = sth*sph;
	lam[1][1] = cth*sph;
	lam[1][2] = sth*cph;


	// z
	lam[2][0] = cth;
	lam[2][1] = -sth;
	lam[2][2] = 0.;

	v[0] = v[1] = v[2] = 0.;

	for(l=0;l<SPACEDIM;l++) {
		for(m=0;m<SPACEDIM;m++) {
			v[l] += lam[l][m]*vp[m]*ND_ELEM(geom,i,j,k).scale[0][m];
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

	int ii;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		redge[ii] = r_of_x(fabs(ii*dx[0]+startx[0]));
		//rcenter[ii] = fabs(r_of_x((ii+0.5)*dx[0]+startx[0]));
		dr[ii] = r_of_x(fabs((ii+1)*dx[0]+startx[0])) - redge[ii];
		rcenter[ii] = redge[ii]+0.5*dr[ii];//0.5*(redge[ii] + 0.5*dr[ii] + r_of_x((ii+0.5)*dx[0]+startx[0]));
		//fprintf(stderr,"%d %g\n", ii, redge[ii]);
		spider_fact[ii] = 1;
	}

}

void init_interp_vol() {

	int ii,jj,kk;
	double th0,th1;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		interp_vol[0][ii] = 4./3.*M_PI*fabs((pow(redge[ii]+dr[ii],3.) - pow(redge[ii],3.)));
	}
	if(startx[0] == 0.) {
		for(ii=istart[0]-NG;ii<istart[0];ii++) {
			interp_vol[0][ii] = interp_vol[0][istart[0] - ii - 1];
		}
	}

	#if(NDIM>1)
	for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
		th0 = jj*dx[1] + startx[1];
		th1 = th0 + dx[1];
		interp_vol[1][jj] = fabs(cos(th1)-cos(th0));
	}
	#endif

	#if(NDIM==3)
	for(kk=istart[2]-NG;kk<istop[3]+NG;k++) {
		interp_vol[2][kk] = 1.;
	}
	#endif

	return;
}

void init_volume() {

	int ii,jj,kk;
	double th0,th1;

	#if(NDIM==1)
	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		ND_ELEM(geom,ii,jj,kk).volume = 4./3.*M_PI*(pow(redge[ii]+dr[ii],3.) - pow(redge[ii],3.));
	}
	#elif(NDIM==2)
	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
			th0 = jj*dx[1] + startx[1];
			th1 = th0+dx[1];
			ND_ELEM(geom,ii,jj,kk).volume = fabs(2.*M_PI/3.*(pow(redge[ii]+dr[ii],3.) - pow(redge[ii],3.))*(cos(th1)-cos(th0)));
		}
	}
	#else
	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
			th0 = jj*dx[1] + startx[1];
			th1 = th0+dx[1];
			for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
				ND_ELEM(geom,ii,jj,kk).volume = dx[2]/3.*(pow(redge[ii]+dr[ii],3.) - pow(redge[ii],3.))*(cos(th1)-cos(th0));
			}
		}
	}
	#endif

	return;

}

void init_gdet() {

	int ii,jj,kk,loc;
	double x[SPACEDIM];
	double r,drdx;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		#if(NDIM==1)
		ijk_to_x(ii,0,0,x);
		r = r_of_x(fabs(x[0]));
		drdx = dr_dx(fabs(x[0]));
		//ND_ELEM(geom,ii,jj,kk).gdet[0] = r*r*drdx;
		ND_ELEM(geom,ii,jj,kk).gdet[0] = fabs(pow(rtrans,3.)*(pow(sinh(((ii+1)*dx[0]+startx[0])/rtrans),3.) - pow(sinh((ii*dx[0]+startx[0])/rtrans),3.))/(3.*dx[0]));
		r = r_of_x(fabs(x[0]-0.5*dx[0]));
		drdx = dr_dx(fabs(x[0]-0.5*dx[0]));
		ND_ELEM(geom,ii,jj,kk).gdet[1] = r*r*drdx;
		#endif
		#if(NDIM>1)
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
			ijk_to_x(ii,jj,0,x);			
			//ijk_to_r_dr(ii,jj,0,r);
			#if(NDIM==3)
			for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
			#endif
			r = r_of_x(fabs(x[0]));
			drdx = dr_dx(fabs(x[0]));
			ND_ELEM(geom,ii,jj,kk).gdet[0] = fabs(pow(rtrans,3.)*(pow(sinh(((ii+1)*dx[0]+startx[0])/rtrans),3.) - pow(sinh((ii*dx[0]+startx[0])/rtrans),3.))/(3.*dx[0])) * 0.5*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));//r*r*drdx*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			r = r_of_x(fabs(x[0]+0.5*dx[0]));
			drdx = dr_dx(fabs(x[0]+0.5*dx[0]));
			ND_ELEM(geom,ii,jj,kk).gdet[1] = r*r*drdx*0.5*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			r = r_of_x(fabs(x[0]));
			drdx = dr_dx(fabs(x[0]));
			x[1] += 0.5*dx[1];
			ND_ELEM(geom,ii,jj,kk).gdet[2] = fabs(pow(rtrans,3.)*(pow(sinh(((ii+1)*dx[0]+startx[0])/rtrans),3.) - pow(sinh((ii*dx[0]+startx[0])/rtrans),3.))/(3.*dx[0]))*0.5*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));//r*r*drdx*fabs(sin(x[1]));//0.5*x[0]*x[0]*(fabs(sin(x[1])) + fabs(sin(M_PI-x[1])));
			#if(NDIM==3)
			ND_ELEM(geom,ii,jj,kk).gdet[3] = ND_ELEM(geom,ii,jj,kk).gdet[0];
			}
			#endif
		}
		#endif
		//fprintf(stderr,"gdet: %d %g %g\n", ii, ND_ELEM(geom,ii,jj,kk).gdet[0], ND_ELEM(geom,ii,jj,kk).gdet[1]);
	}

	//ND_ELEM(geom,-1,jj,kk).gdet[0] = ND_ELEM(geom,0,jj,kk).gdet[0];
	//ND_ELEM(geom,-1,jj,kk).gdet[1] = 0;

	return;
}

void init_conn() {

	int ii,jj,kk;
	double x[SPACEDIM],sth,cth,discrete_cth;
	double x0,xp;

	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
		#if(NDIM>1)
		for(jj=istart[1]-NG;jj<istop[1]+NG;jj++) {
		#endif
		#if(NDIM==3)
		for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
		#endif
		ijk_to_x(ii,jj,kk,x);
		x0 = ii*dx[0];
		xp = (ii+1)*dx[0];
		memset(&ND_ELEM(geom,ii,jj,kk).conn[0][0][0], 0, 27*sizeof(double));
		sth = 0.5*(sin(x[1]) + sin(M_PI-x[1]));
		cth = 0.5*(cos(x[1]) - cos(x[1]+M_PI));
		discrete_cth = (0.5*(sin(x[1]+0.5*dx[1]) + sin(M_PI-x[1]-0.5*dx[1]))-0.5*(sin(x[1]-0.5*dx[1]) + sin(M_PI - x[1] + 0.5*dx[1])))/dx[1];
		/*ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = tanh(x[0]/rtrans)/rtrans;
		ND_ELEM(geom,ii,jj,kk).conn[0][1][1] = -rtrans*tanh(x[0]/rtrans);
		ND_ELEM(geom,ii,jj,kk).conn[0][2][2] = -sth*sth * rtrans*tanh(x[0]/rtrans);
		ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = 1./(tanh(x[0]/rtrans)*rtrans);
		ND_ELEM(geom,ii,jj,kk).conn[1][1][0] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[1][2][2] = -sth*cth;
		ND_ELEM(geom,ii,jj,kk).conn[2][0][2] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[2][1][2] = discrete_cth/sth;
		ND_ELEM(geom,ii,jj,kk).conn[2][2][0] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[2][2][1] = discrete_cth/sth;*/

		ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = (-9.*cosh(x0/rtrans) + cosh(3.*x0/rtrans) + 9.*cosh(xp/rtrans) - cosh(3.*xp/rtrans))/(4.*rtrans*(pow(sinh(x0/rtrans),3.)-pow(sinh(xp/rtrans),3.)));
		ND_ELEM(geom,ii,jj,kk).conn[0][1][1] = rtrans*(9.*cosh(x0/rtrans) - cosh(3.*x0/rtrans) - 9.*cosh(xp/rtrans) + cosh(3.*xp/rtrans))/(4.*(pow(sinh(x0/rtrans),3.)-pow(sinh(xp/rtrans),3.)));
		ND_ELEM(geom,ii,jj,kk).conn[0][2][2] = sth*sth*ND_ELEM(geom,ii,jj,kk).conn[0][1][1];
		ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = (pow(cosh(x0/rtrans),3.) - pow(cosh(xp/rtrans),3.))/(rtrans*(pow(sinh(x0/rtrans),3.) - pow(sinh(xp/rtrans),3.)));
		ND_ELEM(geom,ii,jj,kk).conn[1][1][0] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[1][2][2] = -sth*cth;
		ND_ELEM(geom,ii,jj,kk).conn[2][0][2] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[2][1][2] = discrete_cth/sth;
		ND_ELEM(geom,ii,jj,kk).conn[2][2][0] = ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
		ND_ELEM(geom,ii,jj,kk).conn[2][2][1] = discrete_cth/sth;
		//fprintf(stderr,"conn: %g %g %g %g %g %g %g %g %g %g\n", ND_ELEM(geom,ii,jj,kk).conn[0][0][0], ND_ELEM(geom,ii,jj,kk).conn[0][1][1], ND_ELEM(geom,ii,jj,kk).conn[0][2][2], ND_ELEM(geom,ii,jj,kk).conn[1][0][1], ND_ELEM(geom,ii,jj,kk).conn[1][1][0], ND_ELEM(geom,ii,jj,kk).conn[1][2][2], ND_ELEM(geom,ii,jj,kk).conn[2][0][2], ND_ELEM(geom,ii,jj,kk).conn[2][1][2], ND_ELEM(geom,ii,jj,kk).conn[2][2][0], ND_ELEM(geom,ii,jj,kk).conn[2][2][1] = discrete_cth/sth);

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

	double x[SPACEDIM],r[2],drdx[2];
	int loc,i,j,ii,jj,kk;

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
		r[0] = r_of_x(fabs(x[0]));
		drdx[0] = dr_dx(fabs(x[0]));
		r[1] = r_of_x(fabs(x[0]-0.5*dx[0]));
		drdx[1] = dr_dx(fabs(x[0]-0.5*dx[0]));
		
		// centered
		ND_ELEM(geom,ii,jj,kk).gcov[0][0][0] = drdx[0]*drdx[0];
		ND_ELEM(geom,ii,jj,kk).gcov[0][1][1] = r[0]*r[0];
		ND_ELEM(geom,ii,jj,kk).gcov[0][2][2] = pow(r[0]*sin(x[1]),2.);
		// + 0.5*dx[0]
		ND_ELEM(geom,ii,jj,kk).gcov[1][0][0] = drdx[1]*drdx[1];
		ND_ELEM(geom,ii,jj,kk).gcov[1][1][1] = r[1]*r[1];
		ND_ELEM(geom,ii,jj,kk).gcov[1][2][2] = pow(r[1]*sin(x[1]),2.);
		// + 0.5*dx[1]
		x[1] += 0.5*dx[1];
		ND_ELEM(geom,ii,jj,kk).gcov[2][0][0] = drdx[0]*drdx[0];
		ND_ELEM(geom,ii,jj,kk).gcov[2][1][1] = r[0]*r[0];
		ND_ELEM(geom,ii,jj,kk).gcov[2][2][2] = pow(r[0]*sin(x[1]),2.);
		
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
	//fprintf(stderr,"gcov: %d %g\n", ii, ND_ELEM(geom,ii,jj,kk).gcov[0][0][0]);
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
	//	for(loc=0;loc<1+NDIM;loc++) {
		for(i=0;i<SPACEDIM;i++) {
			for(j=0;j<SPACEDIM;j++) {
				if(fabs(ND_ELEM(geom,ii,jj,kk).gcov[0][i][j]) > 1.e-50) {
					ND_ELEM(geom,ii,jj,kk).gcon[i][j] = 1./ND_ELEM(geom,ii,jj,kk).gcov[0][i][j];
				} else {
					ND_ELEM(geom,ii,jj,kk).gcon[i][j] = 0.;
				}
			}
		}
	//	}
	#if(NDIM==3)
	}
	#endif
	#if(NDIM>1)
	}
	#endif
	//fprintf(stderr,"gcon: %d %g\n", ii, ND_ELEM(geom,ii,jj,kk).gcon[0][0]);
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
