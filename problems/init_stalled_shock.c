
#include "../decs.h"
#include "../constants.h"

void init_grid() {

	#if(0)
	startx[0] = log(5.e6);
	startx[1] = 0.;
	dx[0] = (log(1.e8)-log(5.e6))/N1;
	dx[1] = M_PI/N2;
	#endif

	#if(0)
	startx[0] = 5.e6;
	startx[1] = 0.;
	dx[0] = (1.e8-5.e6)/N1;
	dx[1] = M_PI/N2;
	#endif

	startx[0] = 1.e5;
	startx[1] = 0.;
	dx[0] = 5.e4;
	dx[1] = M_PI/N2;
	rtrans = 4.8e6;

	return;
}

static double f_freefall,Mdot,heat_constant,rshock,compression,entK;

void init_problem() {

	double x[SPACEDIM],r;
	double vr;
	void set_rshock();
	double find_vr(double r);

	gam = 4./3.;
	rho_floor = 1.e-6;
	e_floor = 1./(gam-1.);

	istep = 0;
 	t = 0.;

	dt_dump = 0.001;
	dt_restart = 1000.;
	dt_log = 0.01;
	dt_diag = 0.1;

	t_next_dump = dt_dump;
	t_next_restart = dt_restart;
	t_next_log = dt_log;
	t_next_diag = dt_diag;

	tmax = 0.5;
	dtmin = 1.e-8;

	bound_flag[0][0] = PROB;
	bound_flag[0][1] = PROB;
	bound_flag[1][0] = REFLECT;
	bound_flag[1][1] = REFLECT;

	M_prescribed = 1.5*MSUN;
	Mdot = -0.25*MSUN;
	f_freefall = 1.;///sqrt(2.);
	//rshock = 2.e7; 
	compression = 0.999*(gam + 1.)/(gam - 1.);
	//entK = exp(36.03);
	entK = exp(35.2);
	set_rshock();
	fprintf(stderr,"rshock = %g\n", rshock);

	ZLOOP {
		ijk_to_x(ii,jj,kk,x);
		r = r_of_x(x[0]);
		//r = exp((ii+0.5)*dx[0] + startx[0]);
		//vr = -sqrt(2.*GNEWT*M_prescribed/r);
		vr = find_vr(r);
		/*vr = -f_freefall * sqrt(2.*GNEWT*M_prescribed/r);
		if(r < rshock) {
			vr = find_vr(r);//-f_freefall * sqrt(2.*GNEWT*M_prescribed/rshock)/compression * (r/rshock);
		}*/
		NDP_ELEM(sim.p,ii,jj,kk,RHO) = Mdot/(4.*M_PI*r*r*vr);//*(1 + 0.001*(2.*rand() / (double)RAND_MAX - 1.));//*(1.+0.1*sin(8.*(jj+0.5)*dx[1]));
		NDP_ELEM(sim.p,ii,jj,kk,U1) = vr/ND_ELEM(geom,ii,jj,kk).scale[0][0];;//*(1 + 0.0*(2.*rand() / (double)RAND_MAX - 1.));
		NDP_ELEM(sim.p,ii,jj,kk,U2) = 0.;//vr*0.0*(2.*rand() / (double)RAND_MAX - 1.);
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;
		if(r < rshock) {
			NDP_ELEM(sim.p,ii,jj,kk,UU) = (GNEWT*M_prescribed/r - 0.5*vr*vr)*NDP_ELEM(sim.p,ii,jj,kk,RHO)/gam;
		} else {
			NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e22*pow(r/rshock,-3.*gam/2.);//1.e-2*0.5*NDP_ELEM(sim.p,ii,jj,kk,RHO)*vr*vr;//1.e-10*GNEWT*M_prescribed/r;
//entK*pow(Mdot/(4.*M_PI*rshock*rshock*(-f_freefall*sqrt(2.*GNEWT*M_prescribed/rshock)))*compression,gam)/(gam-1.) * ((gam+1.) - compression*(gam-1.))/(compression*(gam+1.)-(gam-1.)) * pow(r/rshock,-3.*gam/2.);
		}

	}

	heat_constant = 0.;//-2.27e3*NDP_ELEM(sim.p,0,0,0,UU)/NDP_ELEM(sim.p,0,0,0,RHO)*NDP_ELEM(sim.p,0,0,0,U1)/startx[0];
	fprintf(stderr,"heat_constant = %g\n", heat_constant);

	reset_boundaries(sim.p);
	update_eos(sim.p, sim.eos);
	//exit(1);	
	return;
}

void prob_bounds(int i, int j, int k, double *p) {

	double r,vr;
	double find_vr(double r);

	r = exp((i+0.5)*dx[0] + startx[0]);
	if(r > 5.e7){//rshock) {
		vr = -f_freefall * sqrt(2.*GNEWT*M_prescribed/r);
		p[RHO] = Mdot/(4.*M_PI*r*r*vr);
		p[U1] = vr/ND_ELEM(geom,i,j,k).scale[0][0];;
		p[UU] = 1.e22*pow(r/rshock,-3.*gam/2.);//1.e-2*0.5*p[RHO]*vr*vr;
	} else {
		vr = find_vr(r);
		p[RHO] = Mdot/(4.*M_PI*r*r*vr);
		p[U1] = vr/ND_ELEM(geom,i,j,k).scale[0][0];;
		p[UU] = (GNEWT*M_prescribed/r - 0.5*vr*vr)*p[RHO]/gam;
		//fprintf(stderr,"%g %g\n", r, vr);
	}
	p[U2] = 0.;
	p[U3] = 0.;


	return;
}

void ext_src(double NDP_PTR p) {

	double r,rhat[SPACEDIM];
	double ijk_to_r(int i, int j, int k, double rhat[]);

	ZLOOP {
		r = ijk_to_r(ii,jj,kk,rhat);
		if(r < rshock) {
			NDP_ELEM(sim.src,ii,jj,kk,UU) += heat_constant*NDP_ELEM(p,ii,jj,kk,RHO)*(startx[0]*startx[0]/(r*r));
		}
	}

}

double find_vr_poly(double r) {

	int cnt;
	double vr,cst,f,df,delta,vrinit;
	double vrtop,vrbot,ftop,fbot,fmid;

	vr = f_freefall*sqrt(2.*GNEWT*M_prescribed/rshock)/compression;

	cst = 0.;//0.5*vr*vr + gam/(gam-1.)*entK*pow(-Mdot/(4.*M_PI*rshock*rshock*vr),gam-1.) - GNEWT*M_prescribed/rshock;

	//vr *= pow(r/rshock,2.);

	vr = -Mdot/(4.*M_PI*r*r)*pow(GNEWT*M_prescribed/r * (gam-1.)/(entK*gam), -1./(gam-1.));

	vrinit = vr;

	cnt = 0;
	do {
		f = 0.5*vr*vr + gam/(gam-1.)*entK*pow(-Mdot/(4.*M_PI*r*r*vr),gam-1.) - GNEWT*M_prescribed/r - cst;
		df = vr - gam*entK*pow(-Mdot/(4.*M_PI*r*r),gam-1.)*pow(vr,-gam);
		delta = - f/df;

		//fprintf(stderr,"vr,f,df,delta: %d %g %g %g %g %g\n", cnt, r, vr, f, df, delta);
		if(isnan(f) || isnan(delta)) {
			//cnt = 100;
			//break;
			exit(1);
		}

		vr += delta;
		cnt++;
	}while(cnt < 100 && (fabs(delta/vr) > 1.e-12 || fabs(f) > 1.e-12));
	if(cnt==100) vr = 0/0;

	//exit(1);

	return -vr;
}

double find_vr(double r) {

	static int firstc=1;
	int cnt;
	double vr,cst,f,df,delta,vr0,vrinit;
	double vrtop,vrbot,ftop,fbot,fmid;
	double r0,rhat[SPACEDIM];

	r0 = startx[0];//ijk_to_r(0,0,0,rhat);

	vr0 = -1.2e6;

	cst = 2.*((2.*gam-3.)*log(r0) + (gam-1.)*log((5.-3.*gam)*fabs(vr0)) + log(2.*GNEWT*M_prescribed - r0*vr0*vr0))/(3.*gam-5.);

	//vr = vrinit = vr0 * r/startx[0];

	if(r < rshock) {
		vr = vr0 * r/r0;
	} else {
		vr = -sqrt(2.*GNEWT*M_prescribed/r);
	}

	vrbot = 5.*vr;
	vrtop = 1.;

	if(vrbot < -sqrt(2.*GNEWT*M_prescribed/r)) {
		vrbot = -sqrt(2.*GNEWT*M_prescribed/r);
		fbot = 2.*((2.*gam-3.)*log(r) + (gam-1.)*log((5.-3.*gam)*fabs(vrbot)) + log(2.*GNEWT*M_prescribed - r*vrbot*vrbot))/(3.*gam-5.) - cst;
		//fprintf(stdout,"%g %g %d %g %g %g %g\n", r, vrbot, 0, fbot, 2.*(2.*gam-3.)*log(r)/(3.*gam-5), 2.*(gam-1.)*log((5.-3.*gam)*fabs(vrbot))/(3.*gam-5.),cst);
		return vrbot;
	}

/*	fbot = 2.*((2.*gam-3.)*log(r) + (gam-1.)*log((5.-3.*gam)*fabs(vrbot)) + log(fabs(2.*GNEWT*M_prescribed - r*vrbot*vrbot)))/(3.*gam-5.) - cst;
	ftop = 2.*((2.*gam-3.)*log(r) + (gam-1.)*log((5.-3.*gam)*fabs(vrtop)) + log(fabs(2.*GNEWT*M_prescribed - r*vrtop*vrtop)))/(3.*gam-5.) - cst;

	if(fbot*ftop > 0) {
		fprintf(stderr,"oops! %g %g %g\n", r/rshock, fbot,ftop);
		exit(1);
	}

	//fprintf(stderr,"cst = %g %g %g %g\n", cst, fbot, ftop, r);

	cnt = 0;
	do {
		vr = 0.5*(vrbot+vrtop);
		fmid = 2.*((2.*gam-3.)*log(r) + (gam-1.)*log((5.-3.*gam)*fabs(vr)) + log((2.*GNEWT*M_prescribed - r*vr*vr)))/(3.*gam-5.) - cst;
		if(fmid*fbot < 0) {
			ftop = fmid;
			vrtop = vr;
		} else {
			fbot = fmid;
			vrbot = vr;
		}
		cnt++;
	} while(fabs((vrtop-vrbot)/vr) > 1.e-7);
	vr = 0.5*(vrtop+vrbot);
*/

	vr = vr0 * r/rshock;
	cnt = 0;
	do {
		f = 2.*((2.*gam-3.)*log(r) + (gam-1.)*log((5.-3.*gam)*fabs(vr)) + log(2.*GNEWT*M_prescribed - r*vr*vr))/(3.*gam-5.) - cst;
		df = 2.*((gam-1.)/vr - 2.*r*vr/(2.*GNEWT*M_prescribed - r*vr*vr))/(3.*gam-5.);
		delta = -f/df;
		vr += delta;
		//fprintf(stderr,"%g %g %g %g\n", vr, f, df, delta);
		cnt++;
	}while(cnt < 100 && (fabs(delta/vr) > 1.e-12 || fabs(f) > 1.e-12));

	//fprintf(stdout,"%g %g %d %g\n", r, vr, cnt, fmid);

	return vr;
}


void set_rshock() {

	int i,cnt;
	double f,df,delta,cst,vr0,vff;
	double r0,rhat[SPACEDIM];

	vr0 = -1.2e6;

//	vff = -sqrt(2.*GNEWT*M_prescribed/startx[0]);

//	for(i=0;i<1000000;i++) {

//		vr0 = vff * (i+0.5)/1000000.;
	r0 = startx[0];//ijk_to_r(0,0,0,rhat);
	fprintf(stderr,"r0 = %g\n", r0);

		cst = 2.*((2.*gam-3.)*log(r0) + (gam-1.)*log((5.-3.*gam)*fabs(vr0)) + log(2.*GNEWT*M_prescribed - r0*vr0*vr0))/(3.*gam-5.);

		rshock = 2.e7;

		cnt = 0;
		do {
			f = 2.*((2.*gam-3.)*log(rshock) + (gam-1.)*log((5.-3.*gam)*sqrt(2.*GNEWT*M_prescribed/rshock)/compression) + log(2.*GNEWT*M_prescribed - 2.*GNEWT*M_prescribed/(compression*compression)))/(3.*gam-5.) - cst;
			df = 1./rshock;
			delta = -f/df;
			rshock += delta;
			cnt++;
		} while(cnt < 100 && (fabs(delta/rshock) > 1.e-12 || fabs(f) > 1.e-12));

//		fprintf(stdout,"%g %g %d\n", vr0, rshock, (int)!(!(cnt < 100) || isnan(rshock)));
//	}
//	exit(1);
	return;
}

void set_rshock_poly() {

	int cnt;
	double f,df,delta;
	rshock = 2.e7;

	cnt = 0;
	do {
		f = GNEWT*M_prescribed/rshock*(pow(f_freefall/compression,2.) - 1.)
			+ gam/(gam-1.)*entK*pow(-Mdot/(4.*M_PI*f_freefall/compression*sqrt(2.*GNEWT*M_prescribed)*pow(rshock,1.5)),gam-1.);
		df = GNEWT*M_prescribed/(rshock*rshock) * (1. - pow(f_freefall/compression,2.))
			- 1.5*gam*entK*pow(-Mdot/(4.*M_PI*f_freefall/compression*sqrt(2.*GNEWT*M_prescribed)), gam-1.)*pow(rshock,-1.5*gam+0.5);
		delta = -f/df;

		fprintf(stderr,"on iter %d  rs = %g   %g %g\n", cnt, rshock, f, df);
		cnt++;
		rshock += delta;
	} while(fabs(delta/rshock > 1.e-10));

	fprintf(stderr,"rshock converged at %g\n", rshock);

//	exit(1);

}
