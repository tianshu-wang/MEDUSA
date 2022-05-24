#include "../decs.h"

double rtrans_solve(double dx, double rout) {

	double rs;
	double rg = 100*rout;

	// do a little bisection first
	rs = 0.1*dx;
	double fl;
    int iter = 0;
	do {
		rs *= 1.1;
		rtrans = rs;
		fl = r_of_x(n1*dx+startx[0])-rout;
        iter++;
	} while(iter < 100 && (isinf(fl) || isnan(fl)));
	rtrans = rg;
	double fr = r_of_x(n1*dx+startx[0])-rout;
	if(fl*fr > 0) {
		fprintf(stderr,"did not bracket rtrans in rtrans_solve %g->%g  %g->%g\n", rs, fl, rg, fr);
		exit(1);
	}
	do {
		rtrans = 0.5*(rs + rg);
		double fm = r_of_x(n1*dx+startx[0])-rout;
		if(fm*fl <= 0) {
			fr = fm;
			rg = rtrans;
		} else {
			fl = fm;
			rs = rtrans;
		}
	} while(fabs((rs-rg)/rg) > 0.1);
	rg = 0.5*(rs+rg);

	do {
		rs = rg;
		rtrans = rs;
		double f = r_of_x(n1*dx+startx[0]) - rout;
		double fp = sinh(n1*dx/rg) - n1*dx*cosh(n1*dx/rg)/rg;
		rg -= f/fp;
	} while(fabs((rs-rg)/rg) > 1.e-10);

	return fabs(rg);
}

int spider_factor1(int i) {

	double r0 = r_of_x(i*dx[0] + startx[0]);

	if(r0 > r_full_res) return 1;

	double r1 = r_of_x((i+1)*dx[0] + startx[0]);
	double dr = r1-r0;
	double rc = r_of_x((i+0.5)*dx[0] + startx[0]);

	double th_size = rc*dx[1];
	double dth_target = dr/rc;
	double nth_target = 2.*M_PI/dth_target;
	int nth=2;
	while(2*nth < nth_target && nth <= n2) nth *= 2;
	int sj = n2/nth;
	sj = MAX(sj,1);
	sj = MIN(sj,n2/2);

	return sj;
}

void ijk_to_x(int i, int j, int k, double x[SPACEDIM]){

	x[0] = (i+0.5)*dx[0] + startx[0];
	#if(NDIM>1)
    x[1] = (j-istart[1] + 0.5)*spider_fact[i]*dx[1] + startx[1];
	#else
	x[1] = 0.;
	#endif
	#if(NDIM==3)
    x[2] = (k+0.5)*dx[2] + startx[2];
	#else
	x[2] = 0.;
	#endif

	return;
}

void ijk_to_Cart(int i, int j, int k, double *x) {

	double r = r_of_x(i*dx[0] + startx[0]);
    double th = j*dx[1];
	double z = z_of_x(k*dx[2] + startx[2]);

	#if(NDIM==1)
	x[0] = r;
	#elif(NDIM==2)
	x[0] = r*cos(th);
	x[1] = r*sin(th);
	#elif(NDIM==3)
	x[0] = r*cos(th);
	x[1] = r*sin(th);
	x[2] = z;
	#endif

	return;
}


void init_spider() {

	spider_fact = (int *)malloc_rank1(n1+2*NG, sizeof(int)) + NG;
    avg_jfilt = (int *)malloc_rank1(n1, sizeof(int));
    avg_fact = (int **)malloc_rank1(n1, sizeof(int*));

    int imax_filter = 0;
	for(int i=0;i<n1+NG;i++) {
		#if(NDIM>1)
		spider_fact[i] = spider_factor1(i);
		//#if(NDIM==3)
		//spider_fact[i] = spider_factor2(i);
		//#endif
		if(i > 0) {
			if(spider_fact[i] < spider_fact[i-1]/2) spider_fact[i] = spider_fact[i-1]/2;
			if(spider_fact[i] > spider_fact[i-1]) spider_fact[i] = spider_fact[i-1];
		}
		if(i > 1) {
			if(spider_fact[i-2] != spider_fact[i-1]) spider_fact[i] = spider_fact[i-1];
		}
		#else
		spider_fact[i] = 1;
		#endif
        spider_fact[i] = 1;
        if(mpi_io_proc() && spider_fact[i] != 1) fprintf(stderr,"spider %d %d\n", i, spider_fact[i]);
	}
	for(int i=-NG;i<0;i++) {
		spider_fact[i] = spider_fact[0];
	}

	max_spider_level = 0;
	int s0 = spider_fact[0];
	while(s0 >>= 1) ++max_spider_level;
	if(mpi_io_proc()) fprintf(stderr,"max_spider_level = %d\n", max_spider_level);

    istop_spider = 0;
    for(int i=0;i<n1;i++) {
        if(spider_fact[i] == 1) {
            istop_spider = i;
            break;
        }
    }

	if(mpi_io_proc()) {
	    fprintf(stderr,"spider grid below cell %d at r = %g\n", istop_spider, r_of_x(startx[0]+istop_spider*dx[0]));
	}

	return;
}

// ************************************* //
//                                       //
//       basic coordinate mappings       //
//                                       //
// ************************************* //

double r_of_x(double x) {

	#if(RCOORD==UNIFORM)
	return fabs(x);
	#endif

	#if(RCOORD==SINH)
	return rtrans*sinh((x-startx[0])/rtrans) + startx[0];
	#endif
}

/*
double x_of_r(double r) {

	#if(RCOORD==UNIFORM)
	return r;
	#endif

	#if(RCOORD==SINH)
	double y = r/rtrans;
	return rtrans*log(y + sqrt(y*y + 1.));
	#endif
}
*/

double dr_dx(double x) {

	#if(RCOORD==UNIFORM)
	return 1;
	#endif

	#if(RCOORD==SINH)
	return cosh((x-startx[0])/rtrans);
	#endif
}

double d2r_dx2(double x) {

	#if(RCOORD==UNIFORM)
	return 0.;
	#endif

	#if(RCOORD==SINH)
	return sinh((x-startx[0])/rtrans)/rtrans;
	#endif
}

double z_of_x(double x) {

	#if(ZCOORD==UNIFORM)
	return x;
	#endif

	#if(ZCOORD==SINH)
	return ztrans*sinh(fabs(x)/ztrans);
	#endif
}

double x_of_z(double z) {

	#if(ZCOORD==UNIFORM)
	return z;
	#endif

	#if(ZCOORD==SINH)
	double y = z/ztrans;
	return ztrans*log(y + sqrt(y*y + 1.));
	#endif
}

double dz_dx(double x) {

	#if(ZCOORD==UNIFORM)
	return 1;
	#endif

	#if(ZCOORD==SINH)
	return cosh(fabs(x)/ztrans);
	#endif
}

double d2z_dx2(double x) {

	#if(ZCOORD==UNIFORM)
	return 0.;
	#endif

	#if(ZCOORD==SINH)
	return sinh(fabs(x)/ztrans)/ztrans;
	#endif
}


// ************************************* //
//                                       //
//     required coordinate integrals     //
//                                       //
// ************************************* //

// ALL DONE NUMERICALLY FOR CONVENIENCE AND GENERALITY

// radial integrals first

double rdr_dx(double x) {

	return r_of_x(x)*dr_dx(x);
}

double Vr(double x) {
	// \int_x^x+dx r dr/dx dx
	return romb(x,x+dx[0],rdr_dx);
}

double rd2r_dx2(double x) {

	return r_of_x(x)*d2r_dx2(x);
}

double int_rd2r_dx2(double x) {
	// \int_x^x+dx r^2 d^2r/dx^2 dx

	return romb(x,x+dx[0],rd2r_dx2);
}

double rdr_dx3(double x) {

	return r_of_x(x)*pow(dr_dx(x),3);
}

double int_rdr_dx3(double x) {
	// \int_x^x+dx r^2 (dr/dx)^3 dx

	return romb(x,x+dx[0],rdr_dx3);
}

double r2(double x) {
    return pow(r_of_x(x),2);
}

double int_r2(double x) {
    return romb(x,x+dx[0],r2);
}

double xr1g(double x) {
	return x*r_of_x(x)*dr_dx(x);
}

double xr2g(double x) {
	return x*x*r_of_x(x)*dr_dx(x);
}

double xr3g(double x) {
	return pow(x,3)*r_of_x(x)*dr_dx(x);
}

double xr4g(double x) {
	return pow(x,4)*r_of_x(x)*dr_dx(x);
}

// now do z integrals

double Vz(double x) {
	// \int_x^x+dx r dr/dx dx
	return romb(x,x+dx[2],dz_dx);
}

double int_d2z_dx2(double x) {
    return romb(x,x+dx[2],d2z_dx2);
}

double dz_dx3(double x) {
    return pow(dz_dx(x),3);
}

double int_dz_dx3(double x) {
    return romb(x,x+dx[2],dz_dx3);
}

double xz1g(double x) {
	return x*z_of_x(x)*dz_dx(x);
}

double xz2g(double x) {
	return x*x*z_of_x(x)*dz_dx(x);
}

double xz3g(double x) {
	return pow(x,3)*z_of_x(x)*dz_dx(x);
}

double xz4g(double x) {
	return pow(x,4)*z_of_x(x)*dz_dx(x);
}

// ************************************* //
//                                       //
//     finished coordinate integrals     //
//                                       //
// ************************************* //

/*
void x_to_rthphi(double x[SPACEDIM], double rvec[SPACEDIM]) {

	rvec[0] = r_of_x(x[0]);
	rvec[1] = ;
	rvec[2] = z_;
}

void ijk_to_rthphi(int i, int j, int k, double rvec[SPACEDIM]) {

	double x[SPACEDIM];

	ijk_to_x(i, j, k, x);

	rvec[0] = r_of_x(x[0]);
	rvec[1] = th_of_x(x[1]);
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
*/


double ijk_to_r(int i, int j, int k, double rhat[SPACEDIM]){

	double r,z,x[SPACEDIM];

	ijk_to_x(i, j, k, x);
	r = r_of_x(x[0]);
    z = z_of_x(x[2]);
	rhat[0] = dr_dx(x[0]);//sqrt(pow(dr_dx(x[0]),2.) + pow(dz_dx(x[2]),2.));
	rhat[1] = 0.;
	rhat[2] = 0.;

	return sqrt(r*r + z*z);
}


// ************************* //
// THIS NEEDS TO BE FIXED!!! //
// ************************* //
void vec_transform_to_xyz(double *vp, int i, int j, int k) {

	int l,m;
	double x[SPACEDIM];
	double v[SPACEDIM];
	double lam[SPACEDIM][SPACEDIM];
	double r,sth,cth,z;
	double drdx,dzdx;

	ijk_to_x(i, j, k, x);

	r = r_of_x(x[0]);
	drdx = dr_dx(x[0]);
    sth = sin(x[1]);
    cth = cos(x[1]);
    //z = z_of_x(x[2]);
    dzdx = dz_dx(x[2]);

	lam[0][0] = drdx*cth;
	lam[0][1] = -r*sth;
	lam[0][2] = 0;

	lam[1][0] = drdx*sth;
	lam[1][1] = r*cth;
	lam[1][2] = 0.;

	lam[2][0] = 0.;
	lam[2][1] = 0.;
	lam[2][2] = dzdx;

	v[0] = v[1] = v[2] = 0.;

	for(l=0;l<SPACEDIM;l++) {
		for(m=0;m<SPACEDIM;m++) {
			v[l] += lam[l][m]*vp[m];
		}
	}

	return;
}

void init_coords() {

    #if(NDIM<2)
    dx[1] = 2.*M_PI;
    #endif
    #if(NDIM<3)
    dx[2] = 1.;
    #endif

}


void init_interp_vol() {

	double x0,x1,r0,r1;

    #pragma omp parallel for private(x0,x1)
	for(int i=istart[0]-NG;i<istop[0]+NG;i++) {
		x0 = startx[0] + i*dx[0];
		x1 = x0 + dx[0];
		Gamma0[i] = Vr(x0);
		beta0[i] = romb(x0,x1,xr1g);
		alpha0[i] = romb(x0,x1,xr2g);
	}

	#if(NDIM>1)
    #pragma omp parallel for private(x0,x1)
	for(int j=istart[1]-NG; j<istop[1]+NG; j++) {
		x0 = startx[1] + j*dx[1];
		x1 = x0 + dx[1];
		Gamma1[j] = dx[1];
		beta1[j] = 0.5*(x1*x1 - x0*x0);
		alpha1[j] = (pow(x1,3.) - pow(x0,3.))/3.;
	}
    #endif

	#if(NDIM==3)
    #pragma omp parallel for private(x0,x1)
	for(int k=istart[2]-NG; k<istop[2]+NG; k++) {
		x0 = startx[2] + k*dx[2];
		x1 = x0 + dx[2];
		Gamma2[k] = Vz(x0);
		beta2[k] = romb(x0,x1,xz1g);
		alpha2[k] = romb(x0,x1,xz2g);
	}
	#endif

	return;
}

void calc_interp_coeffs1(int i, double *alpha, double *beta, double *Gamma, double *x) {

	if(spider_fact[i] == 1) {
		for(int j=istart[1]-NG; j<istop[1]+NG;j++) {
			alpha[j-istart[1]] = alpha1[j];
			beta[j-istart[1]] = beta1[j];
			Gamma[j-istart[1]] = Gamma1[j];
			x[j-istart[1]] = startx[1] + j*dx[1];
		}
	} else {
        for(int j=istart[1]-NG; j < istart[1]+my_grid_dims[1]/spider_fact[i]+NG; j++) {
            double x0 = startx[1] + istart[1]*dx[1] + (j-istart[1])*spider_fact[i]*dx[1];
            double x1 = x0 + spider_fact[i]*dx[1];
			Gamma[j-istart[1]] = x1-x0;
			beta[j-istart[1]] = 0.5*(x1*x1 - x0*x0);
			alpha[j-istart[1]] = (pow(x1,3)-pow(x0,3))/3.;
			x[j-istart[1]] = x0;
        }
	}
	return;
}

void calc_interp_coeffs2(int i, double *alpha, double *beta, double *Gamma, double *x) {

	for(int k=istart[2]-NG; k<istop[2]+NG; k++) {
		alpha[k-istart[2]] = alpha2[k];
		beta[k-istart[2]] = beta2[k];
		Gamma[k-istart[2]] = Gamma2[k];
		x[k-istart[2]] = startx[2] + k*dx[2];
	}
	return;
}


void init_volume() {

	int ii,jj,kk;

    #pragma omp parallel for private(ii,jj,kk) schedule(guided)
	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
    double rint = Vr(ii*dx[0] + startx[0]);
    #if(NDIM>1)
    for(jj=istart[1]-NG;jj<istart[1]+my_grid_dims[1]/spider_fact[ii]+NG;jj++) {
    #endif
    #if(NDIM==3)
    for(kk=istart[2]-NG;kk<istop[2];kk++) {
    double vzf = Vz(kk*dx[2]+startx[2]);
    #else
    double vzf = 1.;
    #endif
        ND_ELEM(geom,ii,jj,kk).volume = rint*vzf*dx[1]*spider_fact[ii];
    #if(NDIM==3)
    }
    #endif
    #if(NDIM>1)
    }
    #endif
    }

/*    double vtot = 0.;
    ZLOOP {
        vtot += ND_ELEM(geom,ii,jj,kk).volume;
    }

    fprintf(stderr,"total volume = %g   should be %g\n", vtot, M_PI*(pow(r_of_x(n1*dx[0]+startx[0]),2) - pow(r_of_x(startx[0]),2)));
*/

    return;
}


void init_area() {

	int ii,jj,kk;

    #pragma omp parallel for private(ii,jj,kk) schedule(guided)
	for(ii=istart[0]-1;ii<=istop[0]+1;ii++) {
		double rint = Vr(ii*dx[0]+startx[0]);
		#if(NDIM>1)
        for(jj=istart[1]-1; jj<=istart[1]+my_grid_dims[1]/spider_fact[ii]+1; jj++) {
		#endif
			#if(NDIM==3)
            for(kk=istart[2]-1; kk<=istop[2]+1; kk++) {
            double zint = Vz(kk*dx[2] + startx[2]);
			#else
            double zint = 1.;
            #endif

				ND_ELEM(geom,ii,jj,kk).area[0] = rdr_dx(ii*dx[0]+startx[0])*zint * dx[1]*spider_fact[ii];
				#if(NDIM>1)
				ND_ELEM(geom,ii,jj,kk).area[1] = rint * zint;
				#endif
				#if(NDIM==3)
				ND_ELEM(geom,ii,jj,kk).area[2] = rint * dx[1] * spider_fact[ii];
				#endif
			#if(NDIM==3)
			}
			#endif
		#if(NDIM>1)
		}
		#endif
	}

/*
	// make sure area on axis is a hard-zero
	#if(NDIM>1)
	if(startx[0] == 0. && istart[0] == 0) {
		for(jj=istart[1];jj<istart[1]+my_grid_dims[1]/spider_fact[0];jj++) {
		#if(NDIM==3)
		for(kk=istart[2];kk<istop[2];kk++) {
		#endif
			ND_ELEM(geom,0,jj,kk).area[0] = 0;
		#if(NDIM==3)
		}
		#endif
		}
	}
	#endif
*/


	return;
}


void init_conn() {

	int ii,jj,kk;

    #pragma omp parallel for private(ii,jj,kk) schedule(guided)
	for(ii=istart[0]-1;ii<istop[0]+1;ii++) {
		double xr0 = ii*dx[0] + startx[0];
		double vrf = Vr(xr0);
        //double ird2r_dx2 = int_rd2r_dx2(xr0);
        double ir2 = int_r2(xr0);
        //double rmid = r_of_x((ii+0.5)*dx[0] + startx[0]);
        //double drmid = dr_dx((ii+0.5)*dx[0] + startx[0]);
        //double d2rmid = d2r_dx2((ii+0.5)*dx[0] + startx[0]);
		#if(NDIM>1)
        for(jj=istart[1]-1; jj<istart[1]+my_grid_dims[1]/spider_fact[ii]+1; jj++) {
		#endif
		#if(NDIM==3)
        for(kk=istart[2]-1; kk<istop[2]+1; kk++) {
        double xz0 = kk*dx[0] + startx[0];
        double vzf = Vz(xz0);
		#else
        double xz0 = 0.;
        double vzf = 1.;
        #endif

		memset(&ND_ELEM(geom,ii,jj,kk).conn[0][0][0], 0, 27*sizeof(double));

		//ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = ird2r_dx2/vrf;
        ND_ELEM(geom,ii,jj,kk).conn[0][1][1] = -ir2/vrf;
        ND_ELEM(geom,ii,jj,kk).conn[2][2][2] = int_d2z_dx2(xz0)/vzf;

        // set for exact centrifugal balance
        ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = rdr_dx(beta0[ii]/Gamma0[ii])/ND_ELEM(geom,ii,jj,kk).gcov[0][1];
        //r_of_x((ii+0.5)*dx[0]+startx[0])*dr_dx((ii+0.5)*dx[0]+startx[0])/ND_ELEM(geom,ii,jj,kk).gcov[0][1];

        // set this to ensure exact angular momentum conservation
        ND_ELEM(geom,ii,jj,kk).conn[1][1][0] = -ND_ELEM(geom,ii,jj,kk).conn[0][1][1] * ND_ELEM(geom,ii,jj,kk).gcov[0][0]/ND_ELEM(geom,ii,jj,kk).gcov[0][1];


        //ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = ND_ELEM(geom,ii,jj,kk).conn[1][1][0];

        // set this to get \partial P\partial r = 0 correct
        if(spider_fact[ii+1] != spider_fact[ii]) {
            int j0 = 2*(jj-istart[1]) + istart[1];
            ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = (ND_ELEM(geom,ii+1,j0,kk).area[0]+ND_ELEM(geom,ii+1,j0+1,kk).area[0] - ND_ELEM(geom,ii,jj,kk).area[0])/ND_ELEM(geom,ii,jj,kk).volume - ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
        } else {
            ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = (ND_ELEM(geom,ii+1,jj,kk).area[0] - ND_ELEM(geom,ii,jj,kk).area[0])/ND_ELEM(geom,ii,jj,kk).volume - ND_ELEM(geom,ii,jj,kk).conn[1][0][1];
        }

/*
        ND_ELEM(geom,ii,jj,kk).conn[0][0][0] = d2rmid/drmid;
        ND_ELEM(geom,ii,jj,kk).conn[0][1][1] = -rmid/drmid;
        ND_ELEM(geom,ii,jj,kk).conn[1][0][1] = drmid/rmid;
        ND_ELEM(geom,ii,jj,kk).conn[1][1][0] = drmid/rmid;
        ND_ELEM(geom,ii,jj,kk).conn[2][2][2] = d2z_dx2((kk+0.5)*dx[2]+startx[2])/dz_dx((kk+0.5)*dx[2]+startx[2]);*/

		#if(NDIM==3)
		}
		#endif
		#if(NDIM>1)
		}
		#endif
	}

	return;
}

void init_gcov() {

	int i,j,k;

    #pragma omp parallel for private(i,j,k) schedule(guided)
	for(int i=istart[0]-NG; i<istop[0]+NG; i++) {
		double x0 = startx[0] + i*dx[0];
		double x1 = x0 + dx[0];
        double r0 = r_of_x(x0);
        double r1 = r_of_x(x1);
		double vrf = Vr(x0);
		double irdr_dx3 = int_rdr_dx3(x0);
        //double th0 = startx[1];
        //double th1 = th0 + dx[1];
        double dr0 = dr_dx(i*dx[0]+startx[0]);
        //double rmid = r_of_x((i+0.5)*dx[0] + startx[0]);
        //double drmid = dr_dx((i+0.5)*dx[0] + startx[0]);

		#if(NDIM>1)
        for(j=istart[1]-NG; j<istart[1]+my_grid_dims[1]/spider_fact[i] + NG; j++) {
		#endif
			#if(NDIM==3)
            for(k=istart[2]-NG; k<istop[2] + NG; k++) {
			#endif
				// initialize to zero
				for(int loc=0;loc<=NDIM;loc++) {
					for(int m=0;m<SPACEDIM;m++) {
						ND_ELEM(geom,i,j,k).gcov[loc][m] = 0.;
					}
				}

                double z0 = k*dx[2] + startx[2];
                double vzf = Vz(z0);
                double dzmid = dz_dx((k+0.5)*dx[2]+startx[2]);

				// volume average
				ND_ELEM(geom,i,j,k).gcov[0][0] = irdr_dx3/vrf;
				ND_ELEM(geom,i,j,k).gcov[0][1] = 0.25*(pow(r1,4) - pow(r0,4))/vrf;
                ND_ELEM(geom,i,j,k).gcov[0][2] = int_dz_dx3(z0)/vzf;

				// 0-direction -- face-averaged
				ND_ELEM(geom,i,j,k).gcov[1][0] = dr0*dr0;//pow(dr_dx(x0),2.);
				ND_ELEM(geom,i,j,k).gcov[1][1] = r0*r0;//pow(r0,2);
                ND_ELEM(geom,i,j,k).gcov[1][2] = ND_ELEM(geom,i,j,k).gcov[0][2];

				#if(NDIM>1)
				// 1-direction -- face-averaged; metric is axisymmetric!
				ND_ELEM(geom,i,j,k).gcov[2][0] = ND_ELEM(geom,i,j,k).gcov[0][0];
				ND_ELEM(geom,i,j,k).gcov[2][1] = ND_ELEM(geom,i,j,k).gcov[0][1];
				ND_ELEM(geom,i,j,k).gcov[2][2] = ND_ELEM(geom,i,j,k).gcov[0][2];
				#endif

				#if(NDIM==3)
				// 2-direction -- face-averaged
				ND_ELEM(geom,i,j,k).gcov[3][0] = ND_ELEM(geom,i,j,k).gcov[0][0];
				ND_ELEM(geom,i,j,k).gcov[3][1] = ND_ELEM(geom,i,j,k).gcov[0][1];
				ND_ELEM(geom,i,j,k).gcov[3][2] = pow(dz_dx(z0),2.);
				#endif
			#if(NDIM==3)
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

    #pragma omp parallel for private(ii,jj,kk,i) //collapse(NDIM+1)
	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
	#if(NDIM>1)
	for(jj=istart[1]-NG;jj<istart[1]+my_grid_dims[1]/spider_fact[ii] + NG;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
	#endif
	//	for(loc=0;loc<1+NDIM;loc++) {
		for(i=0;i<SPACEDIM;i++) {
			if(fabs(ND_ELEM(geom,ii,jj,kk).gcov[0][i]) > 1.e-50) {
				ND_ELEM(geom,ii,jj,kk).gcon[i] = 1./ND_ELEM(geom,ii,jj,kk).gcov[0][i];
			} else {
				ND_ELEM(geom,ii,jj,kk).gcon[i] = 0.;
			}
		}
	//	}
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

    #pragma omp parallel for private(ii,jj,kk,loc,dd) shared(my_grid_dims,spider_fact) //collapse(NDIM+2)
	for(ii=istart[0]-NG;ii<istop[0]+NG;ii++) {
	#if(NDIM>1)
	for(jj=istart[1]-NG;jj<istart[1]+my_grid_dims[1]/spider_fact[ii] + NG;jj++) {
	#endif
	#if(NDIM==3)
	for(kk=istart[2]-NG;kk<istop[2]+NG;kk++) {
	#endif
		for(loc=0;loc<=NDIM;loc++) {
			SLOOP {
				ND_ELEM(geom,ii,jj,kk).scale[loc][dd] = sqrt(ND_ELEM(geom,ii,jj,kk).gcov[loc][dd]);
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
