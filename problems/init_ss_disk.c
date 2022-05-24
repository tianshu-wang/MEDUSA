#include "../decs.h"
#include "../constants.h"

double alpha,Mdot;

void init_grid() {

    M_prescribed = 10*MSUN;
    Rschw = 2*GNEWT*M_prescribed/(CLIGHT*CLIGHT);

	startx[0] = 2.7/2.0*Rschw; 
	startx[1] = 0;
	startx[2] = -0.2*startx[0];

    double rout = 100*Rschw;
    dx[0] = (rout-startx[0])/n1;
    dx[0] *= 0.5;   // if using sinh grid
	dx[1] = 2.*M_PI/n2;
	dx[2] = -2*startx[2]/n3;

    fprintf(stderr,"%g %g %g\n", startx[0], rout, dx[0]);

	rtrans = rtrans_solve(dx[0], rout);
	if(mpi_io_proc()) fprintf(stderr,"rtrans = %g   rout = %g\n", rtrans, r_of_x(n1*dx[0] + startx[0])/Rschw);
    fflush(stderr);

	r_full_res = 1e99;

	periodic[0] = 0;
	periodic[1] = 1;
	periodic[2] = 0;	

	return;
}

double get_kappa(double R) {
    return SIGMA_THOMSON/MP;
}

double get_f(double R) {
    double RISCO = 6*GNEWT*M_prescribed/(CLIGHT*CLIGHT);
    return 1- sqrt(RISCO/R);
}

double get_Omega(double R) {
    #if(PN_POTENTIAL==TRUE)
    return sqrt(GNEWT*M_prescribed/(R*pow(R-Rschw,2)));
    #else
    return sqrt(GNEWT*M_prescribed/pow(R,3));
    #endif
}

double get_T(double R) {
    return pow((3.*CLIGHT*get_Omega(R)/(alpha*ARAD*get_kappa(R))),1./4.);
}

double get_Sigma(double R) {
    return 8*M_PI*CLIGHT*CLIGHT/(Mdot*alpha*pow(get_kappa(R),2) * get_Omega(R));// /get_f(R);
}

double get_vr(double R) {
    return pow(Mdot,2)*alpha*pow(get_kappa(R),2)*get_Omega(R)/(16*M_PI*R*pow(CLIGHT,2));//*get_f(R);
}

double get_vphi(double R) {
    return sqrt(GNEWT*M_prescribed/R);
}

double get_P(double R) {
    return ARAD*pow(get_T(R),4)/3.;
}

double get_H(double R) {
    return 2.*get_P(R)/(get_Sigma(R)*pow(get_Omega(R),2));
}

double get_rho(double R) {
    return get_Sigma(R)/(2.*get_H(R));
}

double get_Qplus(double R, double T) {
    double P = ARAD*pow(T,4.)/3.;
    return 1.5*alpha*P*get_Omega(R);
}

double get_Qminus(double R, double T) {
    return ARAD*CLIGHT*pow(T,4)/(get_kappa(R)*get_H(R)*get_Sigma(R));
}

void init_problem() {

	int ii,jj,kk,vv;
	double xl,xr;
	double rho_amb,e_amb,u_amb;
	double E_expl,R_expl;

	gam = 5./3.;

	dtmin = 1.e-14;

	VLOOP {
		bc[vv].lo[0] = OUTFLOW;//DISK;
		#if(NDIM>1)
		bc[vv].lo[1] = PERIODIC;
		bc[vv].hi[1] = PERIODIC;
		#endif
		#if(NDIM==3)
		bc[vv].lo[2] = OUTFLOW;
		bc[vv].hi[2] = OUTFLOW;
		#endif
	}
    //bc[U2].lo[0] = EXTRAP;
	HLOOP {
		bc[vv].hi[0] = PROB;
	}

    double Ledd = 4*M_PI*GNEWT*M_prescribed*MP*CLIGHT/SIGMA_THOMSON;
    double Mdotedd = Ledd/(0.1*CLIGHT*CLIGHT);

    alpha = 1.e-2;
    Mdot = 0.1*Mdotedd;


    if(restart_from_hdf) {
        restart_read();
        return;
    }

    rho_floor = get_Sigma(4*Rschw);
    e_floor = rho_floor*KBOLTZ*get_T(4*Rschw)/((gam-1)*MP);

	ZLOOP {
        double rcyl = r_of_x(beta0[ii]/Gamma0[ii]);
        double th = (jj-istart[1] + 0.5)*spider_fact[ii]*dx[1] + istart[1]*dx[1] + startx[1];
        double z = z_of_x((kk+0.5)*dx[2] + startx[2]);
        NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.e5*exp(-pow(rcyl-50*Rschw,2)/(900.*Rschw*Rschw));
        NDP_ELEM(sim.p,ii,jj,kk,UU) = 1.e24;//NDP_ELEM(sim.p,ii,jj,kk,RHO)*KBOLTZ*1e4/((gam-1)*MP);
        NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;
    	//	NDP_ELEM(sim.p,ii,jj,kk,RHO) = get_Sigma(rcyl);//*(1 + 1.e-12*(jj%2+1));
            //NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.e-2 + exp(-(pow(rcyl-30*startx[0],2))/pow(3*startx[0],2.));
		//    NDP_ELEM(sim.p,ii,jj,kk,UU) = NDP_ELEM(sim.p,ii,jj,kk,RHO)*KBOLTZ*get_T(rcyl)/((gam-1)*MP);
    	//	NDP_ELEM(sim.p,ii,jj,kk,U1) = -get_vr(rcyl)/dr_dx((ii+0.5)*dx[0]+startx[0]);
	    	NDP_ELEM(sim.p,ii,jj,kk,U2) = get_Omega(rcyl);
        if(rcyl < 3*Rschw) {
            //NDP_ELEM(sim.p,ii,jj,kk,RHO) = get_Sigma(6*Rschw -rcyl);;
            //NDP_ELEM(sim.p,ii,jj,kk,UU) =  NDP_ELEM(sim.p,ii,jj,kk,RHO)*KBOLTZ*get_T(6*Rschw-rcyl)/((gam-1)*MP);
            //NDP_ELEM(sim.p,ii,jj,kk,U1) = 0.;//-1.5*sqrt((gam-1)*NDP_ELEM(sim.p,ii,jj,kk,UU)/NDP_ELEM(sim.p,ii,jj,kk,RHO))*(3*Rschw - rcyl)/Rschw/dr_dx((ii+0.5)*dx[0]+startx[0]);
            //NDP_ELEM(sim.p,ii,jj,kk,U2) = get_Omega(3*Rschw)*(9*Rschw*Rschw)*ND_ELEM(geom,ii,jj,kk).gcon[1];
        } /*else { 
    		NDP_ELEM(sim.p,ii,jj,kk,RHO) = get_Sigma(rcyl);//*(1 + 1.e-12*(jj%2+1));
            //NDP_ELEM(sim.p,ii,jj,kk,RHO) = 1.e-2 + exp(-(pow(rcyl-30*startx[0],2))/pow(3*startx[0],2.));
		    NDP_ELEM(sim.p,ii,jj,kk,UU) = NDP_ELEM(sim.p,ii,jj,kk,RHO)*KBOLTZ*get_T(rcyl)/((gam-1)*MP);
    		NDP_ELEM(sim.p,ii,jj,kk,U1) = -get_vr(rcyl)/dr_dx((ii+0.5)*dx[0]+startx[0]);
	    	NDP_ELEM(sim.p,ii,jj,kk,U2) = get_Omega(rcyl);
        }*/
        //NDP_ELEM(sim.p,ii,jj,kk,U2) = 30*startx[0]*cos(th)/ND_ELEM(geom,ii,jj,kk).scale[0][1] * exp(-(pow(rcyl-30*startx[0],2))/pow(3*startx[0],2.));
        //if(ii==0 && jj==0) fprintf(stderr,"torb inner = %g\n", 2*M_PI/get_Omega(rcyl));
		NDP_ELEM(sim.p,ii,jj,kk,U3) = 0.;//10.;
        NDP_ELEM(sim.p,ii,jj,kk,YE) = (th < M_PI ? 1. : -1.);//sin(th);
        //NDP_ELEM(sim.p,ii,jj,kk,YE) = (th > 0.5*M_PI && th < 1.5*M_PI ? 1. : 0.);
    }



	return;
}

void prob_bounds(int i, int j, int k, double *p) {
    if(half_step_sync) {
    p[RHO] = NDP_ELEM(sim.ph,istop[0]-1,j,k,RHO);
    p[UU] = NDP_ELEM(sim.ph,istop[0]-1,j,k,UU);
    p[U1] = NDP_ELEM(sim.ph,istop[0]-1,j,k,U1);
    p[U2] = get_Omega(r_of_x(beta0[i]/Gamma0[i]));
    p[U3] = NDP_ELEM(sim.ph,istop[0]-1,j,k,U3);
    p[YE] = NDP_ELEM(sim.ph,istop[0]-1,j,k,YE);
    } else {
    p[RHO] = NDP_ELEM(sim.p,istop[0]-1,j,k,RHO);
    p[UU] = NDP_ELEM(sim.p,istop[0]-1,j,k,UU);
    p[U1] = NDP_ELEM(sim.p,istop[0]-1,j,k,U1);
    p[U2] = get_Omega(r_of_x(beta0[i]/Gamma0[i]));
    p[U3] = NDP_ELEM(sim.p,istop[0]-1,j,k,U3);
    p[YE] = NDP_ELEM(sim.p,istop[0]-1,j,k,YE);
    }
	return;
}


void ext_src(double NDP_PTR p, double dtpush) {

    static int firstc = 1;
    static double NDP_PTR Tvisc;
    int ii,jj,kk,jp,jm;
    double gradv[SPACEDIM][SPACEDIM],inv_gradv[SPACEDIM][SPACEDIM];
    double esrc,src[SPACEDIM],Tm,Tp;
    void gridcell_vec_malloc(double NDP_PTR* p, int vlen);

    if(firstc) {
        gridcell_vec_malloc(&Tvisc, SPACEDIM*SPACEDIM);
        firstc = 0;
    }

    for(ii=istart[0]-1; ii < istop[0]+1; ii++) {
    #if(NDIM>1)
    for(jj=istart[1]-1; jj < istart[1]+my_grid_dims[1]/spider_fact[ii]+1; jj++) {
    #endif
        memset(gradv, 0, SPACEDIM*SPACEDIM*sizeof(double));

        // calculate gradv
        #if(NDIM>1)
        if(spider_fact[ii+1] != spider_fact[ii]) {
            jp = (jj-istart[1])*2 + istart[1];
        }
        if(spider_fact[ii-1] != spider_fact[ii]) {
            jm = (jj-istart[1])/2 + istart[1];
        }
        #endif
        for(int i=0;i<SPACEDIM;i++) {
                gradv[i][0] = (NDP_ELEM(p,ii+1,jp,kk,U1+i) - NDP_ELEM(p,ii-1,jm,kk,U1+i))/(2*dx[0]);
                //if(ii>1021 && i==1) {
                //    fprintf(stderr,"%d gradv[%d][0] = %g   %g %g\n", ii, i, gradv[i][0], NDP_ELEM(p,ii-1,jm,kk,U1+i), NDP_ELEM(p,ii+1,jp,kk,U1+i));
                //}
        }
        #if(NDIM>1)
        for(int i=0;i<SPACEDIM;i++) {
                gradv[i][1] = (NDP_ELEM(p,ii,jj+1,kk,U1+i) - NDP_ELEM(p,ii,jj-1,kk,U1+i))/(2*dx[1]*spider_fact[ii]);
        }
        #endif
        double divv = 0.;
        for(int i=0;i<SPACEDIM;i++) {
            for(int j=0;j<NDIM;j++) {
                for(int k=0;k<SPACEDIM;k++) {
                    gradv[i][j] += ND_ELEM(geom,ii,jj,kk).conn[i][j][k]*NDP_ELEM(p,ii,jj,kk,U1+k);
                }
            }
            divv += gradv[i][i];
        }

        for(int i=0;i<SPACEDIM;i++) {
            for(int j=0;j<SPACEDIM;j++) {
                inv_gradv[i][j] = gradv[j][i]*ND_ELEM(geom,ii,jj,kk).gcon[i]*ND_ELEM(geom,ii,jj,kk).gcov[0][j];
                //inv_gradv[i][j] = gradv[i][j]*ND_ELEM(geom,ii,jj,kk).gcov[0][i];
            }
        }

        double rcyl = r_of_x(beta0[ii]/Gamma0[ii]);
        double mus = alpha*pow(NDP_ELEM(sim.eos,ii,jj,kk,CS),2)/get_Omega(rcyl);
        mus *= NDP_ELEM(p,ii,jj,kk,RHO);
        //double Prad = ARAD*pow(NDP_ELEM(sim.eos,ii,jj,kk,TEMP),4)/3.;
        //double Omega = get_Omega(rcyl);
        //double Hdisk = 0.;//NDP_ELEM(sim.eos,ii,jj,kk,CS)/Omega;
        //double mus = -alpha*(NDP_ELEM(sim.eos,ii,jj,kk,PRESS) + 2*Hdisk*Prad)/Omega;
        //double mus = alpha*NDP_ELEM(sim.eos,ii,jj,kk,PRESS)/get_Omega(rcyl);

        for(int i = 0; i < SPACEDIM; i++) {
            for(int j = 0; j < SPACEDIM; j++) {
                //NDP_ELEM(Tvisc,ii,jj,kk,i*SPACEDIM+j) = 0.5*mus*(inv_gradv[i][j] + inv_gradv[j][i] - 2./3.*ND_ELEM(geom,ii,jj,kk).gcov[0][i]*divv)*ND_ELEM(geom,ii,jj,kk).gcon[i];
                NDP_ELEM(Tvisc,ii,jj,kk,i*SPACEDIM+j) = 0.5*mus*(gradv[i][j] + inv_gradv[i][j] - 2./3.*DELTA(i,j)*divv);
                //if(ii > 1021 && i == 0 && j==1) {
                //    fprintf(stderr,"%d Tvisc[%d][%d] = %g  %g %g %g\n", ii, i, j, NDP_ELEM(Tvisc,ii,jj,kk,i*SPACEDIM+j), mus, inv_gradv[i][j], inv_gradv[j][i]);
                //}
            }
        }
        //if(p == sim.ph) {
        //fprintf(stderr,"%d %g %g %g\n", ii, gradv[0][1], inv_gradv[0][1], gradv[1][0]);
        //}
    #if(NDIM>1)
    }
    #endif
    }

    //if(p == sim.ph) {
    //fflush(stderr);
    //exit(2);
    //}

    ZLOOP {
        #if(NDIM>1)
        if(spider_fact[ii+1] != spider_fact[ii]) {
            jp = (jj-istart[1])*2 + istart[1];
        }
        if(spider_fact[ii-1] != spider_fact[ii]) {
            jm = (jj-istart[1])/2 + istart[1];
        }
        #endif

        esrc = 0.;
        for(int i=0;i<SPACEDIM;i++) src[i] = 0.;

        double vdotTm = 0.;
        double vdotTp = 0.;
        for(int i=0;i<SPACEDIM;i++) {
            Tm = 0.5*(NDP_ELEM(Tvisc,ii,jj,kk,i) + NDP_ELEM(Tvisc,ii-1,jm,kk,i));
            vdotTm += Tm*NDP_ELEM(sim.vedge[0],ii,jj,kk,i)*ND_ELEM(geom,ii,jj,kk).area[0];
            if(spider_fact[ii+1] != spider_fact[ii]) {
                Tp = 0.5*(NDP_ELEM(Tvisc,ii,jj,kk,i) + NDP_ELEM(Tvisc,ii+1,jp,kk,i));
                double Tps = 0.5*(NDP_ELEM(Tvisc,ii,jj,kk,i) + NDP_ELEM(Tvisc,ii+1,jp+1,kk,i));
                vdotTp += Tp*NDP_ELEM(sim.vedge[0],ii+1,jp,kk,i)*ND_ELEM(geom,ii+1,jp,kk).area[0] + Tps*NDP_ELEM(sim.vedge[0],ii+1,jp+1,kk,i)*ND_ELEM(geom,ii+1,jp+1,kk).area[0];
                src[i] -= (Tp*ND_ELEM(geom,ii+1,jp,kk).area[0]+Tps*ND_ELEM(geom,ii+1,jp+1,kk).area[0] - Tm*ND_ELEM(geom,ii,jj,kk).area[0])/ND_ELEM(geom,ii,jj,kk).volume;
            } else {
                Tp = 0.5*(NDP_ELEM(Tvisc,ii,jj,kk,i) + NDP_ELEM(Tvisc,ii+1,jp,kk,i));
                vdotTp += Tp*NDP_ELEM(sim.vedge[0],ii+1,jp,kk,i)*ND_ELEM(geom,ii+1,jp,kk).area[0];
                src[i] -= (Tp*ND_ELEM(geom,ii+1,jj,kk).area[0] - Tm*ND_ELEM(geom,ii,jj,kk).area[0])/ND_ELEM(geom,ii,jj,kk).volume;
                //if(ii > 1020 && i==1) {
                //    fprintf(stderr,"%d %g %g src[%d]=%g\n", ii, Tm, Tp, i, src[i]);
                //}
            }
        }
        esrc += (vdotTp - vdotTm)/ND_ELEM(geom,ii,jj,kk).volume;

        #if(NDIM>1)
        vdotTm = 0.;
        vdotTp = 0.;
        for(int i=0;i<SPACEDIM;i++) {
            Tm = 0.5*(NDP_ELEM(Tvisc,ii,jj,kk,SPACEDIM+i) + NDP_ELEM(Tvisc,ii,jj-1,kk,SPACEDIM+i));
            vdotTm += Tm*NDP_ELEM(sim.vedge[1],ii,jj,kk,i)*ND_ELEM(geom,ii,jj,kk).area[1];
            Tp = 0.5*(NDP_ELEM(Tvisc,ii,jj,kk,SPACEDIM+i) + NDP_ELEM(Tvisc,ii,jj+1,kk,SPACEDIM+i));
            vdotTp += Tp*NDP_ELEM(sim.vedge[1],ii,jj+1,kk,i)*ND_ELEM(geom,ii,jj+1,kk).area[1];
            src[i] -= (Tp*ND_ELEM(geom,ii,jj+1,kk).area[1] - Tm*ND_ELEM(geom,ii,jj,kk).area[1])/ND_ELEM(geom,ii,jj,kk).volume;
        }
        esrc += (vdotTp - vdotTm)/ND_ELEM(geom,ii,jj,kk).volume;
        #endif

        // now add in geometric corrections
        for(int i=0;i<SPACEDIM;i++) {
            for(int j=0;j<SPACEDIM;j++) {
                for(int k=0;k<SPACEDIM;k++) {
                    src[j] += ND_ELEM(geom,ii,jj,kk).conn[i][j][k]*NDP_ELEM(Tvisc,ii,jj,kk,k*SPACEDIM+i);
                    //if(ii > 1020 && j==1) {
                    //    fprintf(stderr,"%d src[%d]=%g   %g %g\n", ii, j, src[j], ND_ELEM(geom,ii,jj,kk).conn[i][j][k], NDP_ELEM(Tvisc,ii,jj,kk,k*SPACEDIM+i));
                    //}
                }
            }
        }

        // cooling from the surface
        #if(NDIM==1)
/*
        double coeff = NDP_ELEM(p,ii,jj,kk,RHO)*KBOLTZ/(2*SBOLTZ*MP*dtpush);
        double Tguess = NDP_ELEM(sim.eos,ii,jj,kk,TEMP);
        double Tlast = Tguess;
        int iter = 0;
        do {
            double f = pow(Tguess,4) + coeff*(Tguess - NDP_ELEM(sim.eos,ii,jj,kk,TEMP));
            double df = 4*pow(Tguess,3) + coeff;
            Tlast = Tguess;
            Tguess -= f/df;
            iter++;
            //fprintf(stderr,"%d %g %g\n", iter, Tguess, NDP_ELEM(sim.eos,ii,jj,kk,TEMP));
        } while(fabs(Tguess-Tlast)/NDP_ELEM(sim.eos,ii,jj,kk,TEMP) > 1.e-4 && iter < 10);
        if(iter == 10) {
            fprintf(stderr,"failed to converge on T for surface cooling\n");
            exit(3);
        }

        esrc += 2*SBOLTZ*pow(Tguess,4);*/
        #endif

        NDP_ELEM(sim.src,ii,jj,kk,ETOT) -= esrc;
        for(int i=0;i<SPACEDIM;i++) NDP_ELEM(sim.src,ii,jj,kk,U1+i) -= src[i];
        //if(istop[0]==n1) {
        //    fprintf(stderr,"%d %g %g\n", ii, NDP_ELEM(Tvisc,ii,jj,kk,1), src[1]);//src[0], src[1]);
        //    fflush(stderr);
        //}
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //exit(2);




}

void aux_src(double NDP_PTR p, double dtpush) {

}
