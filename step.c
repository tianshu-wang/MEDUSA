#include "decs.h"
#include "constants.h"

// #define RAY_BY_RAY

GPU_PRAGMA(omp declare target)
double u[NSIZE];
double uh[NSIZE];
double freq_flux[NSIZE];
double gr_freq_flux[NSIZE];
double divF[NSIZE];
double uaux[NSIZE];
double eosaux[NSIZE];
GPU_PRAGMA(omp end declare target)

void update_cpu(){
  TIMER_START("update_cpu");
  GPU_PRAGMA(omp target update from(
 	    sim_p[:cell_count_all][:nvars],sim_eos[:cell_count_all][:NEOS],geom[:cell_count_all],\
	    sim_shock_flag[:cell_count_all]
	    ))
  GPU_PRAGMA(omp target update from(gr_lapse_edge[:n1+1],gr_lapse[:n1]))
  TIMER_STOP;
}

/* Driver for taking one full time step of all the equations.
 */
void step()
{
  static int firstc = 1;
  int ii,jj,kk,vv,status,debug,j;

#if (USE_LINEAR_ALLOCATION==TRUE)
  void push(const double dtpush, double ** p, double ** ph, double ** pf, int cycle);
  void push_tracers(double ** p, double dtpush, int stage);
#else
  void push(const double dtpush, double NDP_PTR p, double NDP_PTR ph, double NDP_PTR pf, int cycle);
  void push_tracers(double NDP_PTR p, double dtpush, int stage);
#endif
  #if (USE_LINEAR_ALLOCATION==TRUE)
  void ND_FUNC(prim_to_fluxes,NDIM)(double ** p);
  #else
  void ND_FUNC(prim_to_fluxes,NDIM)(double NDP_PTR p);
  #endif

  TIMER_START("step");

  dt_flag = -1;
  min_dt  = 1.0e100;
  vmax    = 0.0;
  taumax  = 0.0;

  debug = 0;
/*
  if (firstc) {
    {
      u            = malloc_rank1(nvars,                    sizeof *u           );
      uh           = malloc_rank1(nvars,                    sizeof *uh          );
      uaux         = malloc_rank1(nvars,                    sizeof *uaux        );
      eosaux       = malloc_rank1(NEOS,                     sizeof *eosaux      );
      divF         = malloc_rank1(nvars,                    sizeof *divF        );
      freq_flux    = malloc_rank1(((ngroups+1)*(1+NDIM)-1), sizeof *freq_flux   );
      #if (GR_MONOPOLE==TRUE)
      gr_freq_flux = malloc_rank1(((ngroups+1)*(1+NDIM)-1), sizeof *gr_freq_flux);
      #endif
    }
    firstc = 0;
  }
*/
  // At this point, assume everything is on gpu and mpi-related cells are updated between cpu and gpu.
  #if (DO_HYDRO==TRUE)
  set_shock_flag(0);  //on gpu // initialize shock_flag so that no grid aligned shocks are tagged
  reset_boundaries_gpu(sim_p,1); // on gpu
  GPU_PRAGMA(omp target update from(sim_p[:cell_count_send][:nvars],sim_p[cell_count:cell_count_recv][:nvars],\
			            sim_shock_flag[:cell_count_send],sim_shock_flag[cell_count:cell_count_recv]))
  sync_mpi_boundaries(sim_p); //cpu, mpi communication
  update_eos(eostable,sim_p,sim_eos); // on gpu, 0:cell_count
  get_geom_src_terms(sim_p); // on gpu, 0:cell_count
  p_to_phys_src_start(sim_p); // calculation on gpu, 0:cell_count, gravity related data map to cpu
  complete_mpi_communication(0); // on cpu, modification on mpi cells
  GPU_PRAGMA(omp target update to(sim_p[cell_count:cell_count_recv][:nvars],\
  			          sim_shock_flag[cell_count:cell_count_recv]))
  complete_mpi_bc(0); // on gpu
  update_eos_ghost(eostable,sim_p,sim_eos); // on gpu
  ND_FUNC(prim_to_fluxes,NDIM)(sim_p); // on gpu
  p_to_phys_src_finish(sim_p,1*dt); // on gpu
  push(dt, sim_p, sim_p, sim_ph, 0); // on gpu
  mpi_recv_tracers(); // ignored
  push_tracers(sim_p, dt, 0); // ignored
  mpi_send_tracers(0); // ignored
  #if (GEOM==SPHERICAL && NDIM==3 && POLAR_AVG==TRUE)
  avg_poles_3d(sim_ph);
  #endif

  // update p^n to p^n+1
  reset_boundaries_gpu(sim_ph,1); 
  GPU_PRAGMA(omp target update from(sim_ph[:cell_count_send][:nvars],sim_ph[cell_count:cell_count_recv][:nvars],\
			            sim_shock_flag[:cell_count_send],sim_shock_flag[cell_count:cell_count_recv]))
  sync_mpi_boundaries(sim_ph);
  update_eos(eostable,sim_ph,sim_eos);
  get_geom_src_terms(sim_ph); 
  p_to_phys_src_start(sim_ph); 
  complete_mpi_communication(1); 
  GPU_PRAGMA(omp target update to(sim_ph[cell_count:cell_count_recv][:nvars],\
			          sim_shock_flag[cell_count:cell_count_recv]))
  complete_mpi_bc(1);
  update_eos_ghost(eostable,sim_ph,sim_eos);
  ND_FUNC(prim_to_fluxes,NDIM)(sim_ph);
  p_to_phys_src_finish(sim_ph,0.5*dt); 
  push(dt, sim_p, sim_ph, sim_p, 1);
  mpi_recv_tracers();
  push_tracers(sim_p, dt, 1);
  mpi_send_tracers(1);
  #if (GEOM==SPHERICAL && NDIM==3 && POLAR_AVG==TRUE)
  avg_poles_3d(sim_p);
  #endif


  #endif

  

  if (detect_tbounce && tbounce<0.0) set_tbounce();


  TIMER_STOP;

  return;
}

/* Compute bounce time */
/* THIS IS VERY SPECIFIC TO THE CORE-COLLAPSE PROBLEM */
void set_tbounce()
{
  int ii,jj,kk;
  double rho_cent = 0.0;

  if (istart[0]==0) {
    ii = istart[0];
    jj = JS(ii,istart[1]);
    kk = KS(ii,jj,istart[2]);
    rho_cent = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
  }

  rho_cent = mpi_max(rho_cent);

  if (rho_cent > rho_cent_max) rho_cent_max = rho_cent;

  if (rho_cent < 0.95*rho_cent_max) {
    tbounce = t;
    if (mpi_io_proc()) fprintf(stderr,"[set_tbounce]:  Bounce detected! tbounce=%1.6f\n",tbounce);
  }

  return;
}

/* Compute/sync an estimate for the next time step.  Note that min_dt is, except
 *   for the first step, already computed based on CFL/accuracy considerations in
 *   the course of updating the variables in the last step.  For the first step,
 *   we just estimate it using the hydro CFL condition and limit it based on the
 *   user-specified initial_dt and dtmax
 */
void estimate_dt()
{
  static int firstc = 1;
  int ii,jj,kk,dd;
  double local_dt,local_dx;

  TIMER_START("estimate_dt");

  if (firstc) {
    min_dt = 1.0e100;
    firstc = 0;
  }

  #if (DO_HYDRO==TRUE)
  ZLOOP {
    local_dt = 0.0;
    DLOOP {
      local_dx = dx[dd];
      if (dd==1) local_dx *= DJS(ii);
      if (dd==2) local_dx *= DKS(ii,jj);
      local_dt += (fabs(NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1+dd))+NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,CS)/ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][dd])/local_dx;

      min_dt = MIN(min_dt,1.0/local_dt);
    }
  }
  #else
  min_dt = dtmax;
  #endif

  if (istep==0) {
    min_dt = MIN(min_dt, initial_dt);
  }
  min_dt = MIN(min_dt, dtmax);

  min_dt *= cfl;
  dt = mpi_min(min_dt);

  TIMER_STOP;

  return;
}

double dt_source_control(double *p, double *u, double *s)
{
  int dd;
  const double alpha = 0.02;
  double local_dt,momterm;
  double mdt = 1.0e100;

  local_dt = alpha*u[0]/(fabs(s[0])+1.0e-16);
  mdt      = MIN(mdt,local_dt);
  local_dt = alpha*u[1]/(fabs(s[1])+1.0e-16);
  mdt      = MIN(mdt,local_dt);

  momterm  = 0.0;
  SLOOP { momterm += p[U1+dd]*s[U1+dd]; }
  local_dt = 100*alpha*p[UU]/(fabs(fabs(s[1])-fabs(momterm)) + 1.0e-16);
  mdt      = MIN(mdt,local_dt);

  return mdt;
}

void interp_tracer_velocities(double x[], int index[], double v[])
{
    double gradv[SPACEDIM][SPACEDIM], vavg[SPACEDIM], xavg[NDIM];
    calc_dvdx(index[0], index[1], index[2], gradv, vavg, xavg);
    for(int d=0;d<NDIM;d++) v[d] = vavg[d];
    for(int d=0;d<NDIM;d++) {
        v[d] = vavg[d];
        for(int dd=0;dd<NDIM;dd++) {
            v[d] += (x[dd]-xavg[dd])*gradv[d][dd];
        }
    }
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void push_tracers(double ** p, double dtpush, int stage)
#else
void push_tracers(double NDP_PTR p, double dtpush, int stage)
#endif
{
    //return;
    int index[] = {0, 0, 0};
    int flag1, flag2, flag3;
    double vel[NDIM];
    tracer_t *tr = tracers;
    while(tr != NULL) {
        // first figure out which cell the tracer lives in
        if(!tr->active) {
            tr = tr->next;
            continue;
        }
        if(stage==0) {
            x_to_ijk(tr->x, index);
        } else {
            x_to_ijk(tr->xh, index);
        }
        int outside = is_outside_rank(index, myrank);
        flag1 = flag2 = flag3 = 0;
        if(outside) {
            if(outside == 1) {
                if(index[0] == istop[0]) index[0]--;
                if(index[0] == istart[0]-1) index[0]++;
                outside = is_outside_rank(index,myrank);
                flag1 = 1;
            }
            if(outside == 2) {
                if(index[1] == JS(index[0],istop[1])) index[1]--;
                if(index[1] == JS(index[0],istart[1])-1) index[1]++;
                outside = is_outside_rank(index,myrank);
                flag2 = 1;
            }
            if(outside == 3) {
                if(index[2] == KS(index[0],index[1],istop[2])) index[2]--;
                if(index[2] == KS(index[0],index[1],istart[2])-1) index[2]++;
                outside = is_outside_rank(index,myrank);
                flag3 = 1;
            }
            if(outside) {
                fprintf(stderr,"out of bounds in stage %d on proc %d!!!   %d %d %d  %d   %d\n", stage, myrank, istart[outside-1], index[outside-1], istop[outside-1], tr->id, DJS(index[0]));
                #if (NDIM>1)
                fprintf(stderr,"outside = %d   %g %g\n", outside, tr->x[1], tr->xh[1]);
                #endif
                fprintf(stderr,"flags = %d  %d  %d\n", flag1, flag2, flag3);
                tr->active = 0;
                tr = tr->next;
                continue;
            }
        }
        if(stage == 0) {
            interp_tracer_velocities(tr->x, index, vel);
        } else {
            interp_tracer_velocities(tr->xh, index, vel);
        }

        if(stage == 0) {
            for(int j = 0; j < NDIM; j++) {
                tr->xh[j] = tr->x[j] + dtpush*vel[j];
            }
        } else {
            for(int j = 0; j < NDIM; j++) {
                tr->x[j] = 0.5*(tr->x[j] + tr->xh[j] + dtpush*vel[j]);
            }
        }
        // deal with axis
        #if(NDIM>1)
        if(stage==1) {
            if(tr->x[1] < startx[1]) tr->x[1] = startx[1] + (startx[1] - tr->x[1]);
            if(tr->x[1] > startx[1]+n2*dx[1]) tr->x[1] = startx[1]+n2*dx[1] - (tr->x[1] - startx[1]-n2*dx[1]);
        } else {
            if(tr->xh[1] < startx[1]) tr->xh[1] = startx[1] + (startx[1] - tr->xh[1]);
            if(tr->xh[1] > startx[1]+n2*dx[1]) tr->xh[1] = startx[1]+n2*dx[1] - (tr->xh[1] - startx[1]-n2*dx[1]);
        }

        #endif
        tr = tr->next;
    }
}

/* A single stage of Runge-Kutta
 * Pushes primitives p to pf using ph to calculate fluxes/source terms/etc
 * This is where all the physics is advanced.  There are no external operator split updates
 */
#if (USE_LINEAR_ALLOCATION==TRUE)
void push(const double dtpush, double ** p, double ** ph, double ** pf, int stage)
#else
void push(const double dtpush, double NDP_PTR p, double NDP_PTR ph, double NDP_PTR pf, int stage)
#endif
{
  int ii=0,jj=0,kk=0,vv,dd,status,g;
  int j,jp,njp,k,kp,nkp;
  zone_geom *gm;
  double gradv[SPACEDIM][SPACEDIM];
  double pfact,phfact,ffact,local_dt;
  double min_dt_src;
  double ddt;
  double min_dt_inel,Etot_old,Ye_old,T_old;
  //#if (USE_LINEAR_ALLOCATION==TRUE)
  //void ND_FUNC(prim_to_fluxes,NDIM)(double ** p);
  //#else
  //void ND_FUNC(prim_to_fluxes,NDIM)(double NDP_PTR p);
  //#endif

  TIMER_START("push");

  if (stage == 0) {
    pfact  = 1;
    phfact = 0;
    ffact  = 1;
  } else if (stage == 1) {
    pfact  = 0.5;
    phfact = 0.5;
    ffact  = 0.5;
  } else {
    pfact  = 1.0/3.0;
    phfact = 2.0/3.0;
    ffact  = phfact;
  }


  //move outside push()
  //// compute fluxes
  //ND_FUNC(prim_to_fluxes,NDIM)(ph);

  //// compute physical source terms
  //p_to_phys_src_finish(ph,ffact*dtpush); 

  min_dt_src = 1.0e100;


  GPU_PRAGMA(omp target data map(tofrom:gr_grav[:n1])){
  GPU_PRAGMA(omp target teams distribute parallel for \
    firstprivate(ffact,pfact,phfact,nhydro,rho_floor,e_floor) \
    private(ii,jj,kk,dd,vv,gm,gradv,local_dt,\
            ddt,jp,njp,kp,nkp,j,k,g,status,u,uh,uaux,eosaux,freq_flux,gr_freq_flux,divF)\
    shared(p,ph,pf) reduction(min: min_dt_src) \
  )
  for(int II=0;II<cell_count;++II){
  //ZLOOP_LINEAR {
    //#if (USE_LINEAR_ALLOCATION==TRUE)
    //get_ijk_from_I(II,&ii,&jj,&kk);
    //#endif
    GET_IJK_FROM_I(II,ii,jj,kk)
    // make a local reference to the zone's geometry info
    gm = &ND_ELEM_LINEAR_REORDERED(geom,ii,jj,kk);
    // get the conserved variables
    p_to_u(ND_ELEM_LINEAR_REORDERED(p,ii,jj,kk),u,gm,nhydro);
    p_to_u(ND_ELEM_LINEAR_REORDERED(ph,ii,jj,kk),uh,gm,nhydro);

    // estimate next time step limits from source terms
    local_dt = dt_source_control(ND_ELEM_LINEAR_REORDERED(p,ii,jj,kk), u, ND_ELEM_LINEAR_REORDERED(sim_src,ii,jj,kk));
    min_dt_src = MIN(min_dt_src,local_dt);

    VLOOP { divF[vv] = 0.0; }

    /* Flux differences in 0-direction */
    jp = jj*DJS(ii)/DJS(ii+1);
    njp = (DJS(ii+1) < DJS(ii)) ? DJS(ii)/DJS(ii+1) : 1;
    for (j=jp; j<jp+njp; j++) {
      nkp = DKS(ii,jj)/DKS(ii+1,j);
      kp = kk*nkp;
      for (k=kp; k<kp+nkp; k++) {
        VLOOP { divF[vv] += NDP_ELEM_LINEAR_F(sim_fdir0,0,ii+1,j,k,vv); }
      }
    }
    VLOOP { divF[vv] -= NDP_ELEM_LINEAR_F(sim_fdir0,0,ii,jj,kk,vv); }

    #if (NDIM>1)
    /* Flux differences in 1-direction */
    /* Outer 1-face */
    if (DKS(ii,jj+1) > DKS(ii,jj)) {
      // Coarsening boundary, 1 outer neighbor, shifted area index, half area, original flux index
      VLOOP { divF[vv] += 0.5*NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,kk,vv); }
    } else if (DKS(ii,jj+1) < DKS(ii,jj)) {
      // Refinement boundary, 2 outer neighbors, shifted area index, original area, shifted flux index
      VLOOP { divF[vv] += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,2*kk,vv) + NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,2*kk+1,vv); }
    } else {
      // Regular boundary, 1 outer neighbor, original area index, original area, original flux index
      VLOOP { divF[vv] += NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj+1,kk,vv); }
    }

    /* Inner 1-face */
    if (DKS(ii,jj-1) < DKS(ii,jj)) {
      // Coarsening boundary, 2 inner neighbors, original area index, half area, shifted flux index
      VLOOP { divF[vv] -= 0.5*NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,2*kk,vv) + 0.5*NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,2*kk+1,vv); }
    } else {
      // Refinement boundary, 1 inner neighbor, original area index, original area, original flux index
      // Regular boundary, 1 inner neighbor, original area index, original area, original flux index
      VLOOP { divF[vv] -= NDP_ELEM_LINEAR_F(sim_fdir1,1,ii,jj,kk,vv); }
    }
    #endif

    #if (NDIM>2)
    /* Flux differences in 2-direction */
    VLOOP { divF[vv] += NDP_ELEM_LINEAR_F(sim_fdir2,2,ii,jj,kk+1,vv) - NDP_ELEM_LINEAR_F(sim_fdir2,2,ii,jj,kk,vv); }
    #endif


    VLOOP {
      ddt  = divF[vv]/gm->volume - NDP_ELEM_LINEAR_REORDERED(sim_src,ii,jj,kk,vv);
      u[vv] = pfact*u[vv] + phfact*uh[vv] - ffact*dtpush * ddt;
    }

    // operator split stiff source terms for energy equation
    Etot_old = 0;
    GLOOP{Etot_old += u[irad1+g];}
    T_old = NDP_ELEM_LINEAR_REORDERED(sim_eos,ii,jj,kk,TEMP);
    //GLOOP{if(u[irad1+g]<-1e10)fprintf(stderr,"test %d %g %g\n",g,u[irad1+g],u[irad1+g]*spec_factor[g]);}

#ifdef RAY_BY_RAY
GLOOP {
  DLOOP {
    if (dd>0) {
      u [ifrad1+g*NDIM+dd] = 0.0;
    }
  }
}
#endif

    //if (status < 1) {
    //  fprintf(stderr,"implicit update failed with error %d in zone (%d,%d,%d) on step %d\n", status, ii,jj,kk, istep);
    //  exit(43);
    //}

    // transform back to primitive variables
    u_to_p(u,ND_ELEM_LINEAR_REORDERED(pf,ii,jj,kk),gm,nhydro,rho_floor,e_floor);

  } // AS:  end of ZLOOP
}//end of omp target data map

  min_dt = MIN(min_dt, min_dt_src);

  TIMER_STOP;

  return;
}


#if (NDIM==1)
void prim_to_fluxes1(double NDP_PTR p)
{
  int i,vv,g,dd,isshock;
  double sig_speed,l_zone,last_vol, vriemann;
  double pl[NSIZE], pr[NSIZE];
  double alpha[2*NG],beta[2*NG],Gamma[2*NG],x[2*NG];
  int firstc = 1;
  double kappa[ngroups];
  double sc[ngroups];
  double delta[ngroups];
  double chil[ngroups];
  double chir[ngroups];
  double flux[NSIZE];
  double vedge[NSIZE];

  TIMER_START("prim_to_fluxes1");

  GPU_PRAGMA(omp target data map(alloc:pl,pr,alpha,beta,Gamma,x,kappa,sc,delta,chil,chir,flux,vedge)){
	GPU_PRAGMA(omp target teams distribute parallel for private(alpha,beta,Gamma,x,sig_speed,pl,pr,chir,chil,flux,vedge) shared(sim_fdir0,sim_vedgedir0) reduction(min:min_dt))
	for (i=istart[0]; i<=istop[0]; i++) {
                double pencil[NSIZE][2*NG];
                double pleft [NSIZE][2*NG];
                double pright[NSIZE][2*NG];
	        for (int itemp=i-NG; itemp<i+NG; itemp++) {
	        	alpha[NG+itemp-i] = alpha0[itemp];
	        	beta [NG+itemp-i] = beta0 [itemp];
	        	Gamma[NG+itemp-i] = Gamma0[itemp];
	        	x    [NG+itemp-i] = startx[0] + itemp*dx[0];
	        	VLOOP {
	        		pencil[vv][NG+itemp-i] = NDP_ELEM_LINEAR(p,itemp,jj,kk,vv);
	        	}
	        	pencil[vv][NG+itemp-i] = NDP_ELEM_LINEAR(sim_eos,itemp,jj,kk,PRESS);
	        	vv++;
	        	pencil[vv][NG+itemp-i] = NDP_ELEM_LINEAR(sim_eos,itemp,jj,kk,GAMMA);
	        }

	        // reconstruct variables
	        //interp(pencil, pleft, pright, alpha, beta, Gamma, x, my_grid_dims[0]);
	        int size = 0;
	        ILOOP {
	        	switch(interp_order[vv]) {
	        		case 3:
	        			para_recon(pencil[vv], pleft[vv], pright[vv], alpha, beta, Gamma, x, -1, size);
	        			break;
	        		case 2:
	        			lin_recon(pencil[vv], pleft[vv], pright[vv], beta, Gamma, x, -1, size);
	        			break;
	        		default:
	        			for (int itemp=-1; itemp<=size; itemp++) {
	        				pleft[vv][itemp+NG] = pright[vv][itemp+NG] = pencil[vv][itemp+NG];
	        			}
	        			break;
	        	}
	        }
		// copy out reconstructed left/right states
		// note that left state is really right state of cell i-1
		// and right state is really left state of cell i
		VLOOP {
			pl[vv] = pright[vv][NG-1];
			pr[vv] = pleft [vv]  [NG];
		}
		pl[vv] = pright[vv][NG-1];
		pr[vv] = pleft [vv]  [NG];
		vv++;
		pl[vv] = pright[vv][NG-1];
		pr[vv] = pleft [vv]  [NG];

		#if (DO_HYDRO==TRUE)
		// get fluxes of hydro variables at this interface
                double flux[NSIZE],vriemann[NSIZE];
		//isshock = riemann_flux_LLF(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 0, flux, &sig_speed, vriemann);
		isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 0, flux, &sig_speed, vriemann);
		HLOOP{
			NDP_ELEM_LINEAR_F(sim_fdir0,0,i,j,k,vv) = flux[vv];}
		SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k,dd) = vriemann[dd];}
		//isshock = riemann_flux_LLF(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 0, ND_ELEM_LINEAR_F(sim_fdir0,0,i,j,k), &ss, ND_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k));

		// update estimate for next time step, if necessary
		min_dt = MIN(min_dt,dx[0]/fabs(sig_speed));
		#endif
	}
  }
    TIMER_STOP;

	return;
}
#endif


#if (NDIM==2)
/* Perform a prolongation in the x0-direction using linear interpolation */
/* (ip,jp) is cell where the flux is to be computed, and ii is the current
 * radial level being prolongated.  Thus, the prolongation point will be at
 * (x0p,x1p), where x0p=x0(ii) and x1p=x1(jp)
 */

#if (USE_LINEAR_ALLOCATION==TRUE)
void prolongate0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,int *dj,double **alpha1s,double **beta1s,double **Gamma1s,double ** p, double ** e, int ip, int jp, int ii,
  double* pp)
#else
void prolongate0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,int *dj,double **alpha1s,double **beta1s,double **Gamma1s,double NDP_PTR p, double NDP_PTR e, int ip, int jp, int ii,
  double* pp)
#endif
{
  int jj,jlev,s,vv;
  double x1p,x10,x11,x12,idx1m,idx1p,idx1c,sm,sp,sc,sg;

  // Get local indices on level ii
  jj = jp*DJS(ip)/DJS(ii);

  // Spatial coordinates at the prolongation point
  jlev = 0;
  s = DJS(ip);
  while (s >>= 1) jlev++;
  x1p = beta1s[jlev][jp]/Gamma1s[jlev][jp];

  // Spatial coordinates of coarse zones
  jlev = 0;
  s = DJS(ii);
  while (s >>= 1) jlev++;
  x10 = beta1s[jlev][jj-1]/Gamma1s[jlev][jj-1];
  x11 = beta1s[jlev][jj  ]/Gamma1s[jlev][jj  ];
  x12 = beta1s[jlev][jj+1]/Gamma1s[jlev][jj+1];
  idx1m = 1.0/(x11-x10);
  idx1p = 1.0/(x12-x11);
  idx1c = 1.0/(x12-x10);

  // Now compute limited slopes in the 1-direction
  VLOOP {
    sm = (NDP_ELEM_LINEAR(p,ii,jj  ,0,vv) - NDP_ELEM_LINEAR(p,ii,jj-1,0,vv))*idx1m;
    sp = (NDP_ELEM_LINEAR(p,ii,jj+1,0,vv) - NDP_ELEM_LINEAR(p,ii,jj  ,0,vv))*idx1p;
    sc = (NDP_ELEM_LINEAR(p,ii,jj+1,0,vv) - NDP_ELEM_LINEAR(p,ii,jj-1,0,vv))*idx1c;
    sg = MINMOD(2.0*MINMOD(sm,sp),sc);
    pp[vv] = NDP_ELEM_LINEAR(p,ii,jj,0,vv) + sg*(x1p-x11);
  }
  sm = (NDP_ELEM_LINEAR(e,ii,jj  ,0,PRESS) - NDP_ELEM_LINEAR(e,ii,jj-1,0,PRESS))*idx1m;
  sp = (NDP_ELEM_LINEAR(e,ii,jj+1,0,PRESS) - NDP_ELEM_LINEAR(e,ii,jj  ,0,PRESS))*idx1p;
  sc = (NDP_ELEM_LINEAR(e,ii,jj+1,0,PRESS) - NDP_ELEM_LINEAR(e,ii,jj-1,0,PRESS))*idx1c;
  sg = MINMOD(2.0*MINMOD(sm,sp),sc);
  pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,0,PRESS) + sg*(x1p-x11);
  vv++;
  sm = (NDP_ELEM_LINEAR(e,ii,jj  ,0,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj-1,0,GAMMA))*idx1m;
  sp = (NDP_ELEM_LINEAR(e,ii,jj+1,0,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj  ,0,GAMMA))*idx1p;
  sc = (NDP_ELEM_LINEAR(e,ii,jj+1,0,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj-1,0,GAMMA))*idx1c;
  sg = MINMOD(2.0*MINMOD(sm,sp),sc);
  pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,0,GAMMA) + sg*(x1p-x11);
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void restrict0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,int *dj,double **alpha1s,double **beta1s,double **Gamma1s,double ** p, double ** e, zone_geom * g, int i, int j, int ii,
  double* pp)
#else
void restrict0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,int *dj,double **alpha1s,double **beta1s,double **Gamma1s,double NDP_PTR p, double NDP_PTR e, zone_geom ND_PTR g, int i, int j, int ii,
  double* pp)
#endif
{
  int njp,jp,jj,vv;
  double vol_sum,vol_frac;

  njp = DJS(i)/DJS(ii);
  jp = j*njp;
  vol_sum = 0.0;
  for (jj=jp; jj<jp+njp; jj++) {
    vol_sum += ND_ELEM_LINEAR(g,ii,jj,0).volume;
  }
  vol_sum = 1.0/vol_sum;

  ILOOP { pp[vv] = 0.0; }
  for (jj=jp; jj<jp+njp; jj++) {
    vol_frac = ND_ELEM_LINEAR(g,ii,jj,0).volume*vol_sum;
    VLOOP { pp[vv] += NDP_ELEM_LINEAR(p,ii,jj,0,vv)*vol_frac; }
    pp[vv++] += NDP_ELEM_LINEAR(e,ii,jj,0,PRESS)*vol_frac;
    pp[vv  ] += NDP_ELEM_LINEAR(p,ii,jj,0,GAMMA)*vol_frac;
  }
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void prim_to_fluxes2(double ** p)
#else
void prim_to_fluxes2(double NDP_PTR p)
#endif
{
  static int firstc = 1;
  int i,i0,ii,j,jj,k,kk,dd,vv,g;
  //double kappa[ngroups]; //need statically allocated later
  //double sc[ngroups];
  //double delta[ngroups];
  //double chil[ngroups];
  //double chir[ngroups];

  TIMER_START("prim_to_fluxes2");

  //GPU_PRAGMA(omp target ){
	GPU_PRAGMA(omp target teams distribute parallel for shared(sim_fdir0,sim_vedgedir0,sim_fdir1,sim_vedgedir1) reduction(min:min_dt))
  for(int II=0;II<cell_count_all;II++) {
    double pencil[NSIZE][2*NG];
    double pleft [NSIZE][2*NG];
    double pright[NSIZE][2*NG];
    double flux[NSIZE],vriemann[NSIZE];
    double pl[NSIZE], pr[NSIZE],pp[NSIZE];
    double alpha[2*NG],beta[2*NG],Gamma[2*NG],x[2*NG];
    int jstart,jstop;
    double ss,vol_sum,vol_frac;
    int jm,jp,njp,isshock,jlev,s;
    GET_IJK_FROM_I(II,i,j,k);
    if(i<istart[0]||i>istop[0]||j<JS(i,istart[1])||j>=JS(i,istop[1])){/*do nothing*/;}
    else{
      for (ii=i-NG; ii<i+NG; ii++) {
        alpha[NG+ii-i] = alpha0[ii];
        beta [NG+ii-i] = beta0 [ii];
        Gamma[NG+ii-i] = Gamma0[ii];
        x    [NG+ii-i] = startx[0] + ii*dx[0];

        if (DJS(ii) < DJS(i)) {
          restrict0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,p,sim_eos,geom,i,j,ii,pp);
        } else {
          prolongate0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,p,sim_eos,i,j,ii,pp);
        }
        ILOOP { pencil[vv][NG+ii-i] = pp[vv]; }
      }
      interp(ninterp,interp_order,pencil, pleft, pright, alpha, beta, Gamma, x, 0);
      ILOOP {
        pl[vv] = pright[vv][NG-1];
        pr[vv] = pleft [vv][NG+0];
      }
      isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 0, flux, &ss, vriemann);
      HLOOP{NDP_ELEM_LINEAR_F(sim_fdir0,0,i,j,k,vv) = flux[vv];}
      SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k,dd) = vriemann[dd];}

      if (isshock) {
        ND_ELEM_LINEAR(sim_shock_flag,i,j,0) = 0;
        jm = j*DJS(i)/DJS(i-1);
        ND_ELEM_LINEAR(sim_shock_flag,i-1,jm,0) = 0;
        jm = j*DJS(i)/DJS(i-2);
        ND_ELEM_LINEAR(sim_shock_flag,i-2,jm,0) = 0;
      }
      min_dt = MIN(min_dt,dx[0]/fabs(ss));
    }
  //} // end for II

  //// fluxes in 1-direction
  //for(int II=0;II<cell_count_all;II++) {
  //  GET_IJK_FROM_I(II,i,j,k);
    jstart = JS(i,istart[1]);
    jstop  = JS(i,istop[1]);
    if(i<istart[0]||i>=istop[0]||j<jstart||j>jstop){/*do nothing*/;}
    else{
      jlev = 0;
      s = DJS(i);
      while (s >>= 1) jlev++;
      for (int jtemp=j-NG; jtemp<j+NG; jtemp++) {
        alpha[NG+jtemp-j] = alpha1s[jlev][jtemp];
        beta [NG+jtemp-j] = beta1s [jlev][jtemp];
        Gamma[NG+jtemp-j] = Gamma1s[jlev][jtemp];
        x    [NG+jtemp-j] = startx[1] + jtemp*DJS(i)*dx[1];
        VLOOP { pencil[vv][NG+jtemp-j] = NDP_ELEM_LINEAR(p,i,jtemp,kk,vv); }
        pencil[vv++][NG+jtemp-j] = NDP_ELEM_LINEAR(sim_eos,i,jtemp,kk,PRESS);
        pencil[vv  ][NG+jtemp-j] = NDP_ELEM_LINEAR(sim_eos,i,jtemp,kk,GAMMA);
      }
      interp(ninterp,interp_order,pencil, pleft, pright, alpha, beta, Gamma, x, 0);  // interpolate (jstop-jstart) 4-zone pencils

      ILOOP {
        pl[vv] = pright[vv][NG-1];
        pr[vv] = pleft [vv][NG+0];
      }
      if (transverse_shock(ijk_to_I,sim_shock_flag,dj,dk,i,j,0,1)) {
        riemann_flux_HLLE(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 1, flux, &ss, vriemann);
        HLOOP{NDP_ELEM_LINEAR_F(sim_fdir1,1,i,j,k,vv) = flux[vv];}
        SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k,dd) = vriemann[dd];}
      } else {
        isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 1, flux, &ss, vriemann);
        HLOOP{NDP_ELEM_LINEAR_F(sim_fdir1,1,i,j,k,vv) = flux[vv];}
        SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k,dd) = vriemann[dd];}
      }

      // update estimate for next time step, if necessary
      min_dt = MIN(min_dt,DJS(i)*dx[1]/fabs(ss));
    }
  }
  //}
  TIMER_STOP;
  return;
}
#endif

#if (NDIM==3)
//static double* pl    = NULL;
//static double* pr    = NULL;
//static double* pp    = NULL;
//static double* alpha = NULL;
//static double* beta  = NULL;
//static double* Gamma = NULL;
//static double* x     = NULL;
//static double* kappa = NULL;
//static double* sc    = NULL;
//static double* delta = NULL;
//static double* chil  = NULL;
//static double* chir  = NULL;


#if (USE_LINEAR_ALLOCATION==TRUE)
void kslopes0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
              int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
              double ** p, double ** e, int ip, int jp, int kp, int ii, int jj, double *pp)
#else
void kslopes0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
              int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
              double NDP_PTR p, double NDP_PTR e, int ip, int jp, int kp, int ii, int jj, double *pp)
#endif
{
  int kk,vv,jlev,klev,s;
  double sm,sp,sc,sg,x2p,x20,x21,x22,idx2m,idx2p,idx2c;

  // 2-coordinate of prolongation point
  klev = 0;
  s = DKS(ip,jp);
  while (s >>= 1) klev++;
  x2p = beta2s[klev][kp]/Gamma2s[klev][kp];

  // 2-index of coarse zone at ii,jj
  kk = kp*DKS(ip,jp)/DKS(ii,jj);

  // 2-coordinates of coarse zones at ii,jj and kk-1,kk,kk+1
  klev = 0;
  s = DKS(ii,jj);
  while (s >>= 1) klev++;
  x20 = beta2s[klev][kk-1]/Gamma2s[klev][kk-1];
  x21 = beta2s[klev][kk  ]/Gamma2s[klev][kk  ];
  x22 = beta2s[klev][kk+1]/Gamma2s[klev][kk+1];
  idx2m = 1.0/(x21-x20);
  idx2p = 1.0/(x22-x21);
  idx2c = 1.0/(x22-x20);

  // Compute limited k-slopes at the prolongation point ii,jj,kp
  VLOOP {
    sm = (NDP_ELEM_LINEAR(p,ii,jj,kk  ,vv) - NDP_ELEM_LINEAR(p,ii,jj,kk-1,vv))*idx2m;
    sp = (NDP_ELEM_LINEAR(p,ii,jj,kk+1,vv) - NDP_ELEM_LINEAR(p,ii,jj,kk  ,vv))*idx2p;
    sc = (NDP_ELEM_LINEAR(p,ii,jj,kk+1,vv) - NDP_ELEM_LINEAR(p,ii,jj,kk-1,vv))*idx2c;
    sg = MINMOD(2.0*MINMOD(sm,sp),sc);
    pp[vv] = NDP_ELEM_LINEAR(p,ii,jj,kk,vv) + sg*(x2p-x21);
  }
  sm = (NDP_ELEM_LINEAR(e,ii,jj,kk  ,PRESS) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,PRESS))*idx2m;
  sp = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,PRESS) - NDP_ELEM_LINEAR(e,ii,jj,kk  ,PRESS))*idx2p;
  sc = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,PRESS) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,PRESS))*idx2c;
  sg = MINMOD(2.0*MINMOD(sm,sp),sc);
  pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,kk,PRESS) + sg*(x2p-x21);
  vv++;
  sm = (NDP_ELEM_LINEAR(e,ii,jj,kk  ,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,GAMMA))*idx2m;
  sp = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj,kk  ,GAMMA))*idx2p;
  sc = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,GAMMA))*idx2c;
  sg = MINMOD(2.0*MINMOD(sm,sp),sc);
  pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,kk,GAMMA) + sg*(x2p-x21);

  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void kslopes1(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double ** p, double ** e, int ip, int jp, int kp, int jj, double *pp)
#else
void kslopes1(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double NDP_PTR p, double NDP_PTR e, int ip, int jp, int kp, int jj, double *pp)
#endif
{
  int ii,kk,vv,jlev,klev,s;
  double sm,sp,sc,sg,x2p,x20,x21,x22,idx2m,idx2p,idx2c;

  // 2-coordinate of prolongation point
  klev = 0;
  s = MIN(DKS(ip,jp),DKS(ip,jp-1));
  while (s >>= 1) klev++;
  x2p = beta2s[klev][kp]/Gamma2s[klev][kp];

  // 2-index of coarse zone at level ii,jj
  ii = ip;
  kk = kp*MIN(DKS(ip,jp),DKS(ip,jp-1))/DKS(ii,jj);

  // 2-coordinates of coarse zones at ii,jj and kk-1,kk,kk+1
  klev = 0;
  s = DKS(ii,jj);
  while (s >>= 1) klev++;
  x20 = beta2s[klev][kk-1]/Gamma2s[klev][kk-1];
  x21 = beta2s[klev][kk  ]/Gamma2s[klev][kk  ];
  x22 = beta2s[klev][kk+1]/Gamma2s[klev][kk+1];
  idx2m = 1.0/(x21-x20);
  idx2p = 1.0/(x22-x21);
  idx2c = 1.0/(x22-x20);

  VLOOP {
    sm = (NDP_ELEM_LINEAR(p,ii,jj,kk  ,vv) - NDP_ELEM_LINEAR(p,ii,jj,kk-1,vv))*idx2m;
    sp = (NDP_ELEM_LINEAR(p,ii,jj,kk+1,vv) - NDP_ELEM_LINEAR(p,ii,jj,kk  ,vv))*idx2p;
    sc = (NDP_ELEM_LINEAR(p,ii,jj,kk+1,vv) - NDP_ELEM_LINEAR(p,ii,jj,kk-1,vv))*idx2c;
    sg = MINMOD(2.0*MINMOD(sm,sp),sc);
    pp[vv] = NDP_ELEM_LINEAR(p,ii,jj,kk,vv) + sg*(x2p-x21);
  }
  sm = (NDP_ELEM_LINEAR(e,ii,jj,kk  ,PRESS) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,PRESS))*idx2m;
  sp = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,PRESS) - NDP_ELEM_LINEAR(e,ii,jj,kk  ,PRESS))*idx2p;
  sc = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,PRESS) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,PRESS))*idx2c;
  sg = MINMOD(2.0*MINMOD(sm,sp),sc);
  pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,kk,PRESS) + sg*(x2p-x21);
  vv++;
  sm = (NDP_ELEM_LINEAR(e,ii,jj,kk  ,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,GAMMA))*idx2m;
  sp = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj,kk  ,GAMMA))*idx2p;
  sc = (NDP_ELEM_LINEAR(e,ii,jj,kk+1,GAMMA) - NDP_ELEM_LINEAR(e,ii,jj,kk-1,GAMMA))*idx2c;
  sg = MINMOD(2.0*MINMOD(sm,sp),sc);
  pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,kk,GAMMA) + sg*(x2p-x21);

  return;
}


/* Perform a prolongation in the x0-direction using bilinear interpolation */
/* (ip,jp,kp) is cell where the flux is to be computed, and ii is the current
 * radial level being prolongated.  Thus, the prolongation point will be at
 * (x0p,x1p,x2p), where x0p=x0(ii), x1p=x1(jp), and x2p=x2(kp)
 */
#if (USE_LINEAR_ALLOCATION==TRUE)
void prolongate0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double ** p, double ** e, int ip, int jp, int kp, int ii,
  double* pp)
#else
void prolongate0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double NDP_PTR p, double NDP_PTR e, int ip, int jp, int kp, int ii,
  double* pp)
#endif
{
  int jj,kk0,kk1,kk2,jlev,s,vv;
  double x1p,x10,x11,x12,idx1m,idx1p,idx1c,sm,sp,sc,sg;
  double pp0[NSIZE],pp1[NSIZE],pp2[NSIZE];

  // Get local indices on level ii
  jj = jp*DJS(ip)/DJS(ii);

  // Spatial coordinates at the prolongation point
  jlev = 0;
  s = DJS(ip);
  while (s >>= 1) jlev++;
  x1p = beta1s[jlev][jp]/Gamma1s[jlev][jp];

  // Compute k-slopes at jj-1, jj, and jj+1
  kslopes0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,e,ip,jp,kp,ii,jj-1,pp0);
  kslopes0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,e,ip,jp,kp,ii,jj  ,pp1);
  kslopes0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,e,ip,jp,kp,ii,jj+1,pp2);

  // Spatial coordinates of coarse zones
  jlev = 0;
  s = DJS(ii);
  while (s >>= 1) jlev++;
  x10 = beta1s[jlev][jj-1]/Gamma1s[jlev][jj-1];
  x11 = beta1s[jlev][jj  ]/Gamma1s[jlev][jj  ];
  x12 = beta1s[jlev][jj+1]/Gamma1s[jlev][jj+1];
  idx1m = 1.0/(x11-x10);
  idx1p = 1.0/(x12-x11);
  idx1c = 1.0/(x12-x10);

  // Now compute limited slopes in the 1-direction
  ILOOP {
    sm = (pp1[vv] - pp0[vv])*idx1m;
    sp = (pp2[vv] - pp1[vv])*idx1p;
    sc = (pp2[vv] - pp0[vv])*idx1c;
    sg = MINMOD(2.0*MINMOD(sm,sp),sc);
    pp[vv] = pp1[vv] + sg*(x1p-x11);
  }
}

/* Perform a prolongation in the x1-direction using bilinear interpolation */
/* COARSENING CASE:
 *    +-------------------------------+-------------------------------+-------------------------------+
 *    |                               |                               |                               |
 *    |                               |                               |                               |
 *    |               *               |       o       *       o       |               *               | <--both jj and j
 *    |            (jj,kk-1)          |            (jj,kk)            |            (jj,kk+1)          |
 *  ^ |              x20              |              x21              |              x22              |
 *  | +---------------+---------------+-L-L-L-L-L-L-L-+-R-R-R-R-R-R-R-+---------------+---------------+
 *  j |               |               |  f1(jp-1,kp)  |  f1(jp-1,kp)  |               |               |
 *    |               |               |               |               |               |               |
 *    |               |               |       .       |       .       |               |               |
 *    |               |               |  (jp-1,kp)    |  (jp-1,kp)    |               |               |
 *    |               |               |      x2p      |      x2p      |               |               |
 *    +---------------+---------------+---------------+---------------+---------------+---------------+
 *                                                   k-->
 *
 *  Note that ip is the same for all values used in the prolongation.
 *  Input:  indices ip,jp,kp and data points (ii,jj,kk-1), (ii,jj,kk), and (ii,jj,kk+1)
 *  marked with '*', where ii and ip are the same
 *  Output:  Linear interpolation of profile at prolongation points marked with 'o'
 *  Since we need both fluxes marked 'L' and 'R', this must be called once per side
 */
#if (USE_LINEAR_ALLOCATION==TRUE)
void prolongate1(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double ** p, double ** e, int ip, int jp, int kp, int jj,
  double *pp)
#else
void prolongate1(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double NDP_PTR p, double NDP_PTR e, int ip, int jp, int kp, int jj,
  double *pp)
#endif
{
  int ii,kk,vv;

  // Get local 2-index on level jj
  ii = ip;
  kk = kp*MIN(DKS(ip,jp),DKS(ip,jp-1))/DKS(ii,jj);

  if (MIN(DKS(ip,jp),DKS(ip,jp-1))==DKS(ii,jj)) {
    // No prolongation to be done!
    VLOOP {  pp[vv] = NDP_ELEM_LINEAR(p,ii,jj,kk,vv); }
    pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,kk,PRESS);
    vv++;
    pp[vv] = NDP_ELEM_LINEAR(e,ii,jj,kk,GAMMA);
  }
  else {
    kslopes1(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,e,ip,jp,kp,jj,pp);
  }
}
#if (USE_LINEAR_ALLOCATION==TRUE)
void restrict0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double ** p, double ** e, zone_geom * g, int i, int j, int k, int ii,
  double* pp)
#else
void restrict0(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double NDP_PTR p, double NDP_PTR e, zone_geom ND_PTR g, int i, int j, int k, int ii,
  double* pp)
#endif
{
  int njp,jp,jj,nkp,kp,kk,vv;
  double vol_sum,vol_frac;

  njp = DJS(i)/DJS(ii);
  jp = j*njp;
  vol_sum = 0.0;
  for (jj=jp; jj<jp+njp; jj++) {
    nkp = DKS(i,j)/DKS(ii,jj);
    kp = k*nkp;
    for (kk=kp; kk<kp+nkp; kk++) {
      vol_sum += ND_ELEM_LINEAR(g,ii,jj,kk).volume;
    }
  }
  vol_sum = 1.0/vol_sum;

  ILOOP { pp[vv] = 0.0; }
  for (jj=jp; jj<jp+njp; jj++) {
    nkp = DKS(i,j)/DKS(ii,jj);
    kp = k*nkp;
    for (kk=kp; kk<kp+nkp; kk++) {
      vol_frac = ND_ELEM_LINEAR(g,ii,jj,kk).volume*vol_sum;
      VLOOP { pp[vv] += NDP_ELEM_LINEAR(p,ii,jj,kk,vv)*vol_frac; }
      pp[vv++] += NDP_ELEM_LINEAR(e,ii,jj,kk,PRESS)*vol_frac;
      pp[vv  ] += NDP_ELEM_LINEAR(p,ii,jj,kk,GAMMA)*vol_frac;
    }
  }
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void restrict1(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double ** p, double ** e, zone_geom * g, int i, int j, int k, int jj,
  double *pp)
#else
void restrict1(int ND_PTR ijk_to_I,int nvars,int nhydro,int ninterp,
		              int *dj,double **alpha1s,double **beta1s,double **Gamma1s,
			                    int **dk,double **alpha2s,double **beta2s,double **Gamma2s,
					    double NDP_PTR p, double NDP_PTR e, zone_geom ND_PTR g, int i, int j, int k, int jj,
  double *pp)
#endif
{
  int nkp,kp,kk,vv;
  double vol_sum,vol_frac;

  nkp = MIN(DKS(i,j),DKS(i,j-1))/DKS(i,jj);
  kp = k*nkp;
  vol_sum = 0.0;
  for (kk=kp; kk<kp+nkp; kk++) {
    vol_sum += ND_ELEM_LINEAR(g,i,jj,kk).volume;
  }
  vol_sum = 1.0/vol_sum;

  ILOOP { pp[vv] = 0.0; }
  for (kk=kp; kk<kp+nkp; kk++) {
    vol_frac = ND_ELEM_LINEAR(g,i,jj,kk).volume*vol_sum;
    VLOOP { pp[vv] += NDP_ELEM_LINEAR(p,i,jj,kk,vv)*vol_frac; }
    pp[vv++] += NDP_ELEM_LINEAR(e,i,jj,kk,PRESS)*vol_frac;
    pp[vv  ] += NDP_ELEM_LINEAR(e,i,jj,kk,GAMMA)*vol_frac;
  }
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void prim_to_fluxes3(double ** p)
#else
void prim_to_fluxes3(double NDP_PTR p)
#endif
{
  static int firstc = 1;
  int i,i0,ii,j,jj,k,kk,dd,vv,g,size;
  int jm,jp,njp,km,kp,nkp,isshock,jlev,klev,s,dka;
  int jstart,jstop,kstart,kstop;
  double ss,vol_sum,vol_frac;
  // double u[nhydro];

  TIMER_START("prim_to_fluxes3");

  if (firstc) {
    {
      //pl     = malloc_rank1(ninterp,      sizeof *pl   );
      //pr     = malloc_rank1(ninterp,      sizeof *pr   );
      //pp     = malloc_rank1(ninterp,      sizeof *pp   );
      //alpha  = malloc_rank1(max_grid_dim, sizeof *alpha);
      ////alpha += NG;
      //beta   = malloc_rank1(max_grid_dim, sizeof *beta );
      ////beta  += NG;
      //Gamma  = malloc_rank1(max_grid_dim, sizeof *Gamma);
      ////Gamma += NG;
      //x      = malloc_rank1(max_grid_dim, sizeof *x    );
      ////x     += NG;
    }
		firstc = 0;
  }

//  for(int II=0;II<cell_count_all;II++){
//    double pl[NSIZE],pr[NSIZE],pp[NSIZE];
//    double alpha[NSIZE_GRID],beta[NSIZE_GRID],Gamma[NSIZE_GRID],x[NSIZE_GRID];
//    double flux[NSIZE],vriemann[NSIZE];
//    GET_IJK_FROM_I(II,i,j,k);
//    if(i<istart[0]||i>istop[0]||j<JS(i,istart[1])||j>=JS(i,istop[1])||k<KS(i,j,istart[2])||k>=KS(i,j,istop[2])){/*do nothing*/}
//    else{
//        for (ii=i-NG; ii<i+NG; ii++) {
//          alpha[NG+ii-i] = alpha0[ii];
//          beta [NG+ii-i] = beta0 [ii];
//          Gamma[NG+ii-i] = Gamma0[ii];
//          x    [NG+ii-i] = startx[0] + ii*dx[0];
//
//          if (DJS(ii) < DJS(i)) {
//            restrict0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,sim_eos,geom,i,j,k,ii,pp);
//          } else {
//            prolongate0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,sim_eos,i,j,k,ii,pp);
//          }
//          ILOOP { pencil[vv][NG+ii-i] = pp[vv]; }
//        }
//
//        interp(ninterp,interp_order,pencil, pleft, pright, alpha, beta, Gamma, x, 0);  // interpolate a single 4-zone pencil
//
//        ILOOP {
//          pl[vv] = pright[vv][NG-1];
//          pr[vv] = pleft [vv][NG+0];
//        }
//    }
//  }
//  GPU_PRAGMA(omp target ){
//  for(int II=0;II<cell_count_all;II++){
//	  int Itemp;
//	  GET_IJK_FROM_I(II,i,j,k)
//	  Itemp = ND_ELEM(ijk_to_I,i,j,k);
//	  if(Itemp!=II){printf("error! %d %d %d %d\n",II,i,j,k);}
//	  for(int l=0;l<3;l++){
//	    printf("%d %d %d %d %.17g\n",i,j,k,l,ND_ELEM_LINEAR(geom,i,j,k).scale[3][l]);
//	  }
//  }
//  }
  GPU_PRAGMA(omp target teams distribute parallel for shared(ijk_to_I,I_to_ijk,sim_fdir0,sim_vedgedir0,sim_fdir1,sim_vedgedir1,sim_fdir2,sim_vedgedir2) reduction(min:min_dt))
  for(int II=0;II<cell_count_all;II++){
    double pl[NSIZE],pr[NSIZE],pp[NSIZE];
    double alpha[NSIZE_GRID],beta[NSIZE_GRID],Gamma[NSIZE_GRID],x[NSIZE_GRID];
    double flux[NSIZE],vriemann[NSIZE];
    double pencil[NSIZE][2*NG];
    double pleft [NSIZE][2*NG];
    double pright[NSIZE][2*NG];
    GET_IJK_FROM_I(II,i,j,k);
    if(i<istart[0]||i>istop[0]||j<JS(i,istart[1])||j>=JS(i,istop[1])||k<KS(i,j,istart[2])||k>=KS(i,j,istop[2])){/*do nothing*/}
    else{
        for (ii=i-NG; ii<i+NG; ii++) {
          alpha[NG+ii-i] = alpha0[ii];
          beta [NG+ii-i] = beta0 [ii];
          Gamma[NG+ii-i] = Gamma0[ii];
          x    [NG+ii-i] = startx[0] + ii*dx[0];

          if (DJS(ii) < DJS(i)) {
            restrict0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,sim_eos,geom,i,j,k,ii,pp);
          } else {
            prolongate0(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,sim_eos,i,j,k,ii,pp);
          }
          ILOOP { pencil[vv][NG+ii-i] = pp[vv]; }
        }

        interp(ninterp,interp_order,pencil, pleft, pright, alpha, beta, Gamma, x, 0);  // interpolate a single 4-zone pencil

        ILOOP {
          pl[vv] = pright[vv][NG-1];
          pr[vv] = pleft [vv][NG+0];
        }

        //isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 0, ND_ELEM_LINEAR_F(sim_fdir0,0,i,j,k), &ss, ND_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k));
        isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 0, flux, &ss, vriemann);
        HLOOP{NDP_ELEM_LINEAR_F(sim_fdir0,0,i,j,k,vv) = flux[vv];}
        SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir0,0,i,j,k,dd) = vriemann[dd];}
        if (isshock) {
          ND_ELEM_LINEAR(sim_shock_flag,i,j,k) = 0;
          jm = j*DJS(i)/DJS(i-1);
          km = k*DKS(i,j)/DKS(i-1,jm);
          ND_ELEM_LINEAR(sim_shock_flag,i-1,jm,km) = 0;
          jm = j*DJS(i)/DJS(i-2);
          km = k*DKS(i,j)/DKS(i-2,jm);
          ND_ELEM_LINEAR(sim_shock_flag,i-2,jm,km) = 0;
        }

        // update estimate for next time step, if necessary
        min_dt = MIN(min_dt,dx[0]/fabs(ss));
    }
    if(i<istart[0]||i>=istop[0]||j<JS(i,istart[1])||j>JS(i,istop[1])||k<KS(i,j,istart[2])||k>=KS(i,j,istop[2])){/*do nothing*/}
    else{
          jlev = 0;
          s = DJS(i);
          while (s >>= 1) jlev++;
          nkp = MAX(1,DKS(i,j)/MIN(DKS(i,j),DKS(i,j-1)));
          kp = k*MAX(1,DKS(i,j)/DKS(i,j-1));
          for (kk=kp; kk<kp+nkp; kk++) {
            for (jj=j-NG; jj<j+NG; jj++) {
              // Get interpolation coefficients
              alpha[NG+jj-j] = alpha1s[jlev][jj];
              beta [NG+jj-j] = beta1s [jlev][jj];
              Gamma[NG+jj-j] = Gamma1s[jlev][jj];
              x    [NG+jj-j] = startx[1] + jj*DJS(i)*dx[1];

              //if (DKS(i,jj) < DKS(i,j)) {
              if (DKS(i,jj) < MIN(DKS(i,j),DKS(i,j-1))) {//test required: should do the same thing as the line commented above
                restrict1(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,sim_eos,geom,i,j,kk,jj,pp);
              } else {
                prolongate1(ijk_to_I,nvars,nhydro,ninterp,dj,alpha1s,beta1s,Gamma1s,dk,alpha2s,beta2s,Gamma2s,p,sim_eos,i,j,kk,jj,pp);
              }
              ILOOP { pencil[vv][NG+jj-j] = pp[vv]; }
            }

            interp(ninterp,interp_order,pencil, pleft, pright, alpha, beta, Gamma, x, 0);  // interpolate a single 4-zone pencil

            ILOOP {
              pl[vv] = pright[vv][NG-1];
              pr[vv] = pleft [vv][NG+0];
            }

            if (transverse_shock(ijk_to_I,sim_shock_flag,dj,dk,i,j,k,1)) {
              //riemann_flux_HLLE(pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 1, ND_ELEM_LINEAR_F(sim_fdir1,1,i,j,k), &ss, ND_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k));
              riemann_flux_HLLE(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 1, flux, &ss, vriemann);
              HLOOP{NDP_ELEM_LINEAR_F(sim_fdir1,1,i,j,kk,vv) = flux[vv];}
              SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,kk,dd) = vriemann[dd];}
            } else {
              //riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 1, ND_ELEM_LINEAR_F(sim_fdir1,1,i,j,k), &ss, ND_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,k));
              isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 1, flux, &ss, vriemann);
              HLOOP{NDP_ELEM_LINEAR_F(sim_fdir1,1,i,j,kk,vv) = flux[vv];}
              SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir1,1,i,j,kk,dd) = vriemann[dd];}
            }

            // update estimate for next time step, if necessary
            min_dt = MIN(min_dt,DJS(i)*dx[1]/fabs(ss));
	  }

    }  // end if
    if(i<istart[0]||i>=istop[0]||j<JS(i,istart[1])||j>=JS(i,istop[1])||k<KS(i,j,istart[2])||k>KS(i,j,istop[2])){/*do nothing*/}
    else{
        klev = 0;
        s = DKS(i,j);
        while (s >>= 1) klev++;
        for (int ktemp=k-NG; ktemp<k+NG; ktemp++) {
          x    [NG+ktemp-k] = startx[2] + ktemp*DKS(i,j)*dx[2];
          Gamma[NG+ktemp-k] = Gamma2s[klev][ktemp];
          beta [NG+ktemp-k] = beta2s [klev][ktemp];
          alpha[NG+ktemp-k] = alpha2s[klev][ktemp];
          VLOOP { pencil[vv][NG+ktemp-k] = NDP_ELEM_LINEAR(p,i,j,ktemp,vv); }
          pencil[vv++][NG+ktemp-k] = NDP_ELEM_LINEAR(sim_eos,i,j,ktemp,PRESS);
          pencil[vv  ][NG+ktemp-k] = NDP_ELEM_LINEAR(sim_eos,i,j,ktemp,GAMMA);
        }

        // interp() computes the reconstruction in n+2 zones, where n is the last argument of interp()
        // The reconstruction requires n+4 total zones
        // For every pair of L/R states, we need the reconstruction in 2 zones
        // For n zones, there are n+1 interfaces, requiring reconstructions in n+2 total zones
        // n=kstop-kstart gives the total number of zones
        interp(ninterp,interp_order,pencil, pleft, pright, alpha, beta, Gamma, x, 0);  // interpolate kstop-kstart+1 4-zone pencils

        #if 1
        // EXPERIMENTAL CODE TO TRY TO FIX THE POLAR AXIS PROBLEM

        // if (abs(j) < 2 || abs(j-(JS(i,global_grid_dims[1])-1)) < 2) {
        if (abs(j-(JS(i,global_grid_dims[1])-1)) < 1) {
          //int l,m;
          double r,th,sth,cth,ph,sph,cph,vx,vy;
          double lam[3][3],v[3];

          jlev = 0;
          s = DJS(i);
          while (s >>= 1) jlev++;
          r = r_of_x(rx_info,beta0[i]/Gamma0[i]);

          th = th_of_x(thx_info,beta1s[jlev][j]/Gamma1s[jlev][j]);
          sth = sin(th);
          cth = cos(th);

          for (int ktemp=k-NG; ktemp<k+NG; ktemp++) {
            ph  = startx[2] + (ktemp+0.5)*DKS(i,j)*dx[2];
            sph = sin(ph);
            cph = cos(ph);

            // Compute transformation matrix at zone volume centers
            lam[0][0] = sth*cph;
            lam[0][1] = cth*cph;
            lam[0][2] = -sph;  // /sth ?

            lam[1][0] = sth*sph;
            lam[1][1] = cth*sph;
            lam[1][2] = cph;   // /sth?

            lam[2][0] = cth;
            lam[2][1] = -sth;
            lam[2][2] = 0.0;

            // Convert to Cartesian vector components
            for (int l=0; l<3; l++) {
              v[l] = 0.0;
              for (int m=0; m<3; m++) {
                v[l] += lam[l][m]*NDP_ELEM_LINEAR(p,i,j,ktemp,U1+m)*ND_ELEM_LINEAR(geom,i,j,ktemp).scale[0][m];
              }
            }


            ph  = startx[2] + ktemp*DKS(i,j)*dx[2];
            sph = sin(ph);
            cph = cos(ph);

            // Compute transformation matrix at lower phi-interface centers
            lam[0][0] = sth*cph;
            lam[0][1] = cth*cph;
            lam[0][2] = -sph;

            lam[1][0] = sth*sph;
            lam[1][1] = cth*sph;
            lam[1][2] = cph;

            lam[2][0] = cth;
            lam[2][1] = -sth;
            lam[2][2] = 0.0;

            ILOOP { pleft[vv][NG+ktemp-k] = pencil[vv][NG+ktemp-k]; }
            for (int l=0; l<3; l++) {
              pleft[U1+l][NG+ktemp-k] = 0.0;
              for (int m=0; m<3; m++) {
                // The transformation matrix is orthogonal, so its inverse is its transpose
                pleft[U1+l][NG+ktemp-k] += lam[m][l]*v[m]/ND_ELEM_LINEAR(geom,i,j,ktemp).scale[3][l];
              }
            }

	    if(ktemp>=KS(i,j,istop[2])+NG-1){continue;}

            ph  = startx[2] + (ktemp+1.0)*DKS(i,j)*dx[2];
            sph = sin(ph);
            cph = cos(ph);

            // Compute transformation matrix at upper phi-interface centers
            lam[0][0] = sth*cph;
            lam[0][1] = cth*cph;
            lam[0][2] = -sph;

            lam[1][0] = sth*sph;
            lam[1][1] = cth*sph;
            lam[1][2] = cph;

            lam[2][0] = cth;
            lam[2][1] = -sth;
            lam[2][2] = 0.0;

            ILOOP { pright[vv][NG+ktemp-k] = pencil[vv][NG+ktemp-k]; }
            for (int l=0; l<3; l++) {
              pright[U1+l][NG+ktemp-k] = 0.0;
              for (int m=0; m<3; m++) {
                // The transformation matrix is orthogonal, so its inverse is its transpose
                pright[U1+l][NG+ktemp-k] += lam[m][l]*v[m]/ND_ELEM_LINEAR(geom,i,j,ktemp+1).scale[3][l];
              }
            }
          }
        }
        #endif

        ILOOP {
          pl[vv] = pright[vv][NG-1];
          pr[vv] = pleft [vv][NG+0];
        }

        if (transverse_shock(ijk_to_I,sim_shock_flag,dj,dk,i,j,k,2)) {
          //riemann_flux_HLLE(pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 2, ND_ELEM_LINEAR_F(sim_fdir2,2,i,j,k), &ss, ND_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k));
          riemann_flux_HLLE(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 2, flux, &ss, vriemann);
          HLOOP{NDP_ELEM_LINEAR_F(sim_fdir2,2,i,j,k,vv) = flux[vv];}
          SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k,dd) = vriemann[dd];}
        } else {
          //riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 2, ND_ELEM_LINEAR_F(sim_fdir2,2,i,j,k), &ss, ND_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k));
          isshock = riemann_flux_HLLC(nhydro,nvars,pl, pr, &ND_ELEM_LINEAR(geom,i,j,k), 2, flux, &ss, vriemann);
          HLOOP{NDP_ELEM_LINEAR_F(sim_fdir2,2,i,j,k,vv) = flux[vv];}
          SLOOP{NDP_ELEM_LINEAR_F(sim_vedgedir2,2,i,j,k,dd) = vriemann[dd];}
        }

        // update estimate for next time step, if necessary
        min_dt = MIN(min_dt,DKS(i,j)*dx[2]/fabs(ss));
    }
  }
          
  firstc = 0;

  TIMER_STOP;

  return;
}
#endif
