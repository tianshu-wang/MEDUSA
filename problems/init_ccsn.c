#include "../decs.h"
#include "../constants.h"
#if (NDIM>1 && PERTURB==VELOCITY_SPHERICAL_HARMONIC)
#include "gsl/gsl_sf_legendre.h"
#endif /* PERTURB==VELOCITY_SPHERICAL_HARMONIC */

#if (DO_RADIATION!=FALSE)
#error This problem must be compiled with DO_RADIATION=FALSE
#endif

#if (NEUTRINO!=FALSE)
#error This problem must be compiled with NEUTRINO=FALSE
#endif

#if (GRAV!=SPHERICAL_MONOPOLE_GRAV && GRAV!=SPHERICAL_MULTIPOLE_GRAV)
#error This problem must be compiled with GRAV==SPHERICAL_MONOPOLE_GRAV or GRAV==SPHERICAL_MULTIPOLE_GRAV
#endif

#if (GEOM!=SPHERICAL)
#error This problem must be compiled with GEOM=SPHERICAL
#endif

#if (EOS!=COLLAPSE)
#error This problem must be compiled with EOS=COLLAPSE
#endif

//#if (USE_EXT_SRC!=TRUE)
//#error This problem must be compiled with USE_EXT_SRC=TRUE
//#endif
//
//#if (USE_AUX_SRC!=TRUE)
//#error This problem must be compiled with USE_AUX_SRC=TRUE
//#endif

static double rho_out, ye_out, u_out, u1_out;
double psi(int l, int n, double r, double th, double rmin, double rmax);
void perturb(int l, int n, double vmax, double rmin, double rmax);
double ran2(long int *idum);
const double ep=1.0e-6;
const double km=1.0e5;
double ds_max = 0.0;
static int perturbed=0;
static long int rseed;

#if USE_FORTRAN_CODE
void FORT_EOS_CALL(eos_given_rex)(double *Gamma, double *press, double *cs, double *temp, double *ent, double *dpdr, double *dpde, double *rho, double *e, double *ye, int *pt);
void FORT_EOS_CALL(eos_given_rtx)(double *e, double *press, double *gam, double *cs, double *rho, double *temp, double *ye, int *pt);
endif
#endif

void init_grid()
{
  startx[0] = 0.0;
  startx[1] =-1.0;
  startx[2] = 0.0;
  dx[0]     = ccsn_dr_min;  // This is dr_min at r=0
  dx[1]     = 2.0/n2;  // = (mu_max - mu_min)/N_th, where mu = -cos(th)
  dx[2]     = 2.0*M_PI/n3;  // = (ph_max - ph_min)/N_ph

  rx_info.rtrans = rtrans_solve(dx[0], outer_radius);  // approx transition radius from const to log spacing
  if (mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rx_info.rtrans);

  r_full_res = 1.0e7;

  periodic[0] = 0;
  periodic[1] = 0;
  periodic[2] = 1;

  return;
}

double model_interp(double r, double r0, double r1, double p0, double p1)
{
  double A = (p1 - p0)/(SQR(r1) - SQR(r0));
  double B = p0 - A*SQR(r0);

  return A*SQR(r) + B;
}

void init_problem()
{
  int ii=0,jj=0,kk=0,vv,i,j,idum,lines,g,dd;
  int *pt_index=NULL;
  char ch,dum[1024];
  double r,del,dummy,R,th;
  double rho,ye,temp,v,eint,total_vol,total_mass;
  double* r_model=NULL;
  double* rho_model=NULL;
  double* v_model=NULL;
  double* ye_model=NULL;
  double* temp_model=NULL;
  double* s_model=NULL;
  double* lapse_model=NULL;
  double* lapse_edge_model=NULL;
  double smallt = 1.0e-4;
  double smalld = 1.0e-4;
  double gam = 4.0/3.0;
  double cxlo[3],cxhi[3];
  FILE *fp=NULL;

  double get_model_u(double rho, double T, double ye);
#if USE_FORTRAN_CODE
  void FORT_EOS_CALL(eos_init)(char *name, double *small_temp, double *small_dens, double *gam_in, size_t name_len);

  FORT_EOS_CALL(eos_init)(eos_file, &smallt, &smalld, &gam, strlen(eos_file));
#else
  eos_init(eostable,eos_file);
#endif
  rseed = -(myrank+1);
  // set the boundary conditions
  VLOOP {
    bc[vv].lo[0] = SPHERICAL_ORIGIN;
    #if (NDIM>1)
    bc[vv].lo[1] = SPHERICAL_ORIGIN;
    bc[vv].hi[1] = SPHERICAL_ORIGIN;
    #endif
    #if (NDIM==3)
    bc[vv].lo[2] = PERIODIC;
    bc[vv].hi[2] = PERIODIC;
    #endif
  }
  HLOOP {
    bc[vv].hi[0] = OUTFLOW;
  }

  M_prescribed = 0.0;

  if (restart_from_hdf == TRUE) {
    restart_read();
    ZGLOOP { NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP) = -1.0; }
  }
  #if (NDIM==3)
  else if(restart_from_3d == TRUE) {
    // Read 3D data
    //read3d()
    read3(model_name);

    // Initialize lapse function
    #if (GR_MONOPOLE==TRUE)
    ISGLOOP(ii) {
      gr_lapse[ii] = 1.0;
      gr_lapse_edge[ii] = 1.0;
    }
    #endif

    // Initialize radiation to zero
    ZGLOOP {
      GLOOP {
        NDP_ELEM_LINEAR(sim_p,ii,jj,kk,irad1+g) = 0.0;
        DLOOP {
          NDP_ELEM_LINEAR(sim_p,ii,jj,kk,ifrad1+g*NDIM+dd) = 0.0;
        }
      }
    }
  }
  #endif
  else {
    if (mpi_io_proc()) {
      fp = fopen(model_name,"r");
      if (fp==NULL) {
        fprintf(stderr,"[init_problem]:  problem opening the 1-d model %s!\n",model_name);
        fflush(stderr);
        exit(1);
      }

      // Count & verify the number of lines in the 1d dump
      lines = 0;
      if (restart_from_1d) fscanf(fp,"%[^\n]%*c", &dum[0]);  // skip over the header line
      while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch=='\n') lines++;
      }
      if (restart_from_1d && lines != n1) {
        fprintf(stderr,"[init_problem]:  Found %d lines in file %s, but require %d\n",lines,model_name,n1);
        fflush(stderr);
        exit(1);
      }
      if (mpi_io_proc()) printf("[init_problem]:  Found %d lines in file %s\n",lines,model_name);
    }

    #if (USE_MPI==TRUE)
    MPI_Bcast(&lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
    #endif

    // Allocate 1d model vectors
    r_model     = malloc_rank1(lines, sizeof *r_model    );
    rho_model   = malloc_rank1(lines, sizeof *rho_model  );
    v_model     = malloc_rank1(lines, sizeof *v_model    );
    ye_model    = malloc_rank1(lines, sizeof *ye_model   );
    temp_model  = malloc_rank1(lines, sizeof *temp_model );
    s_model     = malloc_rank1(lines, sizeof *s_model    );
    lapse_model      = malloc_rank1(lines, sizeof *lapse_model     );
    lapse_edge_model = malloc_rank1(lines, sizeof *lapse_edge_model);

    if (mpi_io_proc()) {
      rewind(fp);

      if (restart_from_1d) fscanf(fp,"%[^\n]%*c", &dum[0]);  // skip over the header line
      for (i=0; i<lines; i++) {
        if (restart_from_1d) {
          // r, RHO, UU, U1, U2, U3, YE, PRESS, CS, TEMP [MEV], ENT, GAMMA, MASS COORD
          fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            &r_model[i], &rho_model[i], &dummy, &v_model[i], &dummy, &dummy,
            &ye_model[i], &dummy, &dummy, &temp_model[i], &dummy, &dummy, &dummy);
          fscanf(fp,"%lf %lf ",&lapse_model[i],&lapse_edge_model[i]);
        }
        else {
          // i, mass, r, RHO, TEMP [K], YE, U1, 0
          fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n", &idum, &dummy, &r_model[i],
            &rho_model[i], &temp_model[i], &ye_model[i], &v_model[i], &dummy);
          temp_model[i] *= K_TO_MEV;
          // printf("myrank=%d, i=%d, r_model=%e, rho_model=%e, ye_model=%e, temp_model=%e, v_model=%e\n",
          //   myrank,i,r_model[i],rho_model[i],ye_model[i],temp_model[i],v_model[i]);
        }
      }
      fprintf(stderr, "done!\n");
      fflush(stderr);
      fclose(fp);
    }

    #if (USE_MPI==TRUE)
    MPI_Bcast(&lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(r_model,    lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho_model , lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(temp_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(s_model,    lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ye_model,   lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v_model,    lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lapse_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lapse_edge_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif

    ZGLOOP { NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP) = -1.0; }

    ISLOOP(ii) {
      r = r_of_x(rx_info,startx[0] + (ii+0.5)*dx[0]);
      i = 0;
      while (r > r_model[i] && i<lines) i++;
      i--;
      if (i<0) {
        rho  = model_interp(r, r_model[0], r_model[1],  rho_model[0],  rho_model[1]);
        ye   = model_interp(r, r_model[0], r_model[1],   ye_model[0],   ye_model[1]);
        temp = model_interp(r, r_model[0], r_model[1], temp_model[0], temp_model[1]);
        v    = r/r_model[0]*v_model[0];
      } else {
        del  = (r-r_model[i])/(r_model[i+1]-r_model[i]);
        rho  = (1.0-del)* rho_model[i] + del* rho_model[i+1];
        ye   = (1.0-del)*  ye_model[i] + del*  ye_model[i+1];
        temp = (1.0-del)*temp_model[i] + del*temp_model[i+1];
        v    = (1.0-del)*   v_model[i] + del*   v_model[i+1];
      }
      // printf("myrank=%d, ii=%d, r=%e, rho=%e, ye=%e, temp=%e, v=%e\n",myrank,ii,r,rho,ye,temp,v);
      #if (NDIM>1)
      JSLOOP(ii,jj) {
      #endif
        #if (NDIM==3)
        KSLOOP(ii,jj,kk) {
        #endif
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,RHO ) = rho;
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,UU  ) = get_model_u(rho,temp,ye);
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,U1  ) = v/ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0];
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,U2  ) = 0.0;
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,U3  ) = 0.0;
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,YE  ) = ye;
          NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP) = temp;
          #if (GR_MONOPOLE==TRUE)
          gr_lapse[ii]      = (restart_from_1d) ? lapse_model[ii]      : 1.0;
          gr_lapse_edge[ii] = (restart_from_1d) ? lapse_edge_model[ii] : 1.0;
          #endif
          #if (DO_ROTATION==TRUE)
          /* Rotational profile is Omega(R) = Omega0/(1+(R/A)^2), where R is the cylindrical radius */
          th = th_of_x(startx[1] + (jj+0.5)*DJS(ii)*dx[1]);
          R = r*sin(th);
          NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U3) = R*rotate_Omega0/(1.0 + pow(R/(rotate_A*km),2));
          NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U3) /= ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][2];
          #endif /* DO_ROTATION==TRUE */
        #if (NDIM==3)
        }
        #endif
      #if (NDIM>1)
      }
      #endif
    }

    #if (PERTURB==VELOCITY_SPHERICAL_HARMONIC)
    /* AS:  Perturb the velocity as in Muller & Janka (2014) */
    perturb(perturb_l1, perturb_m1, perturb_n1, perturb_dv1, perturb_r1m, perturb_r1p);
    if (mpi_io_proc()) {
      fprintf(stderr, "Perturbing:  l=%d, m=%d, n=%d, dv_max=%g km/s, rmin=%g km, rmax=%g km\n",
          perturb_l1, perturb_m1, perturb_n1, perturb_dv1/km, perturb_r1m/km, perturb_r1p/km);
      fflush(stderr);
    }
    perturb(perturb_l2, perturb_m2, perturb_n2, perturb_dv2, perturb_r2m, perturb_r2p);
    if (mpi_io_proc()) {
      fprintf(stderr, "Perturbing:  l=%d, m=%d, n=%d, dv_max=%g km/s, rmin=%g km, rmax=%g km\n",
          perturb_l2, perturb_m2, perturb_n2, perturb_dv2/km, perturb_r2m/km, perturb_r2p/km);
      fflush(stderr);
    }
    perturb(perturb_l3, perturb_m3, perturb_n3, perturb_dv3, perturb_r3m, perturb_r3p);
    if (mpi_io_proc()) {
      fprintf(stderr, "Perturbing:  l=%d, m=%d, n=%d, dv_max=%g km/s, rmin=%g km, rmax=%g km\n",
          perturb_l3, perturb_m3, perturb_n3, perturb_dv3/km, perturb_r3m/km, perturb_r3p/km);
      fflush(stderr);
    }
    #endif /* PERTURB==VELOCITY_SPHERICAL_HARMONIC */

    #if (NDIM>1 && PERTURB==DENSITY_RANDOM)
    if (perturb_delay<0.0) {
      /* AS:  Immediately perturb the density as in Summa et al. (2015) */
      if (mpi_io_proc()) fprintf(stderr,"[init_problem]:  Perturbing density...\n");
      ZLOOP { NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO) *= 1.0 + (perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0); }
      perturbed = 1;
    }
    #endif /* PERTURB==DENSITY_RANDOM */

    free(r_model);
    free(rho_model);
    free(temp_model);
    free(ye_model);
    free(v_model);
    free(s_model);
    free(lapse_model);
    free(lapse_edge_model);

    init_tracers_by_mass(&rseed);

  }

  // Calculate outer BCs
  ii = istop[0]-1;
  jj = JS(ii,istart[1]);
  kk = KS(ii,jj,istart[2]);
  rho_out = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
  u_out   = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,UU );
  u1_out  = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1 );
  ye_out  = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,YE );

  if (mpi_io_proc()) {
    fprintf(stderr,"rho_out = %23.16e\n",rho_out);
    fprintf(stderr,"u_out   = %23.16e\n",u_out);
    fprintf(stderr,"u1_out  = %23.16e\n",u1_out*ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0]);
    fprintf(stderr,"ye_out  = %23.16e\n",ye_out);
    fflush(stderr);
  }

  #if (GR_MONOPOLE==TRUE)
  /* copy 1-d lapses across grid */
  gr_lapse_edge[n1] = gr_lapse_edge[n1-1];
  ZLOOP {
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
  }
  ii = istop[0];
  #if (NDIM>1)
  JSLOOP(ii,jj) {
  #endif
    #if (NDIM==3)
    KSLOOP(ii,jj,kk) {
    #endif
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
    #if (NDIM==3)
    }
    #endif
  #if (NDIM>1)
  }
  #endif

  #if (NDIM==2)
  ISLOOP(ii) {
    jj = JS(ii,istop[1]);
    ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
  }
  #endif

  #if (NDIM==3)
  ISLOOP(ii) {
    jj = JS(ii,istop[1]);
    KSLOOP(ii,jj,kk) {
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
    JSLOOP(ii,jj) {
      kk = KS(ii,jj,istop[2]);
      ND_ELEM_LINEAR(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
  }
  #endif
  #endif /* GR_MONOPOLE==TRUE */

  reset_boundaries(sim_p,0);
  sync_mpi_boundaries(sim_p); //cpu, mpi communication
  complete_mpi_communication(0);

  return;
}


void prob_bounds(int i, int j, int k, double *p)
{
  //p[RHO] = rho_out;
  //p[UU]  = u_out;
  //p[U1]  = u1_out;
  //p[YE]  = ye_out;

  return;
}

#if (DEBUG==TRUE)
const int check_ghost = 1;
const int verbose = 1;
const int check_shells = 0;
#endif
void analysis_preloop()
{
  #if (DEBUG==TRUE)
  calculate_total_mass(verbose);
  calculate_total_energy(verbose);

  int status = find_prim_errors(sim_p,check_ghost,"analysis_preloop");
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE,&status,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  #endif
  if (status) {
    if (mpi_io_proc()) fprintf(stderr,"[analysis_preloop]:  Error detected!  Exiting...\n");
    fflush(stderr);
    exit(1);
  } else {
    if (mpi_io_proc()) fprintf(stderr,"[analysis_preloop]:  No errors!\n");
  }
  #endif /* DEBUG */

  return;
}

void analysis_inloop()
{
  #if (DEBUG==TRUE)
  calculate_total_mass(verbose);
  calculate_total_energy(verbose);
  check_flux_differencing(verbose,check_shells);

  int status = find_prim_errors(sim_p,check_ghost,"analysis_inloop");
  #if (USE_MPI==TRUE)
  MPI_Allreduce(MPI_IN_PLACE,&status,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  #endif
  if (status) {
    if (mpi_io_proc()) fprintf(stderr,"[analysis_inloop]:  Error detected!  Exiting...\n");
    fflush(stderr);
    exit(1);
  } else {
    if (mpi_io_proc()) fprintf(stderr,"[analysis_inloop]:  No errors!\n");
  }
  #endif /* DEBUG */

  return;
}

void analysis_postloop()
{
  return;
}

double get_model_u(double rho, double Temp, double ye)
{
  int *pt=NULL;
  double e,press,gam,cs;

  if (Temp < 1.0e-2) fprintf(stderr,"WTF Temp < 1e-2!: %g\n", Temp);
  return u_given_rtx(eostable,rho, Temp, ye);
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void aux_src(double ** p, double dtpush)
#else
void aux_src(double NDP_PTR p, double dtpush)
#endif
{
  int ii,jj,kk;
  double yefloor = 0.27;
  double a0 =  2.120282875020e2;
  double a1 = -1.279083733386e2;
  double a2 =  3.199984467460e1;
  double a3 = -4.238440349145;
  double a4 =  3.133904196302e-1;
  double a5 = -1.226365543366e-2;
  double a6 =  1.983947360151e-4;
  double rho,lrho,ye,yebar,dye;

  //ZLOOP {
  GPU_PRAGMA(omp target teams distribute parallel for)
  for(int II=0;II<cell_count;II++){
    GET_IJK_FROM_I(II,ii,jj,kk);
    rho = NDP_ELEM_LINEAR(p,ii,jj,kk,RHO);
    ye  = NDP_ELEM_LINEAR(p,ii,jj,kk,YE );
    if (rho > 1.0e8 && ye > yefloor) {
      lrho = log10(rho);
      yebar = a0 + lrho*(a1 + lrho*(a2 + lrho*(a3 + lrho*(a4 + lrho*(a5 + lrho*a6)))));
      dye = MIN(0.0, yebar - ye);
      dye = MAX(-0.05*ye, dye);
      if (rho < 3.0e8) {
        dye *= (rho-1.0e8)/2.0e8;
      }
      NDP_ELEM_LINEAR(sim_src,ii,jj,kk,YE) += rho*dye/dtpush;
    }
  }

  if (!perturbed && !restart_from_hdf && !restart_from_1d && !restart_from_3d) {
    if (tbounce>0.0 && perturb_delay>=0.0 && t>=tbounce+perturb_delay) {
      #if (NDIM>1 && PERTURB==DENSITY_RANDOM)
      /* Implements the perturbation scheme of Summa et al. (2015),
       * which perturbs the density once by some relative amount */
      if (mpi_io_proc()) fprintf(stderr,"[init_problem]:  Perturbing density...\n");
      rseed = -myrank;
      //ZLOOP { 
      GPU_PRAGMA(omp target teams distribute parallel for)
      for(int II=0;II<cell_count;II++){
	      GET_IJK_FROM_I(II,ii,jj,kk);
	      NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO) *= 1.0 + (perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0); 
      }
      #endif /* PERTURB==DENSITY_RANDOM */
      perturbed=1;
    }
  }
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void ext_src(double ** p, double dtpush)
#else
void ext_src(double NDP_PTR p, double dtpush)
#endif
{
  int ii,jj,kk;
  static int heaton=0;
  int *pt=NULL;
  double max_rho;
  double yl1,yl2,yl3,yl4,xl1,xl2,xl3,xl4,sl1,sl2,sl3;
  double H0,C0,lrho,x1,y1,s,tlog,tau,eta,xn,xp;
  double heat, cool;
  double *table = eostable;
  //void FORT_EOS_CALL(eos_get_xnxp)(double *xn, double *xp, double *rho, double *T, double *ye, int *pt);

  yl1 = -4.0;
  yl2 = -1.8;
  yl3 = 2.4;
  yl4 = 3.4;
  xl1 = 8.35;
  xl2 = 8.8;
  xl3 = 14.15;
  xl4 = 14.5;
  sl1 = (yl2-yl1)/(xl2-xl1);
  sl2 = (yl3-yl2)/(xl3-xl2);
  sl3 = (yl4-yl3)/(xl4-xl3);

  H0 = 1.544e20 * (L_enu/1.0e52) * pow(T_enu/4.0, 2.0);
  C0 = 1.399e20;

  if (!heaton) {
    ii = istart[0];
    jj = JS(ii,istart[1]);
    kk = KS(ii,jj,istart[2]);
    max_rho = NDP_ELEM_LINEAR(p,ii,jj,kk,RHO);
    max_rho = mpi_max(max_rho);
    if (max_rho > 1.0e14) {
      heaton = 1;
    }
  }

  //ZLOOP {
  GPU_PRAGMA(omp target teams distribute parallel for)
  for(int II=0;II<cell_count;II++) {
    GET_IJK_FROM_I(II,ii,jj,kk);
    lrho = log10(NDP_ELEM_LINEAR(p,ii,jj,kk,RHO));
    if (lrho < xl2) {
      y1 = yl1;
      x1 = xl1;
      s = sl1;
    } else if (lrho > xl3) {
      y1 = yl3;
      x1 = xl3;
      s = sl3;
    } else {
      y1 = yl2;
      x1 = xl2;
      s = sl2;
    }
    tlog = s*(lrho - x1) + y1;
    tau = pow(10.0, tlog) * pow(T_enu/4.0, 2.0);

    eos_get_etaxnxp(table,NDP_ELEM_LINEAR(p,ii,jj,kk,RHO), NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP), NDP_ELEM_LINEAR(p,ii,jj,kk,YE), &eta, &xn, &xp);

    heat = heaton * H0/pow(rcenter[ii]/1.0e7, 2.0);
    cool = C0 * pow(NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP)/2.0, 6.0);

    NDP_ELEM_LINEAR(sim_src,ii,jj,kk,ETOT) += (heat - cool) * NDP_ELEM_LINEAR(p,ii,jj,kk,RHO) * (xn + xp) * exp(-tau);
  }
}

#if (NDIM>1 && PERTURB==VELOCITY_SPHERICAL_HARMONIC)
double psi(int l, int n, double r, double th, double rmin, double rmax)
{
  return (sqrt(sin(th))/r)*sin(n*M_PI*(r-rmin)/(rmax-rmin))*gsl_sf_legendre_sphPlm(l,1,cos(th));
}

void perturb(int l, int m, int n, double v, double rmin, double rmax)
{
  int ii=0,jj=0,kk=0,firstc=1;
  double x[SPACEDIM],r,rm,rp,th,thm,thp,ph,rho,dv,dv_max=0.0,C,tmp;
  double ND_PTR dvr;
  double ND_PTR dvt;

  if (firstc) {
    dvr = dendritic_malloc_double();
    dvt = dendritic_malloc_double();
    firstc = 0;
  }

  ZLOOP {
    ijk_to_x(ii,jj,kk,x);
    r   = r_of_x(rx_info,x[0]);
    th  = th_of_x(x[1]);
    ph  = x[2];
    ND_ELEM(dvr,ii,jj,kk) = 0.0;
    ND_ELEM(dvt,ii,jj,kk) = 0.0;
    if (r>=rmin && r<=rmax) {
      rho = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
      rm  = r*(1.0-ep);
      rp  = r*(1.0+ep);
      thm = MAX(th*(1.0-ep),0.0);
      thp = MIN(th*(1.0+ep),M_PI);
      ND_ELEM(dvr,ii,jj,kk) = (1.0/rho)*(1.0/(r*sin(th)))*(sin(thp)*psi(l,m,n,r,thp,ph,rmin,rmax) - sin(thm)*psi(l,m,n,r,thm,ph,rmin,rmax))/(thp-thm);
      ND_ELEM(dvt,ii,jj,kk) = (1.0/rho)*(-1.0/r)*(rp*psi(l,m,n,rp,th,ph,rmin,rmax) - rm*psi(l,m,n,rm,th,ph,rmin,rmax))/(rp-rm);
      dv = sqrt(SQR(ND_ELEM(dvr,ii,jj,kk)) + SQR(ND_ELEM(dvt,ii,jj,kk)));
      dv_max = MAX(dv_max,dv);
    }
  }
  #if (USE_MPI==TRUE)
  MPI_Allreduce(&dv_max, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  dv_max = tmp;
  #endif

  C = (dv_max) ? v/dv_max : 0.0;
  ZLOOP {
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1) += C*ND_ELEM(dvr,ii,jj,kk)/ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0];
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U2) += C*ND_ELEM(dvt,ii,jj,kk)/ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][1];
  }
}
#endif /* PERTURB==VELOCITY_SPHERICAL_HARMONIC */
