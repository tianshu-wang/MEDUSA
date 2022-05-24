#include "../decs.h"
#include "../constants.h"
#if (NDIM>1 && PERTURB==VELOCITY_SPHERICAL_HARMONIC)
#include "gsl/gsl_sf_legendre.h"
#endif /* PERTURB==VELOCITY_SPHERICAL_HARMONIC */

#if (DO_RADIATION!=TRUE)
#error This problem must be compiled with DO_RADIATION=TRUE
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

static double rho_out, ye_out, u_out, u1_out;
double psi(int l, int m, int n, double r, double th, double ph, double rmin, double rmax);
void perturb(int l, int m, int n, double vmax, double rmin, double rmax);
double ran2(long int *idum);
const double ep=1.0e-6;
const double km=1.0e5;
double ds_max = 0.0;
static int perturbed=0;
static long int rseed;

#if USE_FORTRAN_CODE
void FORT_EOS_CALL(eos_given_rex)(double *Gamma, double *press, double *cs, double *temp, double *ent, double *dpdr, double *dpde, double *rho, double *e, double *ye, int *pt);
void FORT_EOS_CALL(eos_given_rtx)(double *e, double *press, double *gam, double *cs, double *rho, double *temp, double *ye, int *pt);
#else
#include "../eos/burrows/eos_stuff.h"
#endif

void init_grid()
{
  startx[0] = 0.0;
  startx[1] =-1.0;
  startx[2] = 0.0;
  dx[0]     = ccsn_dr_min;  // This is dr_min at r=0
  dx[1]     = 2.0/n2;  // = (mu_max - mu_min)/N_th, where mu = -cos(th)
  dx[2]     = 2.0*M_PI/n3;  // = (ph_max - ph_min)/N_ph

  rtrans = rtrans_solve(dx[0], outer_radius);  // approx transition radius from const to log spacing
  if (mpi_io_proc()) fprintf(stderr,"rtrans = %g\n", rtrans);

  periodic[0] = 0;
  periodic[1] = 0;
  periodic[2] = 1;

  return;
}

double model_interp(double r, double r0, double r1, double p0, double p1)
{
  double A = (p1 - p0)/(r1*r1 - r0*r0);
  double B = p0 - A*r0*r0;

  return A*r*r + B;
}

void init_problem()
{
  int ii=0,jj=0,kk=0,g,dd,vv,i,j,idum,lines,nrad1d;
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
  double* rad_model=NULL;
  double* lapse_model=NULL;
  double* lapse_edge_model=NULL;
  double smallt = 1.0e-4;
  double smalld = 1.0e-4;
  double gam = 4.0/3.0;
  double cxlo[3],cxhi[3];
  FILE *fp=NULL;
  double* Fcov=NULL;
  double* Fred=NULL;

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
    bc[vv].hi[0] = PROB;
  }
  bc[UU].hi[0] = OUTFLOW;
  bc[U2].hi[0] = OUTFLOW;
  bc[U3].hi[0] = OUTFLOW;
  for (vv=irad1; vv<irad1+ngroups; vv++) bc[vv].hi[0] = RADEXTRAP;
  for (vv=ifrad1; vv<ifrad1+ngroups*NDIM; vv+=NDIM) bc[vv].hi[0] = RADEXTRAP;
  for (vv=ifrad1; vv<ifrad1+ngroups*NDIM; vv+=NDIM) {
    for (dd=1; dd<NDIM; dd++) bc[vv+dd].hi[0] = OUTFLOW;
  }

  M_prescribed = 0.0;

  if (restart_from_hdf == TRUE) {
    restart_read();
    ZGLOOP { NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP) = -1.0; }
  }
  #if (NDIM==3)
  else if(restart_from_3d == TRUE) {
    // Read 3D data
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
  else {  // load 1-d progenitor or restart from evolved 1-d model
    if (mpi_io_proc()) {
      fp = fopen(model_name,"r");
      if (fp==NULL) {
        fprintf(stderr,"[init_problem]:  problem opening the 1-d model %s!\n",model_name);
        fflush(stderr);
        exit(1);
      }

      // Count & verify the number of lines in the 1-d model
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

    nrad1d = 2*ngroups;  // 2 = 1 (for rad energy density) + 1 (for rad 1-flux)

    // Allocate 1d model vectors
    r_model          = malloc_rank1(lines, sizeof *r_model         );
    rho_model        = malloc_rank1(lines, sizeof *rho_model       );
    v_model          = malloc_rank1(lines, sizeof *v_model         );
    ye_model         = malloc_rank1(lines, sizeof *ye_model        );
    temp_model       = malloc_rank1(lines, sizeof *temp_model      );
    s_model          = malloc_rank1(lines, sizeof *s_model         );
    rad_model        = malloc_rank1(lines*nrad1d, sizeof *rad_model);
    lapse_model      = malloc_rank1(lines, sizeof *lapse_model     );
    lapse_edge_model = malloc_rank1(lines, sizeof *lapse_edge_model);

    Fcov = malloc_rank1(NDIM*ngroups, sizeof *Fcov);
    Fred = malloc_rank1(NDIM*ngroups, sizeof *Fred);

    if (mpi_io_proc()) {
      fprintf(stderr,"reading %d lines from %s...",lines,model_name);
      fflush(stderr);
      rewind(fp);
      if (restart_from_1d) fscanf(fp,"%[^\n]%*c",&dum[0]);  // skip over the header line
      for (i=0; i<lines; i++) {
        if (restart_from_1d) {  // restart from evolved 1-d model
          // r, RHO, UU, U1, U2, U3, YE, [rad stuff], PRESS, CS, TEMP [MEV], ENT, GAMMA, MASS COORD
          fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf ",
            &r_model[i],&rho_model[i],&dummy,&v_model[i],&dummy,&dummy,&ye_model[i]);
          for (vv=0; vv<nrad1d; vv++) fscanf(fp,"%lf ",&rad_model[i*nrad1d + vv]);  // read rad vars in 1d
          fscanf(fp,"%lf %lf %lf %lf %lf %lf",&dummy,&dummy,&temp_model[i],&dummy,&dummy,&dummy);
          fscanf(fp,"%lf %lf ",&lapse_model[i],&lapse_edge_model[i]);
          fscanf(fp,"\n");
        }
        else {  // load 1-d progenitor
          // i, mass, r, RHO, TEMP [K], YE, U1, 0
          fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf\n",&idum,&dummy,&r_model[i],
            &rho_model[i],&temp_model[i],&ye_model[i],&v_model[i],&dummy);
          temp_model[i] *= K_TO_MEV;
        }
      }
      fclose(fp);
      fprintf(stderr, "done!\n");
      fflush(stderr);
    }
    #if (USE_MPI==TRUE)
    MPI_Bcast(&lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(r_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(temp_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(ye_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rad_model, lines*nrad1d, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lapse_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lapse_edge_model, lines, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    #endif

    ZGLOOP { NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP) = -1.0; }

    #pragma omp parallel for private(ii,jj,kk,r,th,R,i,del,rho,ye,temp,v,g,dd) schedule(dynamic)
    ISLOOP(ii) {
      r = r_of_x(startx[0] + (ii+0.5)*dx[0]);
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
      // printf("myrank=%d, ii=%d, i=%d, r=%e, del=%e, rho=%e, ye=%e, temp=%e, v=%e\n",myrank,ii,i,r,del,rho,ye,temp,v);
      #if (NDIM>1)
      JSLOOP(ii,jj) {
      #endif
        #if (NDIM==3)
        KSLOOP(ii,jj,kk) {
        #endif
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,RHO ) = rho;
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,UU  ) = get_model_u(rho,temp,ye);
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,U1  ) = v/ND_ELEM(geom,ii,jj,kk).scale[0][0];
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,U2  ) = 0.0;
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,U3  ) = 0.0;
          NDP_ELEM_LINEAR(sim_p,  ii,jj,kk,YE  ) = ye;
          NDP_ELEM_LINEAR(sim_eos,ii,jj,kk,TEMP) = temp;
          GLOOP {
            NDP_ELEM_LINEAR(sim_p,ii,jj,kk,irad1+g) = 0.0;
            DLOOP { NDP_ELEM_LINEAR(sim_p,ii,jj,kk,ifrad1+g*NDIM+dd) = 0.0; }
            if (restart_from_1d) {
              NDP_ELEM_LINEAR(sim_p,ii,jj,kk,irad1+g      ) = rad_model[ii*nrad1d + g];
              NDP_ELEM_LINEAR(sim_p,ii,jj,kk,ifrad1+g*NDIM) = rad_model[ii*nrad1d + ngroups + g];
            }
          }
          enforce_flux_limit(&NDP_ELEM_LINEAR(sim_p,ii,jj,kk,irad1),&NDP_ELEM_LINEAR(sim_p,ii,jj,kk,ifrad1),Fcov,Fred,&ND_ELEM(geom,ii,jj,kk),0);
          #if (GR_MONOPOLE==TRUE)
          gr_lapse[ii]      = (restart_from_1d) ? lapse_model[ii]      : 1.0;
          gr_lapse_edge[ii] = (restart_from_1d) ? lapse_edge_model[ii] : 1.0;
          #endif
          #if (DO_ROTATION==TRUE)
          /* Rotational profile is Omega(R) = Omega0/(1+(R/A)^2), where R is the cylindrical radius */
          th = th_of_x(startx[1] + (jj+0.5)*DJS(ii)*dx[1]);
          R = r*sin(th);
          NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U3) = R*rotate_Omega0/(1.0 + pow(R/(rotate_A*km),2));
          NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U3) /= ND_ELEM(geom,ii,jj,kk).scale[0][2];
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

    free(Fcov);
    free(Fred);

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
    fprintf(stderr,"u1_out  = %23.16e\n",u1_out*ND_ELEM(geom,ii,jj,kk).scale[0][0]);
    fprintf(stderr,"ye_out  = %23.16e\n",ye_out);
    fflush(stderr);
  }

  #if (GR_MONOPOLE==TRUE)
  /* copy 1-d lapses across grid */
  gr_lapse_edge[n1] = gr_lapse_edge[n1-1];
  ZLOOP {
    ND_ELEM(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    ND_ELEM(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
  }
  ii = istop[0];
  #if (NDIM>1)
  JSLOOP(ii,jj) {
  #endif
    #if (NDIM==3)
    KSLOOP(ii,jj,kk) {
    #endif
      ND_ELEM(geom,ii,jj,kk).lapse[1] = gr_lapse_edge[ii];
    #if (NDIM==3)
    }
    #endif
  #if (NDIM>1)
  }
  #endif

  #if (NDIM==2)
  ISLOOP(ii) {
    jj = JS(ii,istop[1]);
    ND_ELEM(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
  }
  #endif

  #if (NDIM==3)
  ISLOOP(ii) {
    jj = JS(ii,istop[1]);
    KSLOOP(ii,jj,kk) {
      ND_ELEM(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
    JSLOOP(ii,jj) {
      kk = KS(ii,jj,istop[2]);
      ND_ELEM(geom,ii,jj,kk).lapse[0] = gr_lapse[ii];
    }
  }
  #endif
  #endif /* GR_MONOPOLE==TRUE */

  reset_boundaries(sim_p);
  complete_mpi_communication(0);

  return;
}

void prob_bounds(int i, int j, int k, double *p)
{
  p[RHO] = rho_out;
  p[UU]  = u_out;
  p[U1]  = u1_out;
  p[YE]  = ye_out;

  return;
}

#if (DEBUG==TRUE)
const int check_ghost = 1;
const int verbose = 1;
const int check_shells = 1;
#endif
void analysis_preloop()
{
  #if (DEBUG==TRUE)
  calculate_total_mass(verbose);

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

  if (!perturbed && !restart_from_hdf && !restart_from_1d && !restart_from_3d) {
    if (tbounce>0.0 && perturb_delay>=0.0 && t>=tbounce+perturb_delay) {
      #if (NDIM>1 && PERTURB==DENSITY_RANDOM)
      /* Implements the perturbation scheme of Summa et al. (2015),
       * which perturbs the density once by some relative amount */
      if (mpi_io_proc()) fprintf(stderr,"[init_problem]:  Perturbing density...\n");
      rseed = -(myrank+1);
      ZLOOP { NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO) *= 1.0 + (perturb_level/100.0)*(2.0*ran2(&rseed) - 1.0); }
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
}

#if (PERTURB==VELOCITY_SPHERICAL_HARMONIC)
double psi(int l, int m, int n, double r, double th, double ph, double rmin, double rmax)
{
  // Note that constant coefficients in the Ylm's are unimportant--they get normalized out anyway
  return (sqrt(sin(th))/r)*sin(n*M_PI*(r-rmin)/(rmax-rmin))*gsl_sf_legendre_sphPlm(l,m,cos(th))*cos(m*ph);
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
    r   = r_of_x(x[0]);
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
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1) += C*ND_ELEM(dvr,ii,jj,kk)/ND_ELEM(geom,ii,jj,kk).scale[0][0];
    NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U2) += C*ND_ELEM(dvt,ii,jj,kk)/ND_ELEM(geom,ii,jj,kk).scale[0][1];
  }
}
#endif /* PERTURB==VELOCITY_SPHERICAL_HARMONIC */
