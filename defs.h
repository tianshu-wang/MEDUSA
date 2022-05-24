#include "build.h"

r_of_x_struct_t rx_info;
th_of_x_struct_t thx_info;

#if (USE_LINEAR_ALLOCATION==TRUE)
double ** sim_p;		// primitive variables at integer steps
double ** sim_ph;		// primitive variables at half-steps
double ** sim_eos;
double ** sim_src;		// other source terms (e.g. radiation, etc.)
int     * sim_shock_flag;
double  * sim_Phi;
zone_geom * geom;
#if (USE_LARGER_STEP==TRUE)
double ** sim_dudt;		// primitive variables at integer steps
double ** sim_deosdt;
#endif
//double NDP_PTR sim_f[NDIM];		// flux of conserved variables
double ** sim_fdir0;		// flux of conserved variables
double ** sim_fdir1;		// flux of conserved variables
double ** sim_fdir2;		// flux of conserved variables
double NDP_PTR sim_f1;		// flux of conserved variables at outer 1-edge
//double ND_PTR  sim_grav;		// gravitational acceleration at zone edges
//double NDP_PTR sim_vedge[NDIM];
double ** sim_vedgedir0;
double ** sim_vedgedir1;
double ** sim_vedgedir2;
#else
double NDP_PTR sim_p;		// primitive variables at integer steps
double NDP_PTR sim_ph;		// primitive variables at half-steps
double NDP_PTR sim_eos;
double NDP_PTR sim_src;		// other source terms (e.g. radiation, etc.)
int ND_PTR     sim_shock_flag;
double ND_PTR  sim_Phi;
zone_geom ND_PTR geom;
#if (USE_LARGER_STEP==TRUE)
double NDP_PTR sim_dudt;		// primitive variables at integer steps
double NDP_PTR sim_deosdt;
#endif
//double NDP_PTR sim_f[NDIM];		// flux of conserved variables
double NDP_PTR sim_fdir0;		// flux of conserved variables
double NDP_PTR sim_fdir1;		// flux of conserved variables
double NDP_PTR sim_fdir2;		// flux of conserved variables
double NDP_PTR sim_f1;		// flux of conserved variables at outer 1-edge
//double ND_PTR  sim_grav;		// gravitational acceleration at zone edges
//double NDP_PTR sim_vedge[NDIM];
double NDP_PTR sim_vedgedir0;
double NDP_PTR sim_vedgedir1;
double NDP_PTR sim_vedgedir2;
#endif
double NDP_PTR sim_vedge1;  // velocity at outer 1-edge
//double ND_PTR  sim_pedge[NDIM];
//grid_data sim;
tracer_t *tracers=NULL;

//tianshu
double kom_dtmin,rhocut,kom_epsilon,kom_delta;
double Ecut1,Ecut2,Ecut3;
int use_kom=0;
int Ncurrent,Nskip;
int cell_count_send;
int cell_count_recv;
int cell_count;
int cell_count_all;
int ND_PTR ijk_to_I;
int NDP_PTR face_ijk_to_I;
int face_count_0,face_count_1,face_count_2;
int * I_to_ijk;

// in rad_fluid.c
double *Trad;

// used in reconstruction
//double **pencil;
//double **pleft;
//double **pright;
//double *flatten;
//double pencil[NSIZE][NSIZE_GRID];
//double pleft [NSIZE][NSIZE_GRID];
//double pright[NSIZE][NSIZE_GRID];
//double pencil[NSIZE][2*NG];
//double pleft [NSIZE][2*NG];
//double pright[NSIZE][2*NG];
//double flatten[NSIZE_GRID];
double *pcenter;
double *dq, *d2q;
int *interp_order;
int ninterp,hydro_interp_order,rad_interp_order;
double *interp_vol[3];
double *mu0,*nu0,*alpha0,*beta0,*Gamma0;
double *mu1,*nu1,*alpha1,*beta1,*Gamma1;
double *mu2,*nu2,*alpha2,*beta2,*Gamma2;
double **alpha1s,**beta1s,**Gamma1s;
double **alpha2s,**beta2s,**Gamma2s;


// keep track of space and time
int istep, max_steps;
double max_wtime, tstart, tcurr, telapsed;
double t,dt,dtmin,min_dt,dtmax,tmax,initial_dt;
double dx[SPACEDIM],startx[SPACEDIM];
double tbounce,rho_cent_max;
int dt_flag;
int nstep_log;
int nstep_analysis;

// tracer related
unsigned int n_tracer_target,n_tracer_current,n_tracer_global,consistent_tracer_ids=0;
double mass_inside_tracers, mass_outside_tracers;

// minimum radial resolution for CCSN simulations
double ccsn_dr_min;

// grid index coordinates
int n1,n2,n3;	// problem size
int global_grid_dims[SPACEDIM];
int istart[SPACEDIM], istop[SPACEDIM], my_grid_dims[SPACEDIM];
int max_grid_dim;
int periodic[SPACEDIM];
int *icount;
int* nj=NULL;
int* dj=NULL;
int** nk=NULL;
int** dk=NULL;
int max_jrefine_level;
int istop_jrefine;
int max_krefine_level;
int istop_krefine;
int* jstop_krefine;
int* jstart_kcoarsen;

// radiation related
int ngroups, nr1, nr2, nr3;
int irad1, irad2, irad3;
int ifrad1, ifrad2, ifrad3;
double emin1,emin2,emin3;
double emax1,emax2,emax3;
double *egroup, *egroup_erg, *degroup, *xi;
double *spec_factor;
double implicit_err_tol;
double du_tol, dye_tol, dT_tol;
double chat,vmax,taumax,chat_safety;
int use_chat;
int max_implicit_iter;
int implicit_iter;
int force_flux_limit;

// gravity related
double total_mass;
int multipole_lmax;
double *gr_grav;
double *gr_lapse, *gr_lapse_edge;

// i/o related stuff
double t_next_dump, t_next_restart, t_next_pdump;
int i_next_restart, di_restart;
double dt_dump, dt_restart, dt_pdump;
char model_name[256];
char freq_type[4];
int restart_create_1d;
int restart_from_1d;
int restart_from_3d;
int restart_from_hdf;
int restart_from_last;
double restart_delay;
char restart_file[256];
char restart_dir[256];
char grid_units[10];
int dump_hdf;
int dump_rad_vars;
int detect_tbounce;

// mpi stuff
int numprocs, myrank, iorank;
#if(USE_MPI==TRUE)
proc_info *proc;
MPI_Comm mpi_interface_comm,mpi_active_comm;
int neighbor[NDIM][2];
MPI_Datatype ghost_cells[NDIM];
MPI_Datatype *interface_cells;
int nproc_2, nproc_3, nprocs_spider;
int cells_per_n2, cells_per_n3;
int *interface_ranks;
int *interface_counts;
int *interface_displacements;
int num_requests;
MPI_Request *mpi_requests;
MPI_Status *mpi_status;
#endif


// boundary condition flags
bc_type *bc;
int half_step_sync;

// coordinate/eos/problem specific
double eostable[16*50*300*300]; 
//double rtrans,rsparse_fact,rsparse_trans;
double r_full_res;
double M_prescribed,Rschw;
double g0[SPACEDIM];
double gam,Kpoly;
double rho_floor;
double e_floor;
double temp_guess;
double *rcenter;
double *redge;
double *dr;
char opac_param_file[1024];
char opac_file[1024];
char eos_file[1024];
char inelastic_root[1024];
double L_enu;
double T_enu;
int include_inelastic;

Array * qpole_time = NULL;
Array * qpole_mom = NULL;

// velocity perturbations
double perturb_r1m, perturb_r1p;
int perturb_l1, perturb_m1, perturb_n1;
double perturb_dv1;
double perturb_r2m, perturb_r2p;
int perturb_l2, perturb_m2, perturb_n2;
double perturb_dv2;
double perturb_r3m, perturb_r3p;
int perturb_l3, perturb_m3, perturb_n3;
double perturb_dv3;

// density perturbation
double perturb_level;
double perturb_delay;

// rotation
double rotate_Omega0,rotate_A;

// outer radius, only used by init_*_ccsn.c
double outer_radius;

int decomp_only,decomp_link,decomp_npmin,decomp_npmax,decomp_npskip,decomp_from_file;
char decomp_path[256];

double cfl;
int nvars,ncomp,nhydro;
double dt_init_safety;
double dt_grow_safety;

int step_part;

Timer * timers[MAX_NTIMERS] = {NULL};
int ntimers = 0;
int nstep_timing = 100;
