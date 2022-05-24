#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <float.h>
#include <inttypes.h>
#include <time.h>

#include "build.h"
#include "array.h"
#include "timer.h"

#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wabsolute-value"
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#ifndef M_PI
#define M_PI    3.14159265358979323846264338327
#endif
#ifndef M_PI_2
#define M_PI_2  1.57079632679489661923132169164
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242
#endif
#define ONE_3RD 0.33333333333333333333333333333

#define NSIZE 20 // must be greater than ninterp
#define NSIZE_GRID 1280 // must be greater than max_grid_dim

#define FALSE	0
#define TRUE	1

#if (USE_OMP==TRUE)
#include <omp.h>
#endif
#if (USE_MPI==TRUE)
#include <mpi.h>
#endif
#if (USE_GPU==TRUE)
#define GPU_PRAGMA(x) _Pragma(#x)
#else
#define GPU_PRAGMA(x) 
#endif

#define DUMP_RECON	FALSE

#define NPRIM		  5  // number of primitive variables; is always 5 (rho, u, v1, v2, v3)
#define NEOS		  5  // number of eos-related variables (in the order: pressure, cs, temperature, entropy, Gamma...)

#define SPACEDIM 	3  // number of spatial dimensions; is always 3, unless we find ourselves in some other universe
#define NG        2  // number of ghost zones

// Generic boundary condition flags:
#define OUTFLOW			      0  // boundary values copied for last cell
#define REFLECT		      	1  // boundary values reflected over boundary
#define EXTRAP			      2  // linear extrapolation
#define PERIODIC		      3
#define PROB			        4  // boundary values specified in init_*.c file
#define SPHERICAL_ORIGIN	5  // (should this actually differ from reflect?)
#define CCSN_OUTER        6
#define RADEXTRAP         7
#define DISK              8

// convenience definitions.  never change.
#if (NDIM==1)
#define ND_PTR	*
#define NDP_PTR **
#define NPDIM	2
#elif (NDIM==2)
#define ND_PTR	**
#define NDP_PTR	***
#define NPDIM	3
#else
#define ND_PTR 	***
#define NDP_PTR	****
#define NPDIM	4
#endif

#define ERROR_NONE            0
#define ERROR_INVALID         1
#define ERROR_INDEX           2
#define ERROR_NOTFOUND        3
#define ERROR_IO              4
#define ERROR_NOTIMPL         5
#define ERROR_OUTOFMEM        6
#define ERROR_PARFILE         7
#define ERROR_LOGICAL         8
#define ERROR_INTERNAL        9

#define CHKERRQ(IERR) \
    if((IERR)) { \
        fprintf(stdout, "%s:%d", __FILE__, __LINE__); \
        return IERR; \
    }
#define THROW_ERROR(IERR, MSG) \
    fprintf(stdout, "%s:%d %s", __FILE__, __LINE__, MSG); \
    return (IERR)
#define RESET(T,X,Y)          (*((T *)(&(X)))) = (Y)

/*
typedef struct {
        #if (USE_LINEAR_ALLOCATION==TRUE)
	double ** p;		// primitive variables at integer steps
	double ** ph;		// primitive variables at half-steps
	double ** eos;
        #else
	double NDP_PTR p;		// primitive variables at integer steps
	double NDP_PTR ph;		// primitive variables at half-steps
	double NDP_PTR eos;
        #endif
	double NDP_PTR f[NDIM];		// flux of conserved variables
	double NDP_PTR f1;		// flux of conserved variables at outer 1-edge
	double ND_PTR grav;		// gravitational acceleration at zone edges
	double NDP_PTR src;		// other source terms (e.g. radiation, etc.)
	double NDP_PTR vedge[NDIM];
	double NDP_PTR vedge1;  // velocity at outer 1-edge
	double ND_PTR  pedge[NDIM];
	int ND_PTR  shock_flag;
        double ND_PTR Phi;
#if (USE_LARGER_STEP==TRUE)
	double NDP_PTR dudt;		// primitive variables at integer steps
	double NDP_PTR deosdt;
#endif
} grid_data;
*/

/******************************
geometry structure:
  0 -> zone center
  1 -> i-1/2
  2 -> j-1/2
  3 -> k-1/2
*******************************/
typedef struct {
  double area[NDIM];
  double conn[SPACEDIM][SPACEDIM][SPACEDIM];
  double gcov[1+NDIM][SPACEDIM];
  double gcon[SPACEDIM];
  double scale[1+NDIM][SPACEDIM];
  double volume;
  double lapse[2];    // 0->center 1->edge
	//double vpart[NDIM];
} zone_geom;

typedef struct _tracer {
  char active;
  int id;
  double x[NDIM];
  double xh[NDIM];
  double mass;
  //double *hist;
  struct _tracer *next;
} tracer_t;

typedef struct _tracer_mpi_packet {
    double d[2*NDIM+1];
    int id;
    //char active;
} tracer_mpi_packet_t;

typedef struct {
   double rtrans;
   double facta;
   double factb;
   double factc;
   double rsparse_fact;
   double rsparse_trans;
} r_of_x_struct_t;

typedef struct {
   double poly_norm;
   double poly_xt;
   double poly_alpha;
   double nice_norm;
   double nice_alpha;
} th_of_x_struct_t;

/******************************
       global variables
 ******************************/
extern r_of_x_struct_t rx_info;
extern th_of_x_struct_t thx_info;
#if (USE_LINEAR_ALLOCATION==TRUE)
extern double ** sim_p;		// primitive variables at integer steps
extern double ** sim_ph;		// primitive variables at half-steps
extern double ** sim_eos;
extern double ** sim_src;		// other source terms (e.g. radiation, etc.)
extern int     * sim_shock_flag;
extern double  * sim_Phi;
extern zone_geom * geom;
#if (USE_LARGER_STEP==TRUE)
extern double ** sim_dudt;		// primitive variables at integer steps
extern double ** sim_deosdt;
#endif
//extern double NDP_PTR sim_f[NDIM];		// flux of conserved variables
extern double ** sim_fdir0;		// flux of conserved variables
extern double ** sim_fdir1;		// flux of conserved variables
extern double ** sim_fdir2;		// flux of conserved variables
//extern double ND_PTR  sim_grav;		// gravitational acceleration at zone edges
//extern double NDP_PTR sim_vedge[NDIM];
extern double ** sim_vedgedir0;
extern double ** sim_vedgedir1;
extern double ** sim_vedgedir2;
#else
extern double NDP_PTR sim_p;		// primitive variables at integer steps
extern double NDP_PTR sim_ph;		// primitive variables at half-steps
extern double NDP_PTR sim_eos;
extern double NDP_PTR sim_src;		// other source terms (e.g. radiation, etc.)
extern int ND_PTR     sim_shock_flag;
extern double ND_PTR  sim_Phi;
extern zone_geom ND_PTR geom;
#if (USE_LARGER_STEP==TRUE)
extern double NDP_PTR sim_dudt;		// primitive variables at integer steps
extern double NDP_PTR sim_deosdt;
#endif
//extern double NDP_PTR sim_f[NDIM];		// flux of conserved variables
extern double NDP_PTR sim_fdir0;		// flux of conserved variables
extern double NDP_PTR sim_fdir1;		// flux of conserved variables
extern double NDP_PTR sim_fdir2;		// flux of conserved variables
//extern double ND_PTR  sim_grav;		// gravitational acceleration at zone edges
//extern double NDP_PTR sim_vedge[NDIM];
extern double NDP_PTR sim_vedgedir0;
extern double NDP_PTR sim_vedgedir1;
extern double NDP_PTR sim_vedgedir2;
#endif
//extern double NDP_PTR sim_f1;		// flux of conserved variables at outer 1-edge
//extern double NDP_PTR sim_vedge1;  // velocity at outer 1-edge
//extern double ND_PTR  sim_pedge[NDIM];
//extern grid_data sim;
extern tracer_t *tracers;

//tianshu
extern double kom_dtmin,rhocut;
extern double Ecut1,Ecut2,Ecut3;
extern int use_kom;
extern double kom_epsilon,kom_delta;
extern int Ncurrent,Nskip;
extern int cell_count_send;
extern int cell_count_recv;
extern int cell_count;
extern int cell_count_all;
extern int ND_PTR ijk_to_I;
extern int NDP_PTR face_ijk_to_I;
extern int face_count_0,face_count_1,face_count_2;
//GPU_PRAGMA(omp declare target)
extern int * I_to_ijk;
//GPU_PRAGMA(omp end declare target)


// in rad_fluid.c
extern double *Trad;

// used in reconstruction
//GPU_PRAGMA(omp declare target)
//extern double **pencil;
//extern double **pleft;
//extern double **pright;
/******************************
       global variables
 ******************************/
extern r_of_x_struct_t rx_info;
#if (USE_LINEAR_ALLOCATION==TRUE)
extern double ** sim_p;		// primitive variables at integer steps
extern double ** sim_ph;		// primitive variables at half-steps
extern double ** sim_eos;
extern double ** sim_src;		// other source terms (e.g. radiation, etc.)
extern int     * sim_shock_flag;
extern double  * sim_Phi;
extern zone_geom * geom;
#if (USE_LARGER_STEP==TRUE)
extern double ** sim_dudt;		// primitive variables at integer steps
extern double ** sim_deosdt;
#endif
//extern double NDP_PTR sim_f[NDIM];		// flux of conserved variables
extern double ** sim_fdir0;		// flux of conserved variables
extern double ** sim_fdir1;		// flux of conserved variables
extern double ** sim_fdir2;		// flux of conserved variables
//extern double ND_PTR  sim_grav;		// gravitational acceleration at zone edges
//extern double NDP_PTR sim_vedge[NDIM];
extern double ** sim_vedgedir0;
extern double ** sim_vedgedir1;
extern double ** sim_vedgedir2;
#else
extern double NDP_PTR sim_p;		// primitive variables at integer steps
extern double NDP_PTR sim_ph;		// primitive variables at half-steps
extern double NDP_PTR sim_eos;
extern double NDP_PTR sim_src;		// other source terms (e.g. radiation, etc.)
extern int ND_PTR     sim_shock_flag;
extern double ND_PTR  sim_Phi;
extern zone_geom ND_PTR geom;
#if (USE_LARGER_STEP==TRUE)
extern double NDP_PTR sim_dudt;		// primitive variables at integer steps
extern double NDP_PTR sim_deosdt;
#endif
//extern double NDP_PTR sim_f[NDIM];		// flux of conserved variables
extern double NDP_PTR sim_fdir0;		// flux of conserved variables
extern double NDP_PTR sim_fdir1;		// flux of conserved variables
extern double NDP_PTR sim_fdir2;		// flux of conserved variables
//extern double ND_PTR  sim_grav;		// gravitational acceleration at zone edges
//extern double NDP_PTR sim_vedge[NDIM];
extern double NDP_PTR sim_vedgedir0;
extern double NDP_PTR sim_vedgedir1;
extern double NDP_PTR sim_vedgedir2;
#endif
//extern double NDP_PTR sim_f1;		// flux of conserved variables at outer 1-edge
//extern double NDP_PTR sim_vedge1;  // velocity at outer 1-edge
//extern double ND_PTR  sim_pedge[NDIM];
//extern grid_data sim;
extern tracer_t *tracers;

//tianshu
extern double kom_dtmin,rhocut;
extern double Ecut1,Ecut2,Ecut3;
extern int use_kom;
extern double kom_epsilon,kom_delta;
extern int Ncurrent,Nskip;
extern int cell_count_send;
extern int cell_count_recv;
extern int cell_count;
extern int cell_count_all;
extern int ND_PTR ijk_to_I;
extern int NDP_PTR face_ijk_to_I;
extern int face_count_0,face_count_1,face_count_2;
//GPU_PRAGMA(omp declare target)
extern int * I_to_ijk;
//GPU_PRAGMA(omp end declare target)


// in rad_fluid.c
extern double *Trad;

// used in reconstruction
//GPU_PRAGMA(omp declare target)
//extern double **pencil;
//extern double **pleft;
//extern double **pright;
//extern double *flatten;
//extern double pencil[NSIZE][NSIZE_GRID];
//extern double pleft [NSIZE][NSIZE_GRID];
//extern double pright[NSIZE][NSIZE_GRID];
//extern double pencil[NSIZE][2*NG];
//extern double pleft [NSIZE][2*NG];
//extern double pright[NSIZE][2*NG];
//extern double flatten[NSIZE_GRID];
extern double *pcenter;
extern double *dq,*d2q;
extern int *interp_order;
extern int ninterp,hydro_interp_order,rad_interp_order;
extern double *interp_vol[3];
extern double *mu0,*nu0,*alpha0,*beta0,*Gamma0;
extern double *mu1,*nu1,*alpha1,*beta1,*Gamma1;
extern double *mu2,*nu2,*alpha2,*beta2,*Gamma2;
extern double **alpha1s,**beta1s,**Gamma1s;
extern double **alpha2s,**beta2s,**Gamma2s;

// keep track of space and time
extern int istep, max_steps;
extern double max_wtime, tstart, tcurr, telapsed;
extern double t,dt,dtmin,min_dt,dtmax,tmax,initial_dt;
extern double dx[SPACEDIM],startx[SPACEDIM];
extern double tbounce,rho_cent_max;
extern int dt_flag;
extern int nstep_log;
extern int nstep_analysis;
//GPU_PRAGMA(omp end declare target)

// tracer related
extern unsigned int n_tracer_target,n_tracer_current,n_tracer_global,consistent_tracer_ids;
extern double mass_inside_tracers, mass_outside_tracers;

// minimum radial resolution for CCSN simulations
extern double ccsn_dr_min;

// grid index coordinates
extern int n1,n2,n3;	// problem size
extern int global_grid_dims[SPACEDIM];
extern int istart[SPACEDIM], istop[SPACEDIM], my_grid_dims[SPACEDIM];
extern int max_grid_dim;
extern int periodic[SPACEDIM];
extern int *icount;
extern int* nj;
extern int* dj;
extern int** nk;
extern int** dk;
extern int max_jrefine_level;
extern int istop_jrefine;
extern int max_krefine_level;
extern int istop_krefine;
extern int* jstop_krefine;
extern int* jstart_kcoarsen;

// radiation related
extern int ngroups, nr1, nr2, nr3;
extern int irad1, irad2, irad3;
extern int ifrad1, ifrad2, ifrad3;
extern double emin1,emin2,emin3;
extern double emax1,emax2,emax3;
extern double *egroup, *egroup_erg, *degroup, *xi;
extern double implicit_err_tol;
extern double du_tol, dye_tol, dT_tol;
extern double chat,vmax,taumax,chat_safety;
extern int use_chat;
extern int max_implicit_iter;
extern int implicit_iter;
extern int force_flux_limit;
extern double *spec_factor;

// gravity related
extern double total_mass;
extern int multipole_lmax;
extern double *gr_grav;
extern double *gr_lapse, *gr_lapse_edge;

// i/o related stuff
extern double t_next_dump, t_next_restart, t_next_pdump;
extern int i_next_restart, di_restart;
extern double dt_dump, dt_restart, dt_pdump;
extern char model_name[256];
extern char freq_type[4];
extern int restart_create_1d;
extern int restart_from_1d;
extern int restart_from_3d;
extern int restart_from_hdf;
extern int restart_from_last;
extern double restart_delay;
extern char restart_file[256];
extern char restart_dir[256];
extern char grid_units[10];
extern int dump_hdf;
extern int dump_rad_vars;
extern int detect_tbounce;

// mpi stuff
extern int numprocs, myrank, iorank;
#if (USE_MPI==TRUE)
typedef struct proc_info_ {
	int istart[NDIM],istop[NDIM];
	uint64_t ncells;
	uint32_t ix,iy,iz;
	uint32_t nneigh;
	uint32_t *neigh_id;
	int *neigh_dir;
} proc_info;

extern proc_info *proc;
extern MPI_Comm mpi_interface_comm,mpi_active_comm;
extern int neighbor[NDIM][2];
extern MPI_Datatype ghost_cells[NDIM];
extern  MPI_Datatype *interface_cells;
extern int nproc_2, nproc_3, nprocs_spider;
extern int cells_per_n2, cells_per_n3;
extern int *interface_ranks;
extern int *interface_counts;
extern int *interface_displacements;
extern int num_requests;
extern MPI_Request *mpi_requests;
extern MPI_Status *mpi_status;
#endif

// boundary condition flags
typedef struct {
	int lo[NDIM];
	int hi[NDIM];
} bc_type;
extern bc_type *bc;
extern int half_step_sync;

// coordinate/eos/problem specific
extern double eostable[16*50*300*300]; 
extern double rtrans,rsparse_fact,rsparse_trans,r_full_res;
extern double M_prescribed,Rschw;
extern double g0[SPACEDIM];
extern double gam,Kpoly;
//GPU_PRAGMA(omp declare target)
extern double rho_floor;
extern double e_floor;
//GPU_PRAGMA(omp end declare target)
extern double temp_guess;
extern double *rcenter;
extern double *redge;
extern double *dr;
extern char opac_param_file[1024];
extern char opac_file[1024];
extern char eos_file[1024];
extern char inelastic_root[1024];
extern double L_enu;
extern double T_enu;
extern int include_inelastic;

// quadrupole data
extern Array * qpole_time;
extern Array * qpole_mom;

// velocity perturbations
extern double perturb_r1m, perturb_r1p;
extern int perturb_l1, perturb_m1, perturb_n1;
extern double perturb_dv1;
extern double perturb_r2m, perturb_r2p;
extern int perturb_l2, perturb_m2, perturb_n2;
extern double perturb_dv2;
extern double perturb_r3m, perturb_r3p;
extern int perturb_l3, perturb_m3, perturb_n3;
extern double perturb_dv3;

// density perturbation
extern double perturb_level;
extern double perturb_delay;

// rotation
extern double rotate_Omega0,rotate_A;

// outer radius, only used by init_*_ccsn.c
extern double outer_radius;

extern int decomp_only,decomp_link,decomp_npmin,decomp_npmax,decomp_npskip,decomp_from_file;
extern char decomp_path[256];

extern double cfl;
extern int nvars,ncomp,nhydro;
extern double dt_init_safety;
extern double dt_grow_safety;

extern int step_part;

// timers
#define MAX_NTIMERS        (512)
extern Timer * timers[MAX_NTIMERS];
extern int ntimers;
extern int nstep_timing;

/******************************
      macro definitions
 ******************************/
#define GET_IJK_FROM_I(II,ii,jj,kk) ii=I_to_ijk[3*II];jj=I_to_ijk[3*II+1];kk=I_to_ijk[3*II+2];

#if (USE_LINEAR_ALLOCATION==TRUE)
#define ZLOOP_LINEAR for(int II=0;II<cell_count;II++)
#define ZGLOOP_LINEAR for(int II=0;II<cell_count_all;II++)
#define ND_ELEM_LINEAR(arr,i,j,k) arr[ND_ELEM(ijk_to_I,i,j,k)]
//#define ND_ELEM_LINEAR_F(arr,dir,i,j,k) arr[ND_ELEM(face_ijk_to_I[dir],i,j,k)]
#define NDP_ELEM_LINEAR(arr,i,j,k,l) arr[ND_ELEM(ijk_to_I,i,j,k)][l]
#define NDP_ELEM_LINEAR_F(arr,dir,i,j,k,l) arr[l][ND_ELEM(face_ijk_to_I[dir],i,j,k)]
//#define NDP_ELEM_LINEAR_F(arr,dir,i,j,k,l) arr[ND_ELEM(face_ijk_to_I[dir],i,j,k)][l]
#define ND_ELEM_LINEAR_REORDERED(arr,i,j,k) arr[II]
#define NDP_ELEM_LINEAR_REORDERED(arr,i,j,k,l) arr[II][l]
#else
#define ZLOOP_LINEAR ZLOOP
#define ZGLOOP_LINEAR ZGLOOP
#define ND_ELEM_LINEAR(arr,i,j,k)            ND_ELEM(arr,i,j,k)
#define ND_ELEM_LINEAR_F(arr,dir,i,j,k)      ND_ELEM(arr,i,j,k)
#define NDP_ELEM_LINEAR(arr,i,j,k,l)         NDP_ELEM(arr,i,j,k,l)
#define NDP_ELEM_LINEAR_F(arr,dir,i,j,k,l)     NDP_ELEM(arr,i,j,k,l)
#define ND_ELEM_LINEAR_REORDERED(arr,i,j,k)    ND_ELEM(arr,i,j,k)
#define NDP_ELEM_LINEAR_REORDERED(arr,i,j,k,l) NDP_ELEM(arr,i,j,k,l)
#endif

#define ISLOOP(i)      for (i=istart[0]; i<istop[0]; i++)
#define ISGLOOP(i)     for (i=istart[0]-NG; i<istop[0]+NG; i++)
#if ((GEOM==SPHERICAL || GEOM==CYLINDRICAL) && NDIM>1)
#define DJS(i)         (dj[i])
#define JS(i,j)        ((j)/DJS(i))
#define JSLOOP(i,j)    for (j=JS(i,istart[1]); j<JS(i,istop[1]); j++)
#define JSGLOOP(i,j)   for (j=JS(i,istart[1])-NG; j<JS(i,istop[1])+NG; j++)
#else
#define DJS(i)         (1)
#define JS(i,j)        (j)
#define JSLOOP(i,j)    for (j=istart[1]; j<istop[1]; j++)
#define JSGLOOP(i,j)   for (j=istart[1]-NG; j<istop[1]+NG; j++)
#endif
#if (GEOM==SPHERICAL && NDIM==3)
#define DKS(i,j)       (dk[i][j])
#define KS(i,j,k)      ((k)/DKS(i,j))
#define KSLOOP(i,j,k)  for (k=KS(i,j,istart[2]); k<KS(i,j,istop[2]); k++)
#define KSGLOOP(i,j,k) for (k=KS(i,j,istart[2])-NG; k<KS(i,j,istop[2])+NG; k++)
#else
#define DKS(i,j)       (1)
#define KS(i,j,k)      (k)
#define KSLOOP(i,j,k)  for (k=istart[2]; k<istop[2]; k++)
#define KSGLOOP(i,j,k) for (k=istart[2]-NG; k<istop[2]+NG; k++)
#endif

#if (NDIM==1)
#define ND_ELEM(arr,i,j,k)	     arr[i]
#define NDP_ELEM(arr,i,j,k,l)	   arr[i][l]
#define ND_ARRAY_LIST(arr)	     arr[0]
#define ND_ARRAY_SET(arr,x,y,z)	 arr[0]=(x)
#define DOT(arr1,arr2)		       arr1[0]*arr2[0]
#define ZLOOP	                   ISLOOP(ii)
#define ZGLOOP                   ISGLOOP(ii)
#elif (NDIM==2)
#define ND_ELEM(arr,i,j,k)	     arr[i][j]
#define NDP_ELEM(arr,i,j,k,l)	   arr[i][j][l]
#define ND_ARRAY_LIST(arr)	     arr[0],arr[1]
#define ND_ARRAY_SET(arr,x,y,z)	 arr[0]=(x);arr[1]=(y)
#define DOT(arr1,arr2)		       (arr1[0]*arr2[0] + arr1[1]*arr2[1])
#define ZLOOP                    ISLOOP(ii) JSLOOP(ii,jj)
#define ZGLOOP                   ISGLOOP(ii) JSGLOOP(ii,jj)
#else  /* NDIM==3 */
#define ND_ELEM(arr,i,j,k)	     arr[i][j][k]
#define NDP_ELEM(arr,i,j,k,l)    arr[i][j][k][l]
#define ND_ARRAY_LIST(arr)	     arr[0],arr[1],arr[2]
#define ND_ARRAY_SET(arr,x,y,z)	 arr[0]=(x); arr[1]=(y); arr[2]=(z)
#define DOT(arr1,arr2)		       (arr1[0]*arr2[0] + arr1[1]*arr2[1] + arr1[2]*arr2[2])
#define ZLOOP                    ISLOOP(ii) JSLOOP(ii,jj) KSLOOP(ii,jj,kk)
#define ZGLOOP                   ISGLOOP(ii) JSGLOOP(ii,jj) KSGLOOP(ii,jj,kk)
#endif

#define VLOOP    for (vv=0; vv<nvars; vv++)
#define PLOOP    for (pp=0; pp<NPRIM; pp++)
#define DLOOP    for (dd=0; dd<NDIM; dd++)
#define SLOOP    for (dd=0; dd<SPACEDIM; dd++)
#define EOSLOOP  for (vv=0; vv<NEOS; vv++)
#define ILOOP    for (vv=0; vv<ninterp; vv++)
#define HLOOP    for (vv=0; vv<nhydro; vv++)
#define GLOOP    for (g=0; g<ngroups; g++)

#define MIN(a,b)	  ( (a) < (b) ? (a) : (b) )
#define MAX(a,b)	  ( (a) > (b) ? (a) : (b) )
#define SIGN(a)		  ( (a) < 0.0 ? -1.0 : 1.0 )
#define MINMOD(a,b) ( 0.5*(SIGN(a)+SIGN(b))*MIN(fabs(a),fabs(b)) )
#define DELTA(a,b)	( (a) == (b) ? 1 : 0 )
#define SQR(x)      ( (x)*(x) )
#define CUBE(x)     ( (x)*(x)*(x) )
#define FOURTH(x)   ( (x)*(x)*(x)*(x) )

//compiler choices
#define INTEL  0
#define GCC    1

//primitive variable mnemonics
#define RHO	   0
#define UU	   1
#define U1     2
#define U2	   3
#define U3     4
#define YE	   5
#define ESTAR	 6

//conserved variable mnemonics
#define URHO	 0
#define ETOT	 1
#define P1	   2
#define P2	   3
#define P3  	 4

//eos variables
#define PRESS	 0
#define CS	   1
#define TEMP	 2
#define ENT	   3
#define GAMMA	 4

//grav directions
#define GX1	   0
#define GX2	   1
#define GX3  	 2

//geometry definitions
#define CARTESIAN	  0
#define SPHERICAL	  1
#define CYLINDRICAL   2

// eos definitions
#define GAMMA_LAW   0
#define POLYTROPIC  1
#define COLLAPSE    2

//gravity definitions
#define NO_GRAV                   0
#define FIXED_GRAV                1
#define PRESCRIBED_GRAV           2
#define SPHERICAL_MONOPOLE_GRAV	  3
#define SPHERICAL_MULTIPOLE_GRAV  4
#define USER_GRAV                 5

// perturbation definitions
#define NONE                         0
#define VELOCITY_SPHERICAL_HARMONIC  1
#define DENSITY_RANDOM               2

// io definitions
#define SINGLE 0
#define DOUBLE 1

// Fortran-C interface macros
#define GLUE_HELPER(x,y)		x##y
#define GLUE(x,y)			GLUE_HELPER(x,y)
#define ND_FUNC(func_name,number)	GLUE(func_name,number)
#define ND_GRID_MALLOC			GLUE(grid_malloc_rank,NPDIM)

#if USE_FORTRAN_CODE
#if defined(__INTEL_COMPILER)
#define FORT_EOS_CALL(name)		eos_module_mp_##name##_
#define FORT_OPAC_CALL(name)		opacity_table_module_mp_##name##_
#elif defined(__GNUC__)
#define FORT_EOS_CALL(name)		__eos_module_MOD_##name
#define FORT_OPAC_CALL(name)		__opacity_table_module_MOD_##name
//#define FORT_EOS_CALL(name)	name
//#define FORT_OPAC_CALL(name) name
#else
#error [decs.h]:  UNKNOWN COMPILER!
#endif // compiler
#endif // USE_FORTRAN_CODE

#define FORT_CALL(name)			name##_

#if DO_TIMING==TRUE
#if USE_OMP==TRUE
#define TIMER_START(name) \
    static int timer_idx = -1; \
    Timer * timer = NULL; \
    if(0 == omp_get_thread_num()) { \
      if(timer_idx < 0) { \
        timer_idx = new_timer(name); \
      } \
      timer = get_timer(timer_idx); \
      timer_start(timer); \
    }
#define TIMER_STOP \
  if(0 == omp_get_thread_num()) { \
    timer_stop(timer); \
  }
#else // USE_OMP==FALSE
#define TIMER_START(name) \
    static int timer_idx = -1; \
    if(timer_idx < 0) { \
      timer_idx = new_timer(name); \
    } \
    Timer * timer = get_timer(timer_idx); \
    timer_start(timer)
#define TIMER_STOP      timer_stop(timer)
#endif // DO_OMP
#else // DO_TIMING==FALSE
#define TIMER_START(name)
#define TIMER_STOP
#endif // DO_TIMING

/*******************************
    function declarations
 *******************************/


/* step.c */
void update_cpu();
void step();
void estimate_dt();
void set_tbounce();

/* rad_fluid.c */
#if (USE_LINEAR_ALLOCATION==TRUE)
void get_geom_src_terms(double ** p);
void p_to_phys_src_start(double ** p);
void p_to_phys_src_finish(double ** p, double dtpush);
void calc_gradv(double ** p, int i, int j, int k, double gradv[SPACEDIM][SPACEDIM]);
void apply_implicit_sources(double ** p, double dtpush);
#else
void get_geom_src_terms(double NDP_PTR p);
void p_to_phys_src_start(double NDP_PTR p);
void p_to_phys_src_finish(double NDP_PTR p, double dtpush);
void calc_gradv(double NDP_PTR p, int i, int j, int k, double gradv[SPACEDIM][SPACEDIM]);
void apply_implicit_sources(double NDP_PTR p, double dtpush);
#endif
void calc_dvdx(int i, int j, int k, double gradv[SPACEDIM][SPACEDIM], double vavg[SPACEDIM], double xavg[NDIM]);
GPU_PRAGMA(omp declare target)
void p_to_u(double* p, double* u, const zone_geom* restrict g,int nhydro);
void u_to_p(double* u, double* p, const zone_geom* restrict g,int nhydro,double rho_floor,double e_floor);
int implicit_source_update(double *p, double *eos, zone_geom *gm, double dtpush);
int implicit_source_update_tianshu(double *u, double *eos, zone_geom *gm, double dtpush);

/* fluid.c */
void volume_average_pressure(const int i, const int j, const int k, double *pvec);
void hydro_p_to_u(const double* restrict p, double* restrict u, const double* restrict gcov,int nhydro);
void hydro_u_to_p(double* restrict u, double* restrict p, const double* restrict gcon,int nhydro,double rho_floor, double e_floor);
void hydro_stress_tensor(double *p, double press, double T[SPACEDIM][SPACEDIM], const zone_geom *g);
void u_to_p_edge(double *u, double *p, int i, int j, int k, int dir);
GPU_PRAGMA(omp end declare target)
void set_shock_flag(int cycle);
GPU_PRAGMA(omp declare target)
int transverse_shock(int ND_PTR ijk_to_I,int *sim_shock_flag, int *dj, int **dk,int i, int j, int k, int dir);
GPU_PRAGMA(omp end declare target)
void check_shock_flag(int id);
#if (USE_LINEAR_ALLOCATION==TRUE)
void fix_p(double ** p);
double is_shock(double ** p,int i, int j, int k, int dir);
double velocity_divergence(double ** p, int i, int j, int k, int dir);
void apply_lof(double ** p, double ** ph);
#else
//GPU_PRAGMA(omp declare target)
void fix_p(double NDP_PTR p);
double is_shock(double NDP_PTR p,int i, int j, int k, int dir);
double velocity_divergence(double NDP_PTR p, int i, int j, int k, int dir);
void apply_lof(double NDP_PTR p, double NDP_PTR ph);
//GPU_PRAGMA(omp end declare target)
#endif

/* gravity.c */
#if (USE_LINEAR_ALLOCATION==TRUE)
void gravity_start(double ** p);
void gravity_finish(double ** p);
void multipole_gravity(double ** p);
#else
void gravity_start(double NDP_PTR p);
void gravity_finish(double NDP_PTR p);
void multipole_gravity(double NDP_PTR p);
#endif

/* bc.c */
#if (USE_LINEAR_ALLOCATION==TRUE)
void reset_boundaries(double ** p,int on_gpu);
void reset_boundaries_gpu(double ** p,int on_gpu);
#else
void reset_boundaries(double NDP_PTR p,int on_gpu);
void reset_boundaries_gpu(double NDP_PTR p,int on_gpu);
#endif

/* mpi.c */
void init_mpi_setup();
void init_mpi_messaging();
void init_consistent_tracer_id();
int is_outside_rank(int index[], int rank);
void mpi_send_tracers(int stage);
void mpi_recv_tracers();
double mpi_min(double val);
double mpi_max(double val);
void mpi_global_reduce(double *arr, int len);
void mpi_global_reduce_uint(unsigned int *arr, int len);
void mpi_global_ireduce_start(double *arr, int len, int tag);
void mpi_global_ireduce_finish(double *arr, int len, int tag);
void mpi_bcast(double *arr, int len);
int mpi_io_proc();
void complete_mpi_communication(int half_step);
void complete_mpi_bc(int half_step); 
void ND_FUNC(build_types,NDIM)();
void test_boundaries();
void mpi_axis_average(double *arr, int len);
#if (USE_MPI==TRUE)
void mpi_bcast_gen(void *arr, int len, MPI_Datatype datatype);
proc_info *split_domain(int *ideal, int *actual, int *worst);
#endif
#if (USE_LINEAR_ALLOCATION==TRUE)
void sync_mpi_boundaries(double ** p);
#else
void sync_mpi_boundaries(double NDP_PTR p);
#endif

/* geometry.c */
void init_geometry();
GPU_PRAGMA(omp declare target)
double geom_dot(const double* restrict v1, const double* restrict v2, const double* restrict gcov);
void geom_lower(const double* restrict vcon, double* restrict vcov, const double* restrict gcov);
void geom_raise(const double* restrict vcov, double* restrict vcon, const double* restrict gcon);
GPU_PRAGMA(omp end declare target)
void vec_transform_to_xyz(double *vp, int i, int j, int k);
void ijk_to_x(int i, int j, int k, double x[]);
void x_to_rthphi(double x[SPACEDIM], double rvec[SPACEDIM]);

/* riemann solver */
GPU_PRAGMA(omp declare target)
inline void p_uflux(const int nhydro,const double* restrict p, const zone_geom* restrict g, double press, const int dir, double* restrict u, double* restrict f);
inline void pflux(const int nhydro,const double* restrict p, double press, const zone_geom* restrict g, const int dir, double* restrict f);
inline int riemann_flux_LLF(int nhydro, int nvars,double pl[], double pr[], const zone_geom const *g, const int dir, double flux[], double *sig_speed, double *vriemann);
int riemann_flux_HLLE(int nhydro,int nvars,double pl[], double pr[], const zone_geom *g, const int dir, double flux[], double *sig_speed, double *vriemann);
int riemann_flux_HLLC(const int nhydro, const int nvars,const double* restrict pl, const double* restrict pr, const zone_geom* restrict g, const int dir, double* restrict flux, double* restrict sig_speed, double* restrict vriemann);
GPU_PRAGMA(omp end declare target)

/* reconstruction */
GPU_PRAGMA(omp declare target)
void lin_recon(const double* restrict a, double* restrict al, double* restrict ar, 
               const double* restrict const beta, const double* restrict const Gamma,
               const double* restrict const x, const int istart, const int istop);
void para_recon(const double* restrict a, double* restrict al, double* restrict ar,
                       const double* restrict const alpha, const double* restrict const beta, const double* restrict const Gamma, 
                       const double* restrict const x, const int istart, const int istop);
void interp(int ninterp,int *interp_order,double p[NSIZE][2*NG], double pl[NSIZE][2*NG], double pr[NSIZE][2*NG], const double* restrict alpha, const double* restrict beta, const double* restrict Gamma, const double* restrict x, const int size);
GPU_PRAGMA(omp end declare target)
void init_interp_order();
//void interp(double p[NSIZE][NSIZE_GRID], double pl[NSIZE][NSIZE_GRID], double pr[NSIZE][NSIZE_GRID], const double restrict const *alpha, const double restrict const *beta, const double restrict const*Gamma, const double restrict const *x, int size);

/* utils.c */
double romb(double x0, double x1, double (*f)(double x));
void init_fornax_rand();
double fornax_rand();
double my_mod(double a, double b);
void allocate_memory();
int big_endian(void);
void bswap(void *vdat, int len, int cnt);
char *trim_whitespace(char *str);
void* malloc_rank1(uint64_t n1, uint64_t size);
void init_dendritic();
int ND_PTR dendritic_malloc_int();
double ND_PTR dendritic_malloc_double();
double NDP_PTR dendritic_malloc_vec(int vlen);
int get_prims_count();
double ran2(long int *idum);
int check_cons(double *u, char* msg, int ii, int jj, int kk);
double calculate_total_mass(int verbose);
double calculate_grav_energy();
double calculate_total_energy(int verbose);
double check_flux_differencing(int verbose, int check_shells);
int new_timer(char const * name);
Timer * get_timer(int timer_index);
void output_timers();
void free_timers();
double dynamic_root_find(double (*f)(double), double (*dfdx)(double),
                         double x0, double xlow, double xup, const double tol,
                         int max_iter);
double comp_x_from_r(double r, double xlow, double xhigh);
#if (USE_LINEAR_ALLOCATION==TRUE)
void fill_dead(double ** p);
void pack_prims(double ** p, double * out);
void unpack_prims(double const * in, double ** p);
void avg_poles_3d(double ** p);
int check_prim(double ** p, int ii, int jj, int kk, char* msg);
int find_prim_errors(double ** p, int check_ghost, char* msg);
void get_ijk_from_I(int II,int* ii,int* jj,int *kk);
#else
void fill_dead(double NDP_PTR p);
void pack_prims(double NDP_PTR p, double * out);
void unpack_prims(double const * in, double NDP_PTR p);
void avg_poles_3d(double NDP_PTR p);
int check_prim(double NDP_PTR p, int ii, int jj, int kk, char* msg);
int find_prim_errors(double NDP_PTR p, int check_ghost, char* msg);
#endif

/* problem dependent stuff (i.e. initialization, boundaries, sources, etc.) */
void init_grid();
void init_problem();
GPU_PRAGMA(omp declare target)
void prob_bounds(int i, int j, int k, double *p);
GPU_PRAGMA(omp end declare target)
void analysis_preloop();
void analysis_inloop();
void analysis_postloop();
void grav_accel(int ii, int jj, int kk, double *g);
#if (USE_LINEAR_ALLOCATION==TRUE)
void ext_src(double ** p, double dtpush);
void aux_src(double ** p, double dtpush);
#else
void ext_src(double NDP_PTR p, double dtpush);
void aux_src(double NDP_PTR p, double dtpush);
#endif

/* i/o */
void dump1();
void dump2();
void dump3();
void read3(char const * fname);
void dump_tracers();
void restart_dump();
void restart_dump_v1();
void restart_dump_v2();
void dump_ascii1(int is_restart);
void restart_read();
void restart_read_v1();
void restart_read_v2();
void run_log();
void parse_input(char *name);

/* eos */
void eos_init(double *table,char *name);
GPU_PRAGMA(omp declare target)
void eos_get_etaxnxp(const double *table,double rho, double temp, double ye, double *eta, double *xn, double *xp);
GPU_PRAGMA(omp end declare target)
#if (USE_LINEAR_ALLOCATION==TRUE)
void update_eos(const double *table,double ** p, double ** eos);
void update_eos_ghost(const double *table,double ** p, double ** eos);
int eos_fill_tianshu(const double *table,double **p, double **eos);
int eos_fill_ghost_tianshu(const double *table,double **p, double **eos);
#else
void update_eos(const double *table,double NDP_PTR p, double NDP_PTR eos);
void update_eos_ghost(const double *table,double NDP_PTR p, double NDP_PTR eos);
int eos_fill_tianshu(const double *table,double NDP_PTR p, double NDP_PTR eos);
int eos_fill_ghost_tianshu(const double *table,double NDP_PTR p, double NDP_PTR eos);
#endif

/* eos specifc routines found in eos/eos_xxxx.c */
GPU_PRAGMA(omp declare target)
int eos_fill(const double *table,double *p, double *eos);
GPU_PRAGMA(omp end declare target)
void eos_given_temp(const double *table,double *p, double *eos);
void eos_p_cs(const double *table,double *p, double *press, double *cs, double *Gamma);
void eos_p(const double *table,double *p, double *press);
double eos_cs(const double *table,double *p);
double u_given_rpx(const double *table,double rho, double press, double ye, double Tguess);
void get_lr_eos(const double *table,double rho, double temp, double ye, double *u, double *press, double *cs, double *Gamma);
double eos_get_Temp(const double *table,double rho, double e, double cmp, double Tguess);
double eos_rpx(const double *table,double rho, double press, double Ye, double Tguess, double *Gamma);
double eos_zero_point(const double *table,double rho, double ye);
double u_given_rtx(const double *table,double rho, double temp, double ye);

/* metric_x.c */
void init_interp_vol();
void init_volume();
void init_coords();
void init_area();
void init_conn();
void init_gcov();
void init_gcon();
void init_scale();
GPU_PRAGMA(omp declare target)
double r_of_x(r_of_x_struct_t rx_info,double x);
double th_of_x(th_of_x_struct_t thx_info,double x);
GPU_PRAGMA(omp end declare target)
double x_of_r(double r);
double dr_dx(double x);
double dth_dx(double x);
double z_of_x(double x);
double dz_dx(double x);
void x_to_ijk(double x[], int index[]);
void ijk_to_x(int i, int j, int k, double x[SPACEDIM]);
void ijk_to_Cart(int i, int j, int k, double *x);
double ijk_to_r(int i, int j, int k, double rhat[SPACEDIM]);
void calc_interp_coeffs1(int i, double *alpha, double *beta, double *Gamma, double *x);
void calc_interp_coeffs2(int i, int j, double *alpha, double *beta, double *Gamma, double *x);
double rtrans_solve(double dx, double rout);

/* radiation */
void init_rad_setup();
void rad_p_to_u(double *p, double *u, const zone_geom *gm);
void rad_u_to_p(double *u, double *p, const zone_geom *gm);
void enforce_flux_limit(double* restrict E, double* restrict Fcon, double* restrict Fcov, double* restrict Fred, const zone_geom* restrict gm, const int loc);
double planck(double freq, double T);
double compute_coeff(double kT, double rho_N, double Y_e);
int inelastic_nucleons(double* restrict u, const double dt, const zone_geom* restrict gm, double *initial_temp);
int rad_implicit_update(double* restrict p, const double dt, const zone_geom *gm, double *initial_temp, double intial_ye, int call);
int rad_inelastic_implicit_update(double* restrict u, const double dt, const zone_geom* restrict gm, double *initial_temp);
int rad_inelastic_explicit_update(double* restrict u, const double dt, const zone_geom* restrict gm, double *initial_temp);
int rad_inelastic_explicit_update_tianshu(double* restrict u, const double dt, const zone_geom* restrict gm, double *initial_temp,int scattype);
void rad_fluxes(double *pl, double *pr, double *flux, int dir, double *chi, double vadv, zone_geom *gmp, double *sig_speed, int dxfact);
void set_chat();
void freq_advection(const double dtpush, double *u, const zone_geom* restrict gm, double gradv[SPACEDIM][SPACEDIM], double* restrict fadv_flux
#if (GR_MONOPOLE==TRUE)
, double* restrict gr_fadv_flux, const double grg
#endif
);
void freq_advection_new(double *p, const zone_geom* restrict gm, double gradv[SPACEDIM][SPACEDIM], double* restrict fadv_flux);
void rad_src_update(double *u, double temperature, zone_geom *gm, double dtpush, int ii, int jj, int kk);
#if (USE_LINEAR_ALLOCATION==TRUE)
void rad_src(double ** p);
#else
void rad_src(double NDP_PTR p);
#endif

/* opacities */
void init_opac_emis();
void opac_emis(double rho, double T, double comp, double *kappa, double *jg, double *sc, double *delta, double *dk, double *dj, double *ds);
void inelastic_opac_emis(double *E, double *F, double rho, double T, double Ye, double *kappa, double *jg, double *sc);
void inelastic_opac_emis_tianshu(double *E, double *F, double rho, double T, double Ye, double *kappa, double *jg, double *sc,int scattype);
void opac_emis_freq_deriv(double rho, double T, double Ye, double *dkdlnu, double *djdlnu, double *dscdlnu);

/* quadrupole.c */
void quadrupole_init();
void quadrupole_free();
void quadrupole_start();
void quadrupole_end();
void quadrupole_dump();

// tracer.c
void init_tracers_by_mass(long int *rseed);
