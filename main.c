#include "decs.h"
#include "defs.h"

#if (DEBUG==TRUE)
#include <fenv.h>
int feenableexcept(int);
#endif

int main(int argc, char *argv[])
{

  short dump_flag=FALSE,restart_flag=FALSE,pdump_flag=FALSE;
  int dd;
  long long int zcount,z[3];
  double dtold;

  short terminate_flag = FALSE;
  #if (USE_MPI==TRUE)
  MPI_Request terminate_req;
  #endif

  #if (USE_MPI==TRUE)
  #if (USE_OMP==TRUE)
  int provided;
  int required = MPI_THREAD_FUNNELED;
  MPI_Init_thread(&argc, &argv, required, &provided);
  if (provided < required) {
    fprintf(stderr,"Failed to acquire sufficient thread safety from MPI implementation\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  #else
  MPI_Init(&argc, &argv);
  #endif
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  #endif
  get_wctime(&tstart);

  if (mpi_io_proc()) fprintf(stderr, BUILD_STRING);

  #if (DEBUG==TRUE)
  feenableexcept(FE_INVALID);
  #endif

  #if (USE_GPU==TRUE)
GPU_PRAGMA(omp target teams num_teams(1))
  {
    if (omp_is_initial_device ()) {
      printf("!!! NO OFFLOAD !!!\n");
      //abort ();
    }
    printf("!!!~~~:)HACKATHON =====> OMP_NUM_THREADS=%d\n",omp_get_num_threads());
  }	
  #endif

  // set up the simulation parameters from the input file
  parse_input(argv[1]);

  // initialize a random number generator...just in case we need it
  init_fornax_rand();

  // set problem specific grid parameters, e.g. startx, dx
  init_grid();

  // set parameters of grid
  init_coords();

  // set up dendritic grid so we can figure out how to split the problem up with mpi
  if (mpi_io_proc()) fprintf(stderr,"initializing dendritic grid...");
  init_dendritic();
  if (mpi_io_proc()) fprintf(stderr,"done\n");

  // define the domain decomposition and establish neighbors
  if (mpi_io_proc()) fprintf(stderr,"initializing mpi setup...");
  init_mpi_setup();
  if (mpi_io_proc()) fprintf(stderr,"done\n");

  // allocate memory for this processors chunk of the problem
  allocate_memory(); // GPU finish
  quadrupole_init();

  // now that memory is allocated, we can compute memory displacements to build mpi derived data types
  #if (USE_MPI==TRUE)
  ND_FUNC(build_types,NDIM)();
  #endif

  // set all the geometry related variables, e.g. volume, area, gcov, etc.
  init_geometry();



  // set reasonable defaults
  t = 0.0;
  istep = 0;
  t_next_dump = dt_dump;
  t_next_pdump = dt_pdump;
  t_next_restart = dt_restart;
  i_next_restart = di_restart;
  if (detect_tbounce) {
    tbounce = -1.0;
    rho_cent_max = 0.0;
  }

  // set the problem specific initial data
  if (mpi_io_proc()) fprintf(stderr,"initializing problem...\n");
  #if DO_TIMING==TRUE
  Timer * t_init_pb = get_timer(new_timer("init_problem"));
  timer_start(t_init_pb);
  #endif
  init_problem();
  #if DO_TIMING==TRUE
  timer_stop(t_init_pb);
  #endif
  if (mpi_io_proc()) fprintf(stderr,"... problem done\n");

  init_consistent_tracer_id();

  // set the order of reconstruction for all variables
  init_interp_order();

  if (detect_tbounce && tbounce<0.0) set_tbounce();


  if(NSIZE<=ninterp){printf("NSIZE has to be greater than ninterp!\nNSIZE=%d, ninterp=%d.\nReset NSIZE in decs.h.\n",NSIZE,ninterp);exit(0);}
  if(NSIZE_GRID<=max_grid_dim){printf("NSIZE_GRID has to be greater than max_grid_dim!\nNSIZE=%d, max_grid_dim=%d.\nReset NSIZE_GRID in decs.h.\n",NSIZE_GRID,max_grid_dim);exit(0);}
  if(MULTIPOLE_LMAX<=multipole_lmax+1){printf("MULTIPOLE_LMAX has to be greater than multipole_lmax+1!\nMULTIPOLE_LMAX=%d, multipole_lmax+1=%d.\nReset MULTIPOLE_LMAX in decs.h.\n",MULTIPOLE_LMAX,multipole_lmax+1);exit(0);}



  // Port to GPU 
  GPU_PRAGMA(omp target enter data map(to:dj[-NG:n1+2*NG],nj[-NG:n1+2*NG]))
  GPU_PRAGMA(omp target enter data map(to:dk[-NG:n1+2*NG][-NG:n2+2*NG],nk[-NG:n1+2*NG][-NG:n2+2*NG]))
  GPU_PRAGMA(omp target enter data map(to:rcenter[-NG:n1+2*NG],redge[-NG:n1+2*NG+1],dr[-NG:n1+2*NG]))
  GPU_PRAGMA(omp target enter data map(to:bc[:nvars]))
  //GPU_PRAGMA(omp target enter data map(to:global_grid_dims[:SPACEDIM],istart[:SPACEDIM],dx[:SPACEDIM],istop[:SPACEDIM],startx[:SPACEDIM],my_grid_dims[:SPACEDIM]))

  GPU_PRAGMA(omp target enter data map(to:\
 	    sim_p[:cell_count_all][:nvars],sim_ph[:cell_count_all][:nvars],\
             sim_eos[:cell_count_all][:NEOS],sim_src[:cell_count_all][:nvars],\
 	    sim_shock_flag[:cell_count_all],sim_Phi[:cell_count_all],geom[:cell_count_all],\
 	    sim_fdir0[:nvars][:face_count_0],sim_vedgedir0[:SPACEDIM][:face_count_0]))
 	    //sim_fdir0[:face_count_0][:nvars],sim_vedgedir0[:face_count_0][:SPACEDIM]))
  #if(NDIM>1)
  GPU_PRAGMA(omp target enter data map(to:sim_fdir1[:nvars][:face_count_1],\
 			 sim_vedgedir1[:SPACEDIM][:face_count_1]))
  #endif
  #if(NDIM>2)
  GPU_PRAGMA(omp target enter data map(to:sim_fdir2[:nvars][:face_count_2],\
 			 sim_vedgedir2[:SPACEDIM][:face_count_2]))
  #endif
 
  GPU_PRAGMA(omp target enter data map(to:gr_lapse_edge[:n1+1],gr_lapse[:n1],\
			 alpha0[-NG:n1+2*NG],beta0[-NG:n1+2*NG],Gamma0[-NG:n1+2*NG],\
			 interp_order[:ninterp],Trad[:ngroups*SQR(SPACEDIM)]))
  //GPU_PRAGMA(omp target enter data map(to:pencil[:NSIZE][:2*NG],pleft[:NSIZE][:2*NG],pright[:NSIZE][:2*NG]))

#if(NDIM>1)
  GPU_PRAGMA(omp target enter data map(to:Gamma1s[:max_jrefine_level+1][-NG:(n2+2*NG)],\
			                  beta1s[:max_jrefine_level+1][-NG:(n2+2*NG)],\
					  alpha1s[:max_jrefine_level+1][-NG:(n2+2*NG)]))
#endif
#if(NDIM>2)
  GPU_PRAGMA(omp target enter data map(to:Gamma2s[:max_krefine_level+1][-NG:(n3+2*NG)],\
			                  beta2s[:max_krefine_level+1][-NG:(n3+2*NG)],\
					  alpha2s[:max_krefine_level+1][-NG:(n3+2*NG)]))
#endif
  //end
  int Itest,itest,jtest,ktest;
  int error=0;
  GPU_PRAGMA(omp target teams distribute parallel for reduction(+:error) map(tofrom:error))
  for(int II=0;II<cell_count_all;II++){
    GET_IJK_FROM_I(II,itest,jtest,ktest);
    Itest = ND_ELEM(ijk_to_I,itest,jtest,ktest);
    if(Itest!=II){error+=1;}
  }
  if(error>0){printf("Collapsed Cell-centered Indices Wrong!\n");exit(1);}
  

  // fill the initial eos data
  if (mpi_io_proc()) fprintf(stderr,"initializing eos data...");
  GPU_PRAGMA(omp target update to(sim_p[:cell_count_all][:nvars]))
  update_eos(eostable,sim_p,sim_eos);
  GPU_PRAGMA(omp target update from(sim_eos[:cell_count_all][:NEOS]))
  if (mpi_io_proc()) fprintf(stderr,"done\n");

  // dump initial conditions
  if (!restart_from_hdf) {
    if (mpi_io_proc()) fprintf(stderr,"dumping initial conditions...");
    ND_FUNC(dump,NDIM)();
    dump_tracers();
    if (mpi_io_proc()) fprintf(stderr,"done\n");
  }
  analysis_preloop();

  dtold = 1.0e100;
  dump_flag = FALSE;
  pdump_flag = FALSE;
  #if DO_TIMING==TRUE
  Timer * t_main_loop = get_timer(new_timer("main_loop"));
  timer_start(t_main_loop);
  Timer * t_analysis_inloop = get_timer(new_timer("analysis_inloop"));
  #endif



  Ncurrent = 0;
  while (t < tmax && istep < max_steps && terminate_flag == FALSE) {
    if (mpi_io_proc()) {
      get_wctime(&tcurr);
      telapsed = tcurr - tstart;
      terminate_flag = (telapsed > max_wtime) ? (TRUE) : (FALSE);
      if (terminate_flag) {
        fprintf(stderr,"Elapsed time: %g > %g, triggering termination\n",
            telapsed, max_wtime);
      }
    }
    #if (USE_MPI==TRUE)
    MPI_Ibcast(&terminate_flag, 1, MPI_SHORT, 0, MPI_COMM_WORLD, &terminate_req);
    #endif

    dump_flag = restart_flag = pdump_flag = FALSE;

#if (USE_LARGER_STEP==TRUE)
      Ncurrent = 0;
      Nskip = 1;
#endif

    if(Ncurrent==0){
      estimate_dt();
      if (dt < dtmin) {
        if (mpi_io_proc()) fprintf(stderr,"dt dropped below minimum, aborting on step %d:  %g %g\n", istep, dt, dtmin);
        ND_FUNC(dump,NDIM)();
        dump_tracers();
        quadrupole_dump();
        #if (USE_MPI==TRUE)
        MPI_Finalize();
        #endif
        exit(1);
      }
      if (istep==0) dt *= dt_init_safety;
      else dt = MIN(dtold*dt_grow_safety, dt);
      dtold = dt;
      if (t+Nskip*dt > t_next_dump) {
        dt = 1.000000001*(t_next_dump - t)/Nskip;
      }
      if (t+Nskip*dt > t_next_pdump) {
        dt = 1.000000001*(t_next_pdump - t)/Nskip;
      }
      if (t+Nskip*dt > tmax) {
        dt = 1.000000001*(tmax - t)/Nskip;
      }
      if (dt < dtmin) dt = dtmin;
    }

    if (mpi_io_proc() && istep%nstep_log == 0) {
      run_log();
    }
    if (istep > 0 && istep%nstep_timing == 0) {
      output_timers();
    }

    // evolve the whole system by dt
    step();
    t += dt;
    istep++;
    Ncurrent = (Ncurrent+1)%Nskip;


    #if (NDIM==1)
    if (restart_create_1d && tbounce>0.0 && t>=tbounce+restart_delay) {
      if (mpi_io_proc()) fprintf(stderr,"Writing 1-d ASCII restart file...\n");
      dump_ascii1(1);
      exit(1);
    }
    #endif

    if (istep%nstep_analysis==0) {
      update_cpu();
      #if (DO_TIMING==TRUE)
      timer_start(t_analysis_inloop);
      #endif
      analysis_inloop();
      #if (DO_TIMING==TRUE)
      timer_stop(t_analysis_inloop);
      #endif
      quadrupole_start();
      quadrupole_end();
    }

    if (fabs(t-t_next_dump)/dt < 1.0e-6) {
      update_cpu();
      ND_FUNC(dump,NDIM)();
      quadrupole_dump();
      t_next_dump += dt_dump;
      dump_flag = TRUE;
    }
    if (fabs(t-t_next_pdump)/dt < 1.0e-6) {
      update_cpu();
      dump_tracers();
      t_next_pdump += dt_pdump;
      pdump_flag = TRUE;
    }
    if (t >= t_next_restart || istep >= i_next_restart) {
      update_cpu();
      t_next_restart = t + dt_restart;
      i_next_restart = istep + di_restart;
      restart_dump();
      restart_flag = TRUE;
    }

    #if (USE_MPI==TRUE)
    MPI_Wait(&terminate_req, MPI_STATUS_IGNORE);
    #endif
  }
  #if (USE_MPI==TRUE)
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  #if DO_TIMING==TRUE
  timer_stop(t_main_loop);
  #endif

  if (mpi_io_proc()) {
    run_log();
  }

  if (!dump_flag) {
    ND_FUNC(dump,NDIM)();
    quadrupole_dump();
  }
  if (!pdump_flag) {
    dump_tracers();
  }
  if (!restart_flag) {
    restart_dump();
  }
  analysis_postloop();

  output_timers();

  #if (USE_MPI==TRUE)
  MPI_Finalize();
  #endif

  return 0;
}
