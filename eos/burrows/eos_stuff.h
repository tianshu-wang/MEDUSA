
void collapse_init(double *table,const char *eos_name);

void get_numel(int *num_elements);

void finde(const double *table,double *f,int nelem,double *temp, double rho,double ye,int *jy_out,int *jq_out,int *jr_out);

void findthis(const double *table,double *f,int nelem,double temp,double rho,double ye,int *jy_out, int *jq_out,int *jr_out);

void eos_init(double *table,char *eos_name);

void min_e_given_rx(const double *table,double R, double Y,double *e);

#if(USE_LINEAR_ALLOCATION==TRUE)
void eos_given_rex_tianshu(const double *table,double **p, double **eos);
void eos_given_rex_ghost_tianshu(const double *table,double **p, double **eos);
#else
void eos_given_rex_tianshu(const double *table,double NDP_PTR p, double NDP_PTR eos);
void eos_given_rex_ghost_tianshu(const double *table,double NDP_PTR p, double NDP_PTR eos);
#endif

GPU_PRAGMA(omp declare target)
void eos_given_rex(const double *table,double *G, double *P, double *C, double *T, double *S, double *dpdr, double *dpde, double *R, double *e, double *Y);
void eos_given_rex_per_zone(const double *table, double *p, double *eos);
GPU_PRAGMA(omp end declare target)

void eos_given_rpx(const double *table,double *R, double *P, double *T, double *e, double *Y, double *G);

void eos_given_rsx(const double *table,double *G, double *P, double *C, double *T, double *e, double *R, double *S, double *Y);

void eos_s_given_rex(const double *table,double *S, double *R, double *e, double *T, double *Y);

void eos_given_rtx(const double *table,double *e, double *P, double *G, double *C, double *R, double *T, double *Y);

void eos_get_cv(const double *table,double *cv, double *R, double *T, double *Y);

GPU_PRAGMA(omp declare target)
void eos_get_etaxnxp_(const double *table,double *eta, double *xn, double *xp, double *R, double *T, double *Y);
GPU_PRAGMA(omp end declare target)

GPU_PRAGMA(omp declare target)
void eos(const double *table,int input, double *dens, double *temp, double *ye_in,
                 int npoints, double *pres, double *eint,
                 double *gam1, double *cs, double *entropy, double *cv, double *xn, double *xp,
                 int do_eos_diag);

void bisection(const double *table,double target,int nn,double *temp,double rho, double ye);
GPU_PRAGMA(omp end declare target)
