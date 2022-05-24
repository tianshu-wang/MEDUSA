#include "../../decs.h"
#include "./eos_stuff.h"
// Must define these since decs.h is not included...
//#define TRUE 1
//#define FALSE 0

#define TAB_ELEM(arr,irec,jy,jr,jt) arr[(irec-1)*ny*nr*nt+(jy-1)*nr*nt+(jr-1)*nt+(jt-1)]
#define DEFINE_EOSCONST_1 \
int numel = 16;\
int nt = 300;\
int nr = 300;\
int ny = 50;\
int tablesize = 16*300*300*50;\
double r1  =  1.e0;\
double r2  =  15.5e0;\
double t1  = -2.0e0;\
double t2  =  2.0e0;\
double t12 = -2.0e0;\
double t22 =  2.0e0;\
double y1c  =  0.035e0;\
double y2c  =  0.56e0;

//// number of elements in AB table
//static int numel = 16;
//// table dimensions
//static int nt = 300;
//static int nr = 300;
//// LS stellar collapse requires nt=300, nr=300, ny=50 (100MEV)
//static int ny = 50;
//static double r1  =  1.e0;
//static double r2  =  15.5e0;
//static double t1  = -2.0e0;
//static double t2  =  2.0e0;
//static double t12 = -2.0e0;
//static double t22 =  2.0e0;
//// To avoid name conflicts, y1c and y2c are renamed as y1cc and y2cc
//static double y1c  =  0.035e0;
//static double y2c  =  0.56e0;


// LS stellar collapse requires nt=300, nr=300, ny=50 (HI TEMP)
//  integer, parameter :: ny = 50
//  double precision :: r1  =  1.d0  //  4.0d0
//  double precision :: r2  =  15.0d0
//  double precision :: t1  = -2.0d0
//  double precision :: t2  =  1.82d0 //0.8d0
//  double precision :: t12 = -2.0d0
//  double precision :: t22 =  1.82d0
//  double precision :: y1c  =  0.035d0
//  double precision :: y2c  =  0.56d0

  // LS stellar collapse requires nt=300, nr=300, ny=50
//  integer, parameter :: ny = 50
//  double precision :: r1  =  1.d0  //  4.0d0
//  double precision :: r2  =  15.0d0
//  double precision :: t1  = -2.0d0
//  double precision :: t2  =  1.7d0 //0.8d0
//  double precision :: t12 = -2.0d0
//  double precision :: t22 =  1.7d0
//  double precision :: y1c  =  0.035d0
//  double precision :: y2c  =  0.56d0

  // Shen 300 300 54, t lower 0.01 MeV

//  integer, parameter :: ny = 54
//  double precision :: r1 = 1.0d0
//  double precision :: r2  =  15.0d0
//  double precision :: t1  = -2.0d0  // THIS IS FOR SHEN OCT 16
//  double precision :: t2 = 1.7d0
//  double precision :: t12 = -2.0d0  // THIS IS FOR SHEN OCT 16
//  double precision :: t22 = 1.7d0
//  double precision :: y1c = 0.03d0
//  double precision :: y2c = 0.56d0


void collapse_init(double *table,const char *eos_name){
  DEFINE_EOSCONST_1;
  int i, jy, jr, jt, irec, irecl, tndx, tndx0, ipass;
  int ncpu,mype,info,ierr;
//  double precision, allocatable :: table1d(:)
#if(USE_MPI==TRUE)
  // Init MPI
  info = MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  info = MPI_Comm_size(MPI_COMM_WORLD,&ncpu);
#else
  ncpu = 1;
  mype = 0;
#endif

//--------------------------------------------
//  TAB_ELEM(table,1,y,r,t)  :: energy per baryon      |
//  TAB_ELEM(table,2,y,r,t)  :: pressure               |
//  TAB_ELEM(table,3,y,r,t)  :: entropy per baryon     |
//  TAB_ELEM(table,4,y,r,t)  :: cv                     |
//  TAB_ELEM(table,5,y,r,t)  :: xn                     |
//  TAB_ELEM(table,6,y,r,t)  :: xp                     |
//  TAB_ELEM(table,7,y,r,t)  :: xa                     |
//  TAB_ELEM(table,8,y,r,t)  :: xh                     |
//  TAB_ELEM(table,9,y,r,t)  :: za                     |
//  TAB_ELEM(table,10,y,r,t) :: aw                     |
//  TAB_ELEM(table,11,y,r,t) :: muhat                  |
//  TAB_ELEM(table,12,y,r,t) :: gamma//                 |
//  TAB_ELEM(table,13,y,r,t) :: dhy                    |
//  TAB_ELEM(table,14,y,r,t) :: zht                    |
//  TAB_ELEM(table,15,y,r,t) :: mue                    |
//  TAB_ELEM(table,16,y,r,t) :: dpde                   |
//--------------------------------------------

  irecl = 8 * ny * nr * nt;

  if (mype==0){
    FILE *f = fopen(eos_name,"rb");
    size_t result;
    for(irec=0;irec<numel;irec++){
      fseek(f,irec*irecl,SEEK_SET);
      result = fread(&table[irec*nt*nr*ny],sizeof(double),nt*nr*ny,f);
    }

    for(jt=1;jt<=nt;jt++){
       for(jr=1;jr<=nr;jr++){
          for(jy=1;jy<=ny;jy++){
             TAB_ELEM(table,2,jy,jr,jt) = log(TAB_ELEM(table,2,jy,jr,jt));
          }
       }
    }
    fclose(f);
  }

#if(USE_MPI==TRUE)
    // Broadcast EOS table
  info = MPI_Bcast(table,numel*ny*nr*nt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
}

void get_numel(int *num_elements){
  DEFINE_EOSCONST_1;
  *num_elements = numel;
}


//----------------------------------------------------------------------

void finde(const double *table,double *f,int nelem,double *temp, double rho,double ye,int *jy_out,int *jq_out,int *jr_out){
  DEFINE_EOSCONST_1;
  int jt,jr,jy,nn,jtp,jrp,jyp;
  double rl,tl,drho,dtemp,dye;
  double dr,dt,dy;
  double coeff[4];
  double elo,ehi,tlo,thi;

  rl = log10(rho);
  tl = log10(*temp);

  dr = (r2-r1)/(nr-1);
  dt = (t2-t1)/(nt-1);
  dy = (y2c-y1c)/(ny-1);

  jr = 1 + (int)((rl-r1)/dr);
  jt = 1 + (int)((tl-t1)/dt);
  jy = 1 + (int)((ye-y1c)/dy);

  jyp = jy+1;
  jrp = jr+1;

  drho = (rl - (r1+(jr-1)*dr))/dr;
  dye = (ye - (y1c+(jy-1)*dy))/dy;

  if (jy==ny) {
     jyp = jy;
     dye = 0.0;
  }
  if (jr==nr) {
     jrp = jr;
     drho = 0.0;
  }

  // project into 1-D table
  coeff[0] = (1.0 - drho)*(1.0 - dye);
  coeff[1] = dye*(1.0 - drho);
  coeff[2] = (1.0 - dye)*drho;
  coeff[3] = dye*drho;

  // now find table entries that bound f(1), the energy we're after
  ehi = -1e99;
  while(1){
    elo = coeff[0]*TAB_ELEM(table,1,jy,jr,jt);
    elo = elo + coeff[1]*TAB_ELEM(table,1,jyp,jr,jt);
    elo = elo + coeff[2]*TAB_ELEM(table,1,jy,jrp,jt);
    elo = elo + coeff[3]*TAB_ELEM(table,1,jyp,jrp,jt);
    if (elo<f[0]||jt==1){
      break;
    }
    ehi = elo;
    jt = jt - 1;
  }
  if (elo>f[0]){
    printf("Off the bottom of the table in finde %e %e %e %e %e\n", rho, *temp, ye, f[0], elo);
  }

  if (ehi<-1e98){
    while(1){
      ehi = coeff[0]*TAB_ELEM(table,1,jy,jr,jt+1);
      ehi = ehi + coeff[1]*TAB_ELEM(table,1,jyp,jr,jt+1);
      ehi = ehi + coeff[2]*TAB_ELEM(table,1,jy,jrp,jt+1);
      ehi = ehi + coeff[3]*TAB_ELEM(table,1,jyp,jrp,jt+1);
      if (ehi>f[0]||jt==nt-1){
        break;
      }
      elo = ehi;
      jt = jt + 1;
    }
  }

  if (ehi<f[0]){
    printf("Off the top of the table in finde %e %e %e %e %e\n", rho, *temp, ye, f[0], ehi);
    //exit(1);
  }

  // now use a linear model.  must be linear so it's invertible
  tlo = t1 + (jt-1)*dt;
  thi = tlo + dt;
  // dtemp = de = (e-e0)/(e1-e0)
  dtemp = (f[0] - elo)/(ehi-elo);
  *temp = tlo + (thi-tlo)*dtemp;

  if (dtemp<0.0||dtemp>1.0){
    printf("dtemp out of bounds %e\n", dtemp);
    //exit(1);
  }

  jtp = jt+1;
  //f(:) = ttable(,:,jt)*(1-dtemp) + ttable(,:,jt+1)*dtemp
  for(int irec=1;irec<=nelem;irec++){
    f[irec-1] = coeff[0]*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
    f[irec-1] += coeff[1]*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
    f[irec-1] += coeff[2]*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
    f[irec-1] += coeff[3]*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
    f[irec-1] += coeff[0]*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
    f[irec-1] += coeff[1]*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
    f[irec-1] += coeff[2]*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
    f[irec-1] += coeff[3]*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);
  }

  f[1] = exp(f[1]);
  *temp = pow(10,*temp);

  *jy_out = jy;
  *jq_out = jt;
  *jr_out = jr;
}


void findthis(const double *table,double *f,int nelem,double temp,double rho,double ye,int *jy_out, int *jq_out,int *jr_out){
  DEFINE_EOSCONST_1;
  int jt,jr,jy,nn,jtp,jrp,jyp;
  double rl,tl,drho,dtemp,dye;
  double dr,dt,dy;
  double coeff[8];

  rl = log10(rho);
  tl = log10(temp);

  dr = (r2-r1)/(nr-1);
  dt = (t2-t1)/(nt-1);
  dy = (y2c-y1c)/(ny-1);

  jr = 1 + (int)((rl-r1)/dr);
  jt = 1+ (int)((tl-t1)/dt);
  jy = 1 + (int)((ye-y1c)/dy);

  jyp = jy+1;
  jrp = jr+1;
  jtp = jt+1;

  drho = (rl - (r1+(jr-1)*dr))/dr;
  dtemp = (tl - (t1+(jt-1)*dt))/dt;
  dye = (ye - (y1c+(jy-1)*dy))/dy;

  if (jy==ny){
     jyp = jy;
     dye = 0.0;
  }
  if (jr==nr){
     jrp = jr;
     drho = 0.0;
  }
  if (jt==nt){
     jtp = jt;
     dtemp = 0.0;
  }

  coeff[0] = (1.0 - drho)*(1.0 - dtemp)*(1.0 - dye);
  coeff[1] = (1.0 - drho)*(1.0 - dtemp)*dye;
  coeff[2] = drho*(1.0 - dtemp)*(1.0 - dye);
  coeff[3] = drho*(1.0 - dtemp)*dye;
  coeff[4] = (1.0 - drho)*dtemp*(1.0 - dye);
  coeff[5] = (1.0 - drho)*dtemp*dye;
  coeff[6] = drho*dtemp*(1.0 - dye);
  coeff[7] = drho*dtemp*dye;

  for(int irec=1;irec<=nelem;irec++){
    f[irec-1]  =coeff[0]*TAB_ELEM(table,irec,jy,jr,jt);
    f[irec-1] +=coeff[1]*TAB_ELEM(table,irec,jyp,jr,jt);
    f[irec-1] +=coeff[2]*TAB_ELEM(table,irec,jy,jrp,jt);
    f[irec-1] +=coeff[3]*TAB_ELEM(table,irec,jyp,jrp,jt);
    f[irec-1] +=coeff[4]*TAB_ELEM(table,irec,jy,jr,jtp);
    f[irec-1] +=coeff[5]*TAB_ELEM(table,irec,jyp,jr,jtp);
    f[irec-1] +=coeff[6]*TAB_ELEM(table,irec,jy,jrp,jtp);
    f[irec-1] +=coeff[7]*TAB_ELEM(table,irec,jyp,jrp,jtp);
  }
//  f(1) = exp(f(1))
  f[1] = exp(f[1]);

//  f(12) = (TAB_ELEM(table,12,jy,jr,jt)    + TAB_ELEM(table,12,jyp,jr,jt) &
//          + TAB_ELEM(table,12,jy,jrp,jt)  + TAB_ELEM(table,12,jrp,jrp,jt) &
//          + TAB_ELEM(table,12,jy,jr,jtp)  + TAB_ELEM(table,12,jyp,jr,jtp) &
//          + TAB_ELEM(table,12,jy,jrp,jtp) + TAB_ELEM(table,12,jyp,jrp,jtp))/8.d0

//  do nn=1,nelem
//
//     f(nn) =  coeff(1)*TAB_ELEM(table,nn,jy,jr,jt) &
//            + coeff(2)*TAB_ELEM(table,nn,jyp,jr,jt) &
//            + coeff(3)*TAB_ELEM(table,nn,jy,jrp,jt) &
//            + coeff(4)*TAB_ELEM(table,nn,jyp,jrp,jt) &
//            + coeff(5)*TAB_ELEM(table,nn,jy,jr,jtp) &
//            + coeff(6)*TAB_ELEM(table,nn,jyp,jr,jtp) &
//            + coeff(7)*TAB_ELEM(table,nn,jy,jrp,jtp) &
//            + coeff(8)*TAB_ELEM(table,nn,jyp,jrp,jtp)
//
//     if (nn==5||nn==6||nn==7||nn==8) then
//          f(nn)  = MIN(MAX(f(nn),0.d0),1.d0)
//     end if
//  enddo

  *jy_out = jy;
  *jq_out = jt;
  *jr_out = jr;
}

#define DEFINE_EOSCONST_2 \
int NP = 1;\
int npts = 1;\
int nspec = 0;\
double newton_tol = 1.0e-15;\
double ye_eos;\
double temp_eos;\
double den_eos;\
double e_eos;\
double p_eos;\
double h_eos;\
double cv_eos = 0;\
double cp_eos;\
double xn_eos;\
double xp_eos;\
double xne_eos;\
double eta_eos;\
double pele_eos;\
double dpdt_eos;\
double dpdr_eos;\
double dedr_eos;\
double dedt_eos;\
double gam1_eos;\
double   cs_eos;\
double    s_eos;\
double dsdt_eos;\
double dsdr_eos;\
double table_Tmin, table_Yemin;\
int eos_input_rt = 1;\
int eos_input_rs = 2;\
int eos_input_re = 5;\
int eos_input_cv = 6;\
int eos_input_etaxnxp = 7;\
int eos_input_rp = 8;\
int eos_input_emin = 9;\
double const_hack = 0.0;

/*
static int NP = 1;
static int npts = 1;
static int nspec = 0;
static double newton_tol = 1.0e-15;
static double ye_eos[1];
static double temp_eos[1];
static double den_eos[1];
static double e_eos[1];
static double p_eos[1];
static double h_eos[1];
static double cv_eos[1] = {0};
static double cp_eos[1];
static double xn_eos[1];
static double xp_eos[1];
static double xne_eos[1];
static double eta_eos[1];
static double pele_eos[1];
static double dpdt_eos[1];
static double dpdr_eos[1];
static double dedr_eos[1];
static double dedt_eos[1];
static double gam1_eos[1];
static double   cs_eos[1];
static double    s_eos[1];
static double dsdt_eos[1];
static double dsdr_eos[1];
static double table_Tmin, table_Yemin;
static int eos_input_rt = 1;   // density, temperature are inputs
static int eos_input_rs = 2;   // density, entropy are inputs
static int eos_input_re = 5;   // density, internal energy are inputs
static int eos_input_cv = 6;
static int eos_input_etaxnxp = 7;
static int eos_input_rp = 8;
static int eos_input_emin = 9;
static double const_hack = 0.0;

//never used
static double smallt;
static double smalld;
*/
static int initialized = 0;

// jbb added a constant to eint to have posistive internal energy
//  double,parameter :: const_hack = 9.3d0


void eos_init(double *table,char *eos_name){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
  /*
  if (present(small_temp)) then
     if (small_temp > 0.d0) then
        smallt = small_temp
     else
        // smallt = 5.d6
        // Since we're now in MeV, divide by 10^10
        smallt = 5.d-4
     end if
  else
     // smallt = 5.d6
     // Since we're now in MeV, divide by 10^10
     smallt = 5.d-4
  endif

  if (present(small_dens)) then
     if (small_dens > 0.d0) then
        smalld = small_dens
     else
        smalld = 1.d-5
     end if
  else
     smalld = 1.d-5
  end if
  */

  // initialize collapse
  collapse_init(table,eos_name);

  initialized = 1;

  table_Tmin = pow(10,t1);
  table_Yemin = y1c;

  GPU_PRAGMA(omp target enter data map(to:table[:tablesize]))
}

// never used
/*

  subroutine eos_get_small_temp(small_temp_out)

    double, intent(out) :: small_temp_out

    small_temp_out = smallt

  end subroutine eos_get_small_temp


  subroutine eos_get_small_dens(small_dens_out)

    double, intent(out) :: small_dens_out

    small_dens_out = smalld

  end subroutine eos_get_small_dens

*/

void min_e_given_rx(const double *table,double R, double Y,double *e){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    // local variables
    int do_eos_diag;

    do_eos_diag = 0;

    den_eos  = R;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos   = Y;

    eos(table,eos_input_emin, &den_eos, &temp_eos, &ye_eos,
             npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
             &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *e = e_eos;
}

#if(USE_LINEAR_ALLOCATION==TRUE)
void eos_given_rex_tianshu(const double *table,double **p,double **eos)
#else
void eos_given_rex_tianshu(const double *table,double NDP_PTR p,double NDP_PTR eos)
#endif
{
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
//eos_input_re = 5;
  int jt,jr,jy,nn,jtp,jrp,jyp;
  double rl,tl,drho,dtemp,dye;
  double dr,dt,dy;
  double coeff0,coeff1,coeff2,coeff3;
  double elo,ehi,tlo,thi;
  int irec;

  GPU_PRAGMA(omp target teams distribute parallel for)
  for(int II=0;II<cell_count;II++) {
    temp_eos = NDP_ELEM_LINEAR_REORDERED(eos,ii,jj,kk,TEMP);
    den_eos  = NDP_ELEM_LINEAR_REORDERED(p,ii,jj,kk,RHO);
    e_eos    = NDP_ELEM_LINEAR_REORDERED(p,ii,jj,kk,UU)/den_eos;
    ye_eos   = NDP_ELEM_LINEAR_REORDERED(p,ii,jj,kk,YE);
    if(temp_eos<1.0e-10){temp_eos=1.0;}

    e_eos = e_eos / 0.95655684e18 - const_hack;
    rl = log10(den_eos);
    tl = log10(temp_eos);
    dr = (r2-r1)/(nr-1);
    dt = (t2-t1)/(nt-1);
    dy = (y2c-y1c)/(ny-1);
    jr = 1 + (int)((rl-r1)/dr);
    jt = 1 + (int)((tl-t1)/dt);
    jy = 1 + (int)((ye_eos-y1c)/dy);
    jyp = jy+1;
    jrp = jr+1;
    drho = (rl - (r1+(jr-1)*dr))/dr;
    dye = (ye_eos - (y1c+(jy-1)*dy))/dy;
    if (jy==ny) {
       jyp = jy;
       dye = 0.0;
    }
    if (jr==nr) {
       jrp = jr;
       drho = 0.0;
    }
    
    // project into 1-D table
    coeff0 = (1.0 - drho)*(1.0 - dye);
    coeff1 = dye*(1.0 - drho);
    coeff2 = (1.0 - dye)*drho;
    coeff3 = dye*drho;
    
    // now find table entries that bound e_eos, the energy we're after
    ehi = -1e99;
    while(1){
      elo = coeff0*TAB_ELEM(table,1,jy,jr,jt);
      elo = elo + coeff1*TAB_ELEM(table,1,jyp,jr,jt);
      elo = elo + coeff2*TAB_ELEM(table,1,jy,jrp,jt);
      elo = elo + coeff3*TAB_ELEM(table,1,jyp,jrp,jt);
      if (elo<e_eos||jt==1){
        break;
      }
      ehi = elo;
      jt = jt - 1;
    }
    if (elo>e_eos){
      printf("Off the bottom of the table in finde %d %e %e %e %e %e\n",II,den_eos,temp_eos,ye_eos,e_eos,elo);
      //exit(1);
    }
    
    if (ehi<-1e98){
      while(1){
        ehi = coeff0*TAB_ELEM(table,1,jy,jr,jt+1);
        ehi = ehi + coeff1*TAB_ELEM(table,1,jyp,jr,jt+1);
        ehi = ehi + coeff2*TAB_ELEM(table,1,jy,jrp,jt+1);
        ehi = ehi + coeff3*TAB_ELEM(table,1,jyp,jrp,jt+1);
        if (ehi>e_eos||jt==nt-1){
          break;
        }
        elo = ehi;
        jt = jt + 1;
      }
    }
    
    if (ehi<e_eos){
      printf("Off the top of the table in finde %d %e %e %e %e %e\n",II,den_eos,temp_eos,ye_eos,e_eos,ehi);
      //exit(1);
    }
    
    // now use a linear model.  must be linear so it's invertible
    tlo = t1 + (jt-1)*dt;
    thi = tlo + dt;
    // dtemp = de = (e-e0)/(e1-e0)
    dtemp = (e_eos - elo)/(ehi-elo);
    temp_eos = tlo + (thi-tlo)*dtemp;
    
    if (dtemp<0.0||dtemp>1.0){
      printf("dtemp out of bounds %e\n", dtemp);
    }
    
    jtp = jt+1;

    irec=2;
    p_eos  = coeff0*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
    p_eos += coeff1*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
    p_eos += coeff2*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
    p_eos += coeff3*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
    p_eos += coeff0*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
    p_eos += coeff1*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
    p_eos += coeff2*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
    p_eos += coeff3*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);
    // pressure (MeV/fm^3) -> dyn/cm^2
    p_eos = exp(p_eos)*1.60217733e33;
    
    irec=3;
    s_eos  = coeff0*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
    s_eos += coeff1*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
    s_eos += coeff2*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
    s_eos += coeff3*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
    s_eos += coeff0*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
    s_eos += coeff1*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
    s_eos += coeff2*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
    s_eos += coeff3*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);

    irec=12;
    gam1_eos  = coeff0*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
    gam1_eos += coeff1*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
    gam1_eos += coeff2*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
    gam1_eos += coeff3*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
    gam1_eos += coeff0*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
    gam1_eos += coeff1*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
    gam1_eos += coeff2*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
    gam1_eos += coeff3*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);


    temp_eos = pow(10,temp_eos);
    cs_eos = sqrt((gam1_eos)*(p_eos)/(den_eos));

    NDP_ELEM_LINEAR_REORDERED(eos,ii,jj,kk,PRESS) = p_eos;
    NDP_ELEM_LINEAR_REORDERED(eos,ii,jj,kk,CS)    = cs_eos;
    NDP_ELEM_LINEAR_REORDERED(eos,ii,jj,kk,TEMP)  = temp_eos;
    NDP_ELEM_LINEAR_REORDERED(eos,ii,jj,kk,ENT)   = s_eos;
    NDP_ELEM_LINEAR_REORDERED(eos,ii,jj,kk,GAMMA) = gam1_eos;
  }
}

#if(USE_LINEAR_ALLOCATION==TRUE)
void eos_given_rex_ghost_tianshu(const double *table,double **p,double **eos)
#else
void eos_given_rex_ghost_tianshu(const double *table,double NDP_PTR p,double NDP_PTR eos)
#endif
{
  int ii=0,jj=0,kk=0,success;
  int jstart,jstop,kstart,kstop;

  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
  int jt,jr,jy,nn,jtp,jrp,jyp;
  double rl,tl,drho,dtemp,dye;
  double dr,dt,dy;
  double coeff0,coeff1,coeff2,coeff3;
  double elo,ehi,tlo,thi;
  int irec;

  //GPU_PRAGMA(omp target)
  {
    GPU_PRAGMA(omp target teams distribute parallel for)
    for(int II=cell_count;II<cell_count_all;II++){
      GET_IJK_FROM_I(II,ii,jj,kk);
      int ghost0=0,ghost1=0,ghost2=0;
      if(ii<istart[0]||ii>=istop[0]){ghost0=1;}
#if (NDIM>1)
      if(jj<JS(ii,istart[1])||jj>=JS(ii,istop[1])){ghost1=1;}
#endif
#if (NDIM==3)
      if(kk<KS(ii,jj,istart[2])||kk>=KS(ii,jj,istop[2])){ghost2=1;}
#endif
      if(ghost0+ghost1+ghost2>1){continue;}
      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    }
    //// 0-direction ghost zones
    //for (ii=istart[0]-NG; ii<istart[0]; ii++) {
    //  #if (NDIM>1)
    //  JSLOOP(ii,jj) {
    //  #endif
    //    #if (NDIM==3)
    //    KSLOOP(ii,jj,kk) {
    //    #endif
    //      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    //    #if (NDIM==3)
    //    }
    //    #endif
    //  #if (NDIM>1)
    //  }
    //  #endif
    //}
    //for (ii=istop[0]; ii<istop[0]+NG; ii++) {
    //  #if (NDIM>1)
    //  JSLOOP(ii,jj) {
    //  #endif
    //    #if (NDIM==3)
    //    KSLOOP(ii,jj,kk) {
    //    #endif
    //      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    //    #if (NDIM==3)
    //    }
    //    #endif
    //  #if (NDIM>1)
    //  }
    //  #endif
    //}

    //#if (NDIM>1)
    //// 1-direction ghost zones
    //ISLOOP(ii) {
    //  jstart = JS(ii,istart[1]);
    //  for (jj=jstart-NG; jj<jstart; jj++) {
    //    #if(NDIM==3)
    //    KSLOOP(ii,jj,kk) {
    //    #endif
    //      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    //    #if(NDIM==3)
    //    }
    //    #endif
    //  }
    //}
    //ISLOOP(ii) {
    //  jstop = JS(ii,istop[1]);
    //  for (jj=jstop; jj<jstop+NG; jj++) {
    //    #if(NDIM==3)
    //    KSLOOP(ii,jj,kk) {
    //    #endif
    //      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    //    #if(NDIM==3)
    //    }
    //    #endif
    //  }
    //}
    //#endif // NDIM>1 //

    //#if (NDIM==3)
    //// 2-direction ghost zones
    //ISLOOP(ii) {
    //  JSLOOP(ii,jj) {
    //    kstart = KS(ii,jj,istart[2]);
    //    for (kk=kstart-NG; kk<kstart; kk++) {
    //      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    //    }
    //  }
    //}
    //ISLOOP(ii) {
    //  JSLOOP(ii,jj) {
    //    kstop = KS(ii,jj,istop[2]);
    //    for (kk=kstop; kk<kstop+NG; kk++) {
    //      eos_given_rex_per_zone(table,ND_ELEM_LINEAR(p,ii,jj,kk), ND_ELEM_LINEAR(sim_eos,ii,jj,kk));
    //    }
    //  }
    //}
    //#endif
  } // omp parallel
}

void eos_given_rex_per_zone(const double *table, double *p, double *eos) {
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
  int jt,jr,jy,nn,jtp,jrp,jyp;
  double rl,tl,drho,dtemp,dye;
  double dr,dt,dy;
  double coeff0,coeff1,coeff2,coeff3;
  double elo,ehi,tlo,thi;
  int irec;

  temp_eos = eos[TEMP];
  den_eos  = p[RHO];
  e_eos    = p[UU]/den_eos;
  ye_eos   = p[YE];
  if(temp_eos<1.0e-10){temp_eos=1.0;}

  e_eos = e_eos / 0.95655684e18 - const_hack;
  rl = log10(den_eos);
  tl = log10(temp_eos);
  dr = (r2-r1)/(nr-1);
  dt = (t2-t1)/(nt-1);
  dy = (y2c-y1c)/(ny-1);
  jr = 1 + (int)((rl-r1)/dr);
  jt = 1 + (int)((tl-t1)/dt);
  jy = 1 + (int)((ye_eos-y1c)/dy);
  jyp = jy+1;
  jrp = jr+1;
  drho = (rl - (r1+(jr-1)*dr))/dr;
  dye = (ye_eos - (y1c+(jy-1)*dy))/dy;
  if (jy==ny) {
     jyp = jy;
     dye = 0.0;
  }
  if (jr==nr) {
     jrp = jr;
     drho = 0.0;
  }
  
  // project into 1-D table
  coeff0 = (1.0 - drho)*(1.0 - dye);
  coeff1 = dye*(1.0 - drho);
  coeff2 = (1.0 - dye)*drho;
  coeff3 = dye*drho;
  
  // now find table entries that bound e_eos, the energy we're after
  ehi = -1e99;
  while(1){
    elo = coeff0*TAB_ELEM(table,1,jy,jr,jt);
    elo = elo + coeff1*TAB_ELEM(table,1,jyp,jr,jt);
    elo = elo + coeff2*TAB_ELEM(table,1,jy,jrp,jt);
    elo = elo + coeff3*TAB_ELEM(table,1,jyp,jrp,jt);
    if (elo<e_eos||jt==1){
      break;
    }
    ehi = elo;
    jt = jt - 1;
  }
  if (elo>e_eos){
    printf("Off the bottom of the table in finde %e %e %e %e %e\n",den_eos,temp_eos,ye_eos,e_eos,elo);
  }
  
  if (ehi<-1e98){
    while(1){
      ehi = coeff0*TAB_ELEM(table,1,jy,jr,jt+1);
      ehi = ehi + coeff1*TAB_ELEM(table,1,jyp,jr,jt+1);
      ehi = ehi + coeff2*TAB_ELEM(table,1,jy,jrp,jt+1);
      ehi = ehi + coeff3*TAB_ELEM(table,1,jyp,jrp,jt+1);
      if (ehi>e_eos||jt==nt-1){
        break;
      }
      elo = ehi;
      jt = jt + 1;
    }
  }
  
  if (ehi<e_eos){
    printf("Off the top of the table in finde %e %e %e %e %e\n",den_eos,temp_eos,ye_eos,e_eos,ehi);
  }
  
  // now use a linear model.  must be linear so it's invertible
  tlo = t1 + (jt-1)*dt;
  thi = tlo + dt;
  // dtemp = de = (e-e0)/(e1-e0)
  dtemp = (e_eos - elo)/(ehi-elo);
  temp_eos = tlo + (thi-tlo)*dtemp;
  
  if (dtemp<0.0||dtemp>1.0){
    printf("dtemp out of bounds %e\n", dtemp);
  }
  
  jtp = jt+1;

  irec=2;
  p_eos  = coeff0*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
  p_eos += coeff1*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
  p_eos += coeff2*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
  p_eos += coeff3*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
  p_eos += coeff0*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
  p_eos += coeff1*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
  p_eos += coeff2*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
  p_eos += coeff3*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);
  // pressure (MeV/fm^3) -> dyn/cm^2
  p_eos = exp(p_eos)*1.60217733e33;
  
  irec=3;
  s_eos  = coeff0*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
  s_eos += coeff1*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
  s_eos += coeff2*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
  s_eos += coeff3*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
  s_eos += coeff0*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
  s_eos += coeff1*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
  s_eos += coeff2*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
  s_eos += coeff3*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);

  irec=12;
  gam1_eos  = coeff0*(1-dtemp)*TAB_ELEM(table,irec,jy,jr,jt);
  gam1_eos += coeff1*(1-dtemp)*TAB_ELEM(table,irec,jyp,jr,jt);
  gam1_eos += coeff2*(1-dtemp)*TAB_ELEM(table,irec,jy,jrp,jt);
  gam1_eos += coeff3*(1-dtemp)*TAB_ELEM(table,irec,jyp,jrp,jt);
  gam1_eos += coeff0*dtemp*TAB_ELEM(table,irec,jy,jr,jtp);
  gam1_eos += coeff1*dtemp*TAB_ELEM(table,irec,jyp,jr,jtp);
  gam1_eos += coeff2*dtemp*TAB_ELEM(table,irec,jy,jrp,jtp);
  gam1_eos += coeff3*dtemp*TAB_ELEM(table,irec,jyp,jrp,jtp);


  temp_eos = pow(10,temp_eos);
  cs_eos = sqrt((gam1_eos)*(p_eos)/(den_eos));

  eos[PRESS] = p_eos;
  eos[CS]    = cs_eos;
  eos[TEMP]  = temp_eos;
  eos[ENT]   = s_eos;
  eos[GAMMA] = gam1_eos;
}


void eos_given_rex(const double *table,double *G, double *P, double *C, double *T, double *S, double *dpdr, double *dpde, double *R, double *e, double *Y){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    double e0;
    int do_eos_diag;

    do_eos_diag = 0;

    temp_eos = *T;
    den_eos  = *R;
    e_eos    = *e;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos   = *Y;

    eos(table,eos_input_re, &den_eos, &temp_eos, &ye_eos,
        npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
        &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *G = gam1_eos;
    *P = p_eos;
    *C = cs_eos;
    *T = temp_eos;
    *S = s_eos;

// hack
    *dpdr = 0.0;
    *dpde = 0.0;
// end hack
}

//*************************************************************

void eos_given_rpx(const double *table,double *R, double *P, double *T, double *e, double *Y, double *G){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    // local variables
    int do_eos_diag;
    do_eos_diag = 0;

    temp_eos = *T;
    den_eos  = *R;
    p_eos    = *P;

    if(do_eos_diag) printf("eos_rpx rho, P, T: %e %e %e\n", *R, *P, *T);

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos   = *Y;

    eos(table,eos_input_rp, &den_eos, &temp_eos, &ye_eos,
             npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
             &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *e = e_eos;
    *G = gam1_eos;

    if(do_eos_diag) printf("Energy2c is %e\n", (*e)*(*R));
}

//*************************************************************

void eos_given_rsx(const double *table,double *G, double *P, double *C, double *T, double *e, double *R, double *S, double *Y){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    int do_eos_diag;
    do_eos_diag = 0;

    temp_eos = *T;
    den_eos  = *R;
    s_eos    = *S;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos   = *Y;

    eos(table,eos_input_rs, &den_eos, &temp_eos, &ye_eos,
             npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
             &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *G = gam1_eos;
    *P = p_eos;
    *C = cs_eos;
    *T = temp_eos;
    *e = e_eos;

// hack
    //dpdr = 0.d0
    //dpde = 0.d0
// end hack
}


//*************************************************************

void eos_s_given_rex(const double *table,double *S, double *R, double *e, double *T, double *Y){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    int do_eos_diag;
    do_eos_diag = 0;

    temp_eos = *T;
    den_eos  = *R;
    e_eos    = *e;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos   = *Y;

    eos(table,eos_input_re, &den_eos, &temp_eos, &ye_eos,
          npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
             &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *S = s_eos;
    *T = temp_eos;
}

//*************************************************************

void eos_given_rtx(const double *table,double *e, double *P, double *G, double *C, double *R, double *T, double *Y){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    double e0;
    int do_eos_diag;
    do_eos_diag = 0;

    den_eos  = *R;
    temp_eos = *T;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos   = *Y;

    eos(table,eos_input_rt, &den_eos, &temp_eos, &ye_eos,
             npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
             &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *P = p_eos;
    *e = e_eos;
    *G = gam1_eos;
    *C = cs_eos;

    //call min_e_given_RX(R, Y, e0)
    //e = e - e0
}


//*************************************************************

void eos_get_cv(const double *table,double *cv, double *R, double *T, double *Y){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    int do_eos_diag;
    do_eos_diag = 0;

    den_eos  = *R;
    temp_eos = *T;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos = Y[nspec];

    eos(table,eos_input_cv, &den_eos, &temp_eos, &ye_eos,
         npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
         &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    *cv = cv_eos;
}



//*************************************************************

void eos_get_etaxnxp_(const double *table,double *eta, double *xn, double *xp, double *R, double *T, double *Y){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
    int do_eos_diag;
    do_eos_diag = 0;

    den_eos = *R;
    temp_eos = *T;

    // NOTE: Ye comes in as the (nspec+1) component of Y
    ye_eos = *Y;

    eos(table,eos_input_etaxnxp, &den_eos, &temp_eos, &ye_eos,
         npts, &p_eos, &e_eos, &gam1_eos, &cs_eos, &s_eos,
         &cv_eos, &xn_eos, &xp_eos, do_eos_diag);

    // hack puts eta in cs
    *eta = cs_eos/(*T);
    *xn = xn_eos;
    *xp = xp_eos;

}

//*************************************************************
//*************************************************************

// a wrapper for interacting with Burrows eos
void eos(const double *table,int input, double *dens, double *temp, double *ye_in,
                 int npoints, double *pres, double *eint,
                 double *gam1, double *cs, double *entropy, double *cv, double *xn, double *xp,
                 int do_eos_diag){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
// dens     -- mass density (g/cc)
// temp     -- temperature (K)
// npoints  -- number of elements in i/o arrays
// ye       -- Ye
// pres     -- pressure (dyn/cm**2)
// enthalpy -- enthalpy (erg/g)
// eint     -- internal energy (erg/g)
//   Local variables
    double error;
    double energy_want;
    double entropy_want;
    double pressure_want;
    double temp_row;
    //double cv(npoints)
    double f[numel];
    double tempmev;
    //double xn(npoints)
    //double xp(npoints)
    double ye;

    double ttol = 1.e-8;
    double dtol = 1.e-8;
    double stol = 1.e-8;
    int max_newton = 40;

    int i, k, jy, jq, jr, icount;
    int nel;
    int iter, niter;

    double min_temp, max_temp, temp_hold;
    double temp0,depdt,dppdt,dt,dd;
    double rho,to,tc,epo,epc,ppo,ppc;


// table input:
//   - temp -- temperature (MeV)
//   - dens -- mass density (g/cc)
//   - ye   -- Ye
//
// table output:
//   - f(1)  -- energy (mev/baryon)
//   - f(2)  -- pressure (mev/fm^3)
//   - f(3)  -- entropy (per baryon per Boltzmann's constant)
//   - f(4)  -- specific heat at constant volume (ergs/MeV/gm)
//   - f(5)  -- mass fraction of neutrons
//   - f(6)  -- mass fraction of protons
//   - f(7)  -- mass fraction alpha particles
//   - f(8)  -- mass fraction heavy nuclei
//   - f(9)  -- not needed
//   - f(10) -- not needed
//   - f(11) -- not needed
//   - f(12) -- not needed
//   - f(13) -- not needed
//   - f(14) -- not needed
//   - f(15) -- not needed
//   - f(16) -- not needed [d pressure / d energy (??)]

    min_temp = pow(10.e0,t1);
    max_temp = pow(10.e0,t2);

    // Check if density within bounds
       if ( ((*dens)<pow(10.e0,r1))|| ((*dens)> pow(10.e0,r2)) ){
          printf("DENSITY OUT OF BOUNDS %e\n",(*dens));
          printf("LIMITS ARE            %e %e\n",pow(10.e0,r1),pow(10.e0,r2));
          printf("FROM CALL WITH INPUT  %d\n",input);
       }

    // Enforce ye within bounds
       (ye) = (*ye_in);
       if ( ((ye)<y1c)|| ((ye)>y2c) ){
// Ye out of range can happen sometimes during problems run with radiation.
// This may be an indicator of insufficient resolution or missing physics,
// but for now suppress the per-point printing that was choking the output
// files and set Ye silently to the edge of the valid range:
          (ye) = MAX(MIN((*ye_in),y2c),y1c);
//         stop
       }

    if (input == eos_input_emin){
          findthis(table,f,numel,min_temp+1.e-6*fabs(min_temp),(*dens),(ye),
                        &jy,&jq,&jr);
          (*eint) = (f[0]+const_hack) * 0.95655684e18;
    }
    else if (input == eos_input_rt){

//---------------------------------------------------------------------------
// input = 1: dens, temp, and ye are inputs
//---------------------------------------------------------------------------

       // NOTE: we assume the temp is already in MeV
       // tempmev = temp(k) / 1.160445d10


          if ((*temp) < min_temp || (*temp) > max_temp){
              printf("TEMP OUT OF BOUNDS WITH RTX CALL %e\n", (*temp));
              printf("LIMITS ARE: T1 T2 %e %e\n",min_temp, max_temp);
              //exit(1);
          }

          findthis(table,f,numel,(*temp),(*dens),(ye),&jy,&jq,&jr);

          // energy (MeV) per baryon -> ergs/g
//         eint(k) = f(1) * 0.95655684d18
// jbb hack
          (*eint) = (f[0]+const_hack) * 0.95655684e18;
          //eint(k) = (f(1)+const_hack) * .96485337d18

            (*cv) = f[3]; // * 1.160445d10

          // pressure (MeV/fm^3) -> dyn/cm^2
          (*pres) = f[1] * 1.60217733e33;

          (*gam1) = f[11];
            (*cs) = sqrt((*gam1)*(*pres)/(*dens));

          (*entropy) = f[2];
    }
    else if (input == eos_input_re){
            f[0] = (*eint) / 0.95655684e18 - const_hack;
            finde(table,f,numel,&(*temp),(*dens),(ye),&jy,&jq,&jr);

            // pressure (MeV/fm^3) -> dyn/cm^2
            (*pres) = f[1] * 1.60217733e33;

            (*gam1) = f[11];

            (*cs) = sqrt((*gam1)*(*pres)/(*dens));

            //gam1(k) = f(11)
            //cs(k) = f(15)

            (*entropy) = f[2];
    }
    // THIS IS BASED ON ADAM'S ITERATION SCHEME
    else if (input == 100*eos_input_re){

//---------------------------------------------------------------------------
// input = 5: dens, energy, and ye are inputs
//---------------------------------------------------------------------------

       // NOTE: we assume the temp is already in MeV
          // load initial guess
          temp_row = (*temp);

          if (do_eos_diag==1) printf("T/D INIT %e %e \n",(*temp),(*dens));

          // want to converge to the given energy (erg/g)
          // energy (MeV) per baryon -> ergs/gm
          //   MeV   1.602176462e-6 erg   6.02214199e23
          //   --- x ------------------ x -------------
          //   bar          1 MeV            1 mole

//         energy_want(k) = eint(k) / .95655684d18
// jbb
          energy_want = (*eint) / 0.95655684e18 - const_hack;
          //energy_want(k) = eint(k) / .96485337d18

          if (do_eos_diag) printf("WANT e (erg/g) %e\n",energy_want);

          icount = 0;

          // This is the initial guess
          temp0   = (*temp);

          findthis(table,f,numel,max_temp,(*dens),(ye),&jy,&jq,&jr);
          if (energy_want > f[0]){
            printf( "eos_given_REX: energy too high :        %e\n",energy_want);
            printf( "using rho temp ye %e %e %e\n",(*dens), max_temp, (ye));
            printf( "max possible energy given this density  %e\n",f[0]);
            printf( "Setting temperature to max_temp %e\n",max_temp);
            temp_row = max_temp;
            //exit(1);
            //go to 10
//           call bl_error('EOS: energy too high for table')
          }

          findthis(table,f,numel,min_temp+1.e-6*fabs(min_temp),(*dens),(ye),&jy,&jq,&jr);
          if (energy_want < f[0]){
            printf("eos_given_REX: energy too low :        %e\n",energy_want);
            printf("using rho temp ye %e %e %e \n",(*dens), min_temp+1.e-6*fabs(min_temp), (ye));
            printf("min possible energy given this density  %e\n",f[0]);
            printf("Setting temperature to min_temp %e\n",min_temp*1.000001e0);
            temp_row = min_temp + 1.e-6*fabs(min_temp);
            //exit(1);
            //go to 10
//           call bl_error('EOS: energy too low for table')
          }

          rho = (*dens);

          // Make "old" value = 1.01 * incoming guess
          to = MIN(temp_row * 1.01e0, max_temp);
          findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
          epo = f[0];

          // Make "current" value = incoming guess
          tc  = temp_row;
          depdt  = epo/to;

          for(i = 1;i<=max_newton;i++){

             findthis(table,f,numel,tc,rho,(ye),&jy,&jq,&jr);
             epc = f[0];

             if (do_eos_diag==1){
                printf("ITER %d\n",i);
                printf("TEMP EP OLD %e %e\n",to,epo);
                printf("TEMP EP NEW %e %e\n",tc,epc);
             }

             // Calculate the slope of energy vs temp
             // Be careful//  This is numerically unstable//
             // Just freeze the derivative once the domain is too small
             if(fabs((tc-to)/tc)>1e-8) depdt  = (epc-epo)/(tc-to);


             dd = (energy_want-epc)/depdt;

             // Add that much to the current guess.
             if (do_eos_diag) printf("EPC - EPO %e %e %e\n",epc, epo, epc-epo);
             if (do_eos_diag) printf("TC  -  TO %e %e %e\n",tc,to,    tc-to);
             if (do_eos_diag) printf("DEPDT %e \n",depdt);

             if (do_eos_diag) printf("ADDING DD      %e\n",dd);

             // Reset "old" = "current".
             to = tc;
             epo = epc;

             // Create new "current"
             tc = tc + dd;
             temp_row = tc;

             // Update the iteration counter
             icount = icount + 1;

             // If negative go back to the original guess and add 10%
//            if (temp_row(k)<=0.d0) then
             if (temp_row<=min_temp){
                 temp_row = temp0 + 0.1e0*temp0;
                 to = temp_row;
                 findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
                 epo = f[0];
                 tc = MIN(to * 1.01e0, max_temp);
             }
//            if (temp_row(k)<min_temp) then
//                temp_row(k) = min_temp+1.e-6*abs(min_temp)
//                to = temp_row(k)
//                call findthis(table,f,numel,to,rho,ye(k),jy,jq,jr)
//                epo = f(1)
//                if(epo>energy_want(k))then
//                    print *,'energy at min_temp to high', min_temp
//                    print *, 'energy: ', epo, " want: ", energy_want(k)
//                    stop
//                endif
//                tc = MIN(to * 1.01d0, max_temp)
//            endif
             if (temp_row>max_temp){
                 temp_row = max_temp*0.99e0;
                 to = temp_row;
//                temp_row(k) = temp0 - 0.1d0*temp0
//                to = temp_row(k)
                 findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
                 epo = f[0];
//                if(epo<energy_want(k))then
//                    print *,'energy at max_temp to low', max_temp
//                    print *, 'energy: ', epo, " want: ", energy_want(k)
//                    stop
//                endif
                 tc = MIN(to * 0.99e0, max_temp);
             }


             // If the temperature isn't changing much we're done
             if (fabs(dd/temp_row)<newton_tol) break;

             if (temp_row < min_temp || temp_row > max_temp){
                 printf("TEMP OUT OF BOUNDS IN REX ITER %e %e\n", (*temp),temp_row);
                 printf("LIMITS ARE: T1 T2 %e %e\n",min_temp, max_temp);
                 //exit(1);
             }
          }

          (*temp) = temp_row;

          // If we iterated max_newton times and didn't solve:
          if (icount==max_newton && fabs(dd/(*temp))>1.e-10){
             bisection(table,energy_want,1,&temp0,rho,(ye));
             (*temp) = temp0;
          }

          // pressure (MeV/fm^3) -> dyn/cm^2
          (*pres) = f[1] * 1.60217733e33;

          (*gam1) = f[11];

          (*cs) = sqrt((*gam1)*(*pres)/rho);

          //gam1(k) = f(11)
          //cs(k) = f(15)

          (*entropy) = f[2];

    }
    // THIS IS ALSO BASED ON ADAM'S ITERATION SCHEME
    else if (input == eos_input_rp){

//---------------------------------------------------------------------------
// input = 5: dens, pressure, and ye are inputs
//---------------------------------------------------------------------------

       // NOTE: we assume the temp is already in MeV

          // load initial guess
          temp_row = (*temp);

          if (do_eos_diag) printf("T/D INIT %e %e\n",(*temp),(*dens));

          // want to converge to the given energy (erg/g)
          // energy (MeV) per baryon -> ergs/gm
          //   MeV   1.602176462e-6 erg   6.02214199e23
          //   --- x ------------------ x -------------
          //   bar          1 MeV            1 mole

          pressure_want = (*pres) / 1.60217733e33;

          if (do_eos_diag) printf("WANT SCALED PRESSURE %e\n",pressure_want);

          icount = 0;

          // This is the initial guess
          temp0   = (*temp);

          rho = (*dens);

          // Make "old" value = 1.01 * incoming guess
          to = MIN(temp_row * 1.01e0, max_temp);
          findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
          ppo = f[1];

          // Make "current" value = incoming guess
          tc  = temp_row;

          for(i = 1;i<=max_newton;i++){
             findthis(table,f,numel,tc,rho,(ye),&jy,&jq,&jr);
             ppc = f[1];

             if (do_eos_diag){
                printf("ITER %d\n",i);
                printf("TEMP PP OLD %e %e\n",to,ppo);
                printf("TEMP PP NEW %e %e\n",tc,ppc);
             }

             // Calculate the slope of pressure vs temp
             dppdt  = (ppc-ppo)/(tc-to);

             // How far were we off to start with?
             dd = (pressure_want-ppc)/dppdt;

             // Add that much to the current guess.
             if (do_eos_diag) printf("PPC - PPO %e %e %e\n",ppc, ppo, ppc-ppo);
             if (do_eos_diag) printf("TC  -  TO %e %e %e\n",tc,to,    tc-to);
             if (do_eos_diag) printf("DPPDT %e\n",dppdt);

             if (do_eos_diag) printf("ADDING DD    %e\n",dd);

             // Reset "old" = "current".
             to = tc;
             ppo = ppc;

             // Create new "current"
             tc = tc + dd;
             temp_row = tc;

             // Update the iteration counter
             icount = icount + 1;

             // If negative go back to the original guess and add 10%
             if (temp_row<=min_temp){
                 temp_row = temp0 + 0.1e0*temp0;
                 to = temp_row;
                 findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
                 ppo = f[1];
                 tc = MIN(to * 1.01e0, max_temp);
             }
             if (temp_row>max_temp){
                 temp_row = max_temp*0.99e0;
                 to = temp_row;
                 findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
                 ppo = f[1];
                 tc = MIN(to * 0.99e0, max_temp);
             }

             // If the temperature isn't changing much we're done
             if (abs(dd/temp_row)<1.e-8) break;

             if (temp_row < min_temp || temp_row > max_temp){
                 printf("TEMP OUT OF BOUNDS IN REX ITER %e %e\n", (*temp),temp_row);
                 printf("LIMITS ARE: T1 T2 %e %e\n",min_temp, max_temp);
                 //exit(1);
             }

          }

          (*temp) = temp_row;

          // If we iterated max_newton times and didn't solve:
          if (icount==max_newton){
             //print *, 'no bisection routine for pressure solve in EOS'
             //stop
             bisection(table,pressure_want,2,&temp0,rho,(ye));
             findthis(table,f,numel,temp0,rho,(ye),&jy,&jq,&jr);
             (*temp) = temp0;
          }

          // energy
          (*eint) = (f[0]+const_hack) * 0.95655684e18;
          if(do_eos_diag) printf("Energy returned is %e\n", (*eint)*rho);
    }
    else if (input == eos_input_rs){

       // NOTE: we assume the temp is already in MeV

          // load initial guess
          temp_row = (*temp);

          if (do_eos_diag) printf("T/D INIT %e %e\n",(*temp),(*dens));

          // want to converge to the given energy (erg/g)
          // energy (MeV) per baryon -> ergs/gm
          //   MeV   1.602176462e-6 erg   6.02214199e23
          //   --- x ------------------ x -------------
          //   bar          1 MeV            1 mole

// jbb
          entropy_want = (*entropy);

          if (do_eos_diag) printf("WANT S  %e\n",entropy_want);

          icount = 0;

          // This is the initial guess
          temp0   = (*temp);

          findthis(table,f,numel,max_temp,(*dens),(ye),&jy,&jq,&jr);
          if (entropy_want > f[2]){
            printf("eos_given_REX: entropy is too high given this rho and max_temp: %e\n",entropy_want);
            printf("maximum possible entropy given this density                     %e\n",f[2]);
            printf("density and max_temp and Ye: %e %e %e",(*dens), max_temp, (ye));
            //exit(1);
          }

          rho = (*dens);

          // Make "old" value = 1.01 * incoming guess
          to = MIN(temp_row * 1.01e0, max_temp);
          findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
          epo = f[2];

          // Make "current" value = incoming guess
          tc  = temp_row;

          for(i = 1;i<=max_newton;i++){

             findthis(table,f,numel,tc,rho,(ye),&jy,&jq,&jr);
             epc = f[2];

             if (do_eos_diag){
                printf("ITER %d\n",i);
                printf("TEMP EP OLD %e %e\n",to,epo);
                printf("TEMP EP NEW %e %e\n",tc,epc);
             }

             // Calculate the slope of energy vs temp
             depdt  = (epc-epo)/(tc-to);

             // How far were we off to start with?
             dd = (entropy_want-epc)/depdt;

             // Add that much to the current guess.
             if (do_eos_diag) printf("EPC - EPO %e %e %e\n",epc, epo, epc-epo);
             if (do_eos_diag) printf("TC  -  TO %e %e %e\n",tc,to,    tc-to);
             if (do_eos_diag) printf("DEPDT %e\n",depdt);

             if (do_eos_diag) printf("ADDING DD      %e\n",dd);

             // Reset "old" = "current".
             to = tc;
             epo = epc;

             // Create new "current"
             tc = tc + dd;
             temp_row = tc;

             // Update the iteration counter
             icount = icount + 1;

             // If negative go back to the original guess and add 10%
             if (temp_row<=0.e0){
                 temp_row = temp0 + 0.1e0*temp0;
                 to = temp_row;
                 findthis(table,f,numel,to,rho,(ye),&jy,&jq,&jr);
                 epo = f[2];
                 tc = MIN(to * 1.01e0, max_temp);
             }

             // If the temperature isn't changing much we're done
             if (abs(dd/temp_row)<1.e-8) break;

             if (temp_row < min_temp || (*temp) > max_temp){
                 printf("TEMP OUT OF BOUNDS IN REX ITER %e\n", (*temp));
                 printf("LIMITS ARE: T1 T2 %e %e\n",min_temp, max_temp);
                 //exit(1);
             }

          }
          (*temp) = temp_row;

          // If we iterated max_newton times and didn't solve:
          if (icount==max_newton){
             bisection(table,entropy_want,3,&temp0,rho,(ye));
             findthis(table,f,numel,temp0,rho,(ye),&jy,&jq,&jr);
             (*temp) = temp0;
          }

          // pressure (MeV/fm^3) -> dyn/cm^2
          (*pres) = f[2-1] * 1.60217733e33;

          (*gam1) = f[11];

            (*cs) = sqrt((*gam1)*(*pres)/rho);

//         eint(k) = f(1) * 0.95655684d18
// jbb
          (*eint) = (f[0]+const_hack) * 0.95655684e18;

    }
    else if (input == eos_input_cv) {


          findthis(table,f,numel,(*temp),(*dens),(ye),&jy,&jq,&jr);

          (*cv) = f[3];

    }
    else if (input == eos_input_etaxnxp){

             findthis(table,f,numel,(*temp),(*dens),(ye),&jy,&jq,&jr);
             (*cs) = f[14];
             (*xn) = f[4];
             (*xp) = f[5];
    }
    else {
       printf("Not implemented conversion in EOS\n");
    }
}

void bisection(const double *table,double target,int nn,double *temp,double rho, double ye){
  DEFINE_EOSCONST_1
  DEFINE_EOSCONST_2
      double ep,e,ep1,ep2,dt,depdt,dd;
      double temp1,temp2,dtemp,tmid,temp00;
      double f[numel],f1[numel],f2[numel],fmid[numel];
      double min_temp,max_temp,eps_temp,eps_ener;
      int i,jy,jq,jr,mcount;

      min_temp = pow(10.e0,t1);
      max_temp = pow(10.e0,t2);

      eps_temp = 1.e-12;

      // Compute the energy of max_temp to be used in defining whether we're "close enough"
      findthis(table,f1,numel,max_temp,rho,ye,&jy,&jq,&jr);
      eps_ener = 1.e-12 * f1[nn-1];

      // target is the target value for energy
      // nn is the index in f, i.e. energy = f(1)
      // temp is what we're solving for
      // rho and ye are fixed values of density and Ye

      temp00 = *temp;
      mcount = 0;

//     temp1 = 0.9d0 * temp
//     temp2 = 1.1d0 * temp
      temp1 = MAX(0.8e0 * (*temp),min_temp);
      temp2 = MIN(1.2e0 * (*temp),max_temp);

      findthis(table,f,numel,*temp,rho,ye,&jy,&jq,&jr);
      f[nn-1]=f[nn-1]-target;

      findthis(table,f1,numel,temp1,rho,ye,&jy,&jq,&jr);
      f1[nn-1]=f1[nn-1]-target;

      findthis(table,f2,numel,temp2,rho,ye,&jy,&jq,&jr);
      f2[nn-1]=f2[nn-1]-target;


      while(1){
        if (f1[nn-1]*f2[nn-1]>=0.e0){
           mcount=mcount+1;
           temp1 = MAX(0.8e0 * temp1,min_temp);
           temp2 = MIN(1.2e0 * temp2,max_temp);
           findthis(table,f1,numel,temp1,rho,ye,&jy,&jq,&jr);
           findthis(table,f2,numel,temp2,rho,ye,&jy,&jq,&jr);
           f1[nn-1]=f1[nn-1]-target;
           f2[nn-1]=f2[nn-1]-target;
           if (fabs(f1[nn-1]) <= eps_ener){
              *temp = temp1;
              return;
           }
           if (abs(f2[nn-1]) <= eps_ener){
              *temp = temp2;
              return;
           }
           if (mcount>50){
              findthis(table,f1,numel,min_temp,rho,ye,&jy,&jq,&jr);
              findthis(table,f1,numel,max_temp,rho,ye,&jy,&jq,&jr);
              printf( "BISECTION FAILED in eos_stuff\n");
              printf( "MINTEMP / MAXTEMP %e %e\n",min_temp,max_temp);
              printf( "  TEMP1 / TEMP2   %e %e\n",temp1, temp2);
              printf( "TARGET ENERGY     %e\n",target);
              printf( "ENERGY OF MINTEMP %e\n",f1[nn-1]);
              printf( "ENERGY OF MAXTEMP %e\n",f1[nn-1]);
              printf( "eos_give_ReX: bisection cant bracket energy %e %e %e\n", rho, ye, *temp);
           }
        }
      }
      if (f1[nn-1]<0.e0){
         *temp=temp1;
         dtemp=temp2-temp1;
      }
      else{
         *temp=temp2;
         dtemp=temp1-temp2;
      }

      for(i=1;i<=80;i++){
         dtemp = dtemp*0.5e0;
          tmid = *temp+dtemp;
         findthis(table,fmid,numel,tmid,rho,ye,&jy,&jq,&jr);
         fmid[nn-1] = fmid[nn-1] - target;
         if (fmid[nn-1]<=0.e0) *temp = tmid;
         if (fabs(dtemp)<eps_temp) return;
      }

      printf( "eos_given_ReX: bisection cant get to the energy we want after 50 iterations\n");
}
