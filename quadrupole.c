#include <math.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>

#include "constants.h"
#include "decs.h"

#if USE_MPI==TRUE
static MPI_Request qpole_mpi_req;
#endif

#if NDIM==3
#define QPOLE_NCOMP (6)
#else
#define QPOLE_NCOMP (1)
#endif
#define QPOLE_NELEM (2*QPOLE_NCOMP)
static double qpole_loc[QPOLE_NELEM];
static double qpole[QPOLE_NELEM];

void quadrupole_init() {
    #if NDIM == 1
        return;
    #endif

    int ierr;
    ierr = array_init(sizeof(t),     &qpole_time); assert(!ierr);
    ierr = array_init(sizeof(qpole), &qpole_mom); assert(!ierr);
}

void quadrupole_free() {
    #if NDIM == 1
        return;
    #endif

    array_free(qpole_time); qpole_time = NULL;
    array_free(qpole_mom);  qpole_mom  = NULL;
}

void quadrupole_start() {
    #if NDIM==1
        return;
    #endif

    TIMER_START("quadrupole_start");

    int ii=0, jj=0, kk=0, dd=0;
    int ierr;

    #if NDIM == 2
        double I_zz=0, I_zz_dot=0;
        ZLOOP {
            int jlev = 0;
            int s = DJS(ii);
            while (s >>= 1) jlev++;

            double const r = r_of_x(rx_info,beta0[ii]/Gamma0[ii]);;
            double const th = th_of_x(thx_info,beta1s[jlev][jj]/Gamma1s[jlev][jj]);

            double const cos_th = cos(th);
            double const sin_th = sin(th);
            double const P2p    = 0.5*(3.*SQR(cos_th) - 1.);
            double const dP2dth = -3.*cos_th*sin_th;

            double const vol = ND_ELEM_LINEAR(geom,ii,jj,kk).volume;

            double const rho = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
            double const sr  = ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0];
            double const vr  = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1)*sr;
            double const sth = ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][1];
            double const vth = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U2)*sth;

            I_zz     += 2./3.*(rho/MSUN)*SQR(r)*P2p*vol;
            I_zz_dot += 4./3.*(rho/MSUN)*r*(P2p*vr + 0.5*dP2dth*vth)*vol;
        }

        qpole_loc[0] = I_zz;
        qpole_loc[1] = I_zz_dot;
    #elif NDIM == 3
        double I_xx=0, I_xy=0, I_xz=0, I_yy=0, I_yz=0, I_zz=0;
        double Q_xx=0, Q_xy=0, Q_xz=0, Q_yy=0, Q_yz=0, Q_zz=0;
        ZLOOP {
            int jlev = 0;
            int s = DJS(ii);
            while (s >>= 1) jlev++;

            double const r  = r_of_x(rx_info,beta0[ii]/Gamma0[ii]);;
            double const th = th_of_x(thx_info,beta1s[jlev][jj]/Gamma1s[jlev][jj]);
            double const ph = startx[2] + (kk+0.5)*DKS(ii,jj)*dx[2];

            double const sth = sin(th);
            double const cth = cos(th);
            double const sph = sin(ph);
            double const cph = cos(ph);

            double const xp[3] = {
                r*sth*cph, r*sth*sph, r*cth
            };
            double const np[3] = {
                sth*cph, sth*sph, cth
            };

            double lam[3][3], v[3];
            lam[0][0] = sth*cph;
            lam[0][1] = cth*cph;
            lam[0][2] = -sph;
            lam[1][0] = sth*sph;
            lam[1][1] = cth*sph;
            lam[1][2] = cph;
            lam[2][0] = cth;
            lam[2][1] = -sth;
            lam[2][2] = 0.0;
            for (int l=0; l<3; l++) {
                v[l] = 0.0;
                for (int m=0; m<3; m++) {
                    v[l] += lam[l][m]*NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1+m)*ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][m];
                }
            }

            double const vol = ND_ELEM_LINEAR(geom,ii,jj,kk).volume;
            double const rho = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,RHO);
            double const sr  = ND_ELEM_LINEAR(geom,ii,jj,kk).scale[0][0];
            double const vr  = NDP_ELEM_LINEAR(sim_p,ii,jj,kk,U1)*sr;

            double * I_tensor[3][3] = {
                {&I_xx, &I_xy, &I_xz},
                {&I_xy, &I_yy, &I_yz},
                {&I_xz, &I_yz, &I_zz},
            };
            double * Q_tensor[3][3] = {
                {&Q_xx, &Q_xy, &Q_xz},
                {&Q_xy, &Q_yy, &Q_yz},
                {&Q_xz, &Q_yz, &Q_zz},
            };

            for(int a = 0; a < 3; ++a)
            for(int b = a; b < 3; ++b) {
                *I_tensor[a][b] += vol*(rho/MSUN)*SQR(r)*
                    (np[a]*np[b] - 1.0/3.0*DELTA(a,b));
                *Q_tensor[a][b] += vol*(rho/MSUN)*
                    (v[a]*xp[b] + v[b]*xp[a] - 2.0/3.0*DELTA(a,b)*r*vr);
            }
        }

        qpole_loc[ 0] = I_xx;
        qpole_loc[ 1] = I_xy;
        qpole_loc[ 2] = I_xz;
        qpole_loc[ 3] = I_yy;
        qpole_loc[ 4] = I_yz;
        qpole_loc[ 5] = I_zz;
        qpole_loc[ 6] = Q_xx;
        qpole_loc[ 7] = Q_xy;
        qpole_loc[ 8] = Q_xz;
        qpole_loc[ 9] = Q_yy;
        qpole_loc[10] = Q_yz;
        qpole_loc[11] = Q_zz;
    #endif // NDIM switch

    #if USE_MPI==TRUE
        ierr = MPI_Ireduce(&qpole_loc[0], &qpole[0], QPOLE_NELEM, MPI_DOUBLE,
                MPI_SUM, 0, MPI_COMM_WORLD, &qpole_mpi_req); assert(!ierr);
    #else
        memcpy(&qpole[0], &qpole_loc[0], QPOLE_NELEM*sizeof(double));
    #endif // MPI

    TIMER_STOP;
    return;
}

void quadrupole_dump() {
    int ierr;
    #if NDIM == 1
        return;
    #endif
    TIMER_START("quadrupole_dump");
    if(mpi_io_proc()) {
        assert(qpole_time->size == qpole_mom->size);
        if(qpole_time->size == 0) {
            TIMER_STOP;
            return;
        }

        char const * quadrupole_fname = "quadrupole.dat";
        FILE * fp = NULL;
        if(0 != access(quadrupole_fname, F_OK)) {
            fprintf(stderr, "Creating %s\n", quadrupole_fname);
            fp = fopen(quadrupole_fname, "w");
            #if NDIM == 2
            fprintf(fp, "# 1:time 2:I_zz 3:I_zz_dot\n");
            #else
            fprintf(fp, "# 1:time 2:I_xx 3:I_xy 4:I_xz 5:I_yy 6:I_yz 7:I_zz "
                    "8:I_dot_xx 9:I_dot_xy 10:I_dot_xz 11:I_dot_yy 12:I_dot_yz "
                    "13:I_dot_zz\n");
            #endif
        }
        else {
            fprintf(stderr, "Appending to %s\n", quadrupole_fname);
            fp = fopen(quadrupole_fname, "a");
        }

        for(int idx = 0; idx < qpole_time->size; ++idx) {
            double * p_t;
            ierr = array_at(qpole_time, idx, (void **)&p_t); assert(!ierr);
            fprintf(fp, "%16.15e ", *p_t);

            double * p_qpole;
            ierr = array_at(qpole_mom, idx, (void **)&p_qpole); assert(!ierr);
            for(int el = 0; el < QPOLE_NELEM-1; ++el) {
                fprintf(fp, "%16.15e ", p_qpole[el]);
            }
            fprintf(fp, "%16.15e\n", p_qpole[QPOLE_NELEM-1]);
        }
        ierr = array_clear(qpole_time); assert(!ierr);
        ierr = array_clear(qpole_mom); assert(!ierr);
        fclose(fp);
    }
    TIMER_STOP;
    return;
}

void quadrupole_end() {
    int ierr;
    #if NDIM == 1
        return;
    #endif
    TIMER_START("quadrupole_end");
    #if USE_MPI==TRUE
        ierr = MPI_Wait(&qpole_mpi_req, MPI_STATUS_IGNORE); assert(!ierr);
    #endif
    if(mpi_io_proc()) {
        ierr = array_push_back(qpole_time, &t); assert(!ierr);
        ierr = array_push_back(qpole_mom, &qpole[0]); assert(!ierr);
    }
    TIMER_STOP;
    return;
}
