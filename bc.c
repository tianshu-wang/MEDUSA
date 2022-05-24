
#include "decs.h"
#include "constants.h"

#if (USE_LINEAR_ALLOCATION==TRUE)
void reset_boundaries(double ** p,int on_gpu)
#else
void reset_boundaries(double NDP_PTR p,int on_gpu)
#endif
{
  int i=0,j=0,k=0,vv,dd,ip,jp,kp,njp,nkp,jj,kk;
  int jstart,jstop,kstart,kstop;
  double r,r0,vol_sum,var_sum;

  TIMER_START("reset_boundaries");

  // communicate boundary data between MPI domains
  if (p == sim_p) half_step_sync = 0;
  else            half_step_sync = 1;


  if (half_step_sync) t += dt;

  // set the physical boundary conditions
  // Inner i-boundary
  if (istart[0] == 0) {

    VLOOP {

      switch (bc[vv].lo[0]) {

        case OUTFLOW:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istart[0],j,k,vv);
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case REFLECT:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                ip = istart[0] - 1 - i;
                jp = j*DJS(i)/DJS(ip);
                kp = k*DKS(i,j)/DKS(ip,jp);
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,ip,jp,kp,vv);
                if (vv==U1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==0) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        #if (GEOM==SPHERICAL)
        case SPHERICAL_ORIGIN:
          if (istart[2] == 0 && istop[2] == global_grid_dims[2]) {  // only do this if self-periodic
            for (i=istart[0]-NG; i<istart[0]; i++) {
              ip  = istart[0] - 1 - i;  // i reflected across origin
              #if (NDIM>1)
              JSLOOP(i,j) {
                njp = my_grid_dims[1]/DJS(ip);
                jp  = njp - 1 - (j - JS(ip,istart[1])) + JS(ip,istart[1]);  // j reflected across origin
              #endif
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                  nkp = my_grid_dims[2]/DKS(ip,jp);
                  kp  = (k - KS(ip,jp,istart[2]) + nkp/2) % nkp + KS(ip,jp,istart[2]);  // k periodic across origin
                #endif
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,ip,jp,kp,vv);
                  if (vv==U1 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #if (DO_RADIATION)
                  if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==0 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #endif
                #if (NDIM==3)
                }
                #endif
              #if (NDIM>1)
              }
              #endif
            }
          }
          break;
        #endif /* GEOM==SPHERICAL && NDIM>1 */

        case EXTRAP:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istart[0],j,k,vv) + (i-istart[0])*(NDP_ELEM_LINEAR(p,istart[0]+1,j,k,vv)-NDP_ELEM_LINEAR(p,istart[0],j,k,vv));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;


        case PROB:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case PERIODIC:
          if (istop[0] == global_grid_dims[0]) {
            for (i=istart[0]-NG; i<istart[0]; i++) {
              #if (NDIM>1)
              JSLOOP(i,j) {
              #endif
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i+my_grid_dims[0],j,k,vv);
                #if (NDIM==3)
                }
                #endif
              #if (NDIM>1)
              }
              #endif
            }
          }
          break;

        case DISK:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) =  NDP_ELEM_LINEAR(p,0,j,k,vv) + i*(NDP_ELEM_LINEAR(p,1,j,k,vv)-NDP_ELEM_LINEAR(p,0,j,k,vv));
                if (vv == U1 && NDP_ELEM_LINEAR(p,i,j,k,vv) > 0.0) NDP_ELEM_LINEAR(p,i,j,k,vv) = 0.0;
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].lo[0] = %d\n", vv, bc[vv].lo[0]);
          //exit(1);
      }
    }
  }

  // Outer i-boundary:  must assume there is no j- or k-refinement with i here
  if (istop[0] == global_grid_dims[0]) {

    VLOOP {

      switch (bc[vv].hi[0]) {

        case OUTFLOW:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv);
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case REFLECT:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,2*istop[0]-i-1,j,k,vv);
                if (vv==U1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==0) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case EXTRAP:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv) + (i-istop[0]+1)*(NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv)-NDP_ELEM_LINEAR(p,istop[0]-2,j,k,vv));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        #if (GEOM==SPHERICAL)
        case RADEXTRAP:
          r0 = r_of_x(rx_info,(istop[0]-0.5)*dx[0]+startx[0]);
          for (i=istop[0]; i<istop[0]+NG; i++) {
            r = r_of_x(rx_info,startx[0] + (i+0.5)*dx[0]);
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                if (vv>=irad1 && vv < ifrad1) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = pow(r0/r,2.0)*NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv);
                } else if (vv >= ifrad1 && vv < nvars) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = CLIGHT*NDP_ELEM_LINEAR(p,i,j,k,irad1 + (vv-ifrad1)/NDIM)/ND_ELEM_LINEAR(geom,i,j,k).scale[0][(vv-ifrad1)%NDIM];
                } else {
                  printf("SET RADEXTRAP boundary condition for non-rad variable\n");
                  //exit(1);
                }
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;
        #endif

        case PROB:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case PERIODIC:
          if (istart[0] == 0) {
            for (i=istop[0]; i<istop[0]+NG; i++) {
              #if (NDIM>1)
              JSLOOP(i,j) {
              #endif
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i-my_grid_dims[0],j,k,vv);
                #if (NDIM==3)
                }
                #endif
              #if (NDIM>1)
              }
              #endif
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].hi[0] = %d\n", vv, bc[vv].hi[0]);
          //exit(1);
      }
    }
  }


  #if (NDIM>1)
  // Inner j-boundary
  if (istart[1] == 0) {

    VLOOP {

      switch (bc[vv].lo[1]) {

        case OUTFLOW:
          ISLOOP(i) {
            jstart = JS(i,istart[1]);
            for (j=jstart-NG; j<jstart; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jstart,k,vv);
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case REFLECT:
          ISLOOP(i) {
            jstart = JS(i,istart[1]);
            for (j=jstart-NG; j<jstart; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                jp = jstart - 1 - j;
                kp = k*DKS(i,j)/DKS(i,jp);  // probably unnecessary; ought to be that DKS(i,jp)=DKS(i,j)
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                if (vv==U2) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case PERIODIC:
          if (istop[1] == global_grid_dims[1]) {
            ISLOOP(i) {
              jstart = JS(i,istart[1]);
              for (j=jstart-NG; j<jstart; j++) {
                njp = my_grid_dims[1]/DJS(i);
                jp = j + njp;
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  kp = KS(i,jp,k);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;

        #if (GEOM==SPHERICAL)
        case SPHERICAL_ORIGIN:
          if (istart[2] == 0 && istop[2] == global_grid_dims[2]) {  // only do this if self-periodic
            ISLOOP(i) {
              jstart = JS(i,istart[1]);
              for (j=jstart-NG; j<jstart; j++) {
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  jp = 2*jstart - 1 - j;
                  nkp = my_grid_dims[2]/DKS(i,jp);
                  kp = (k - KS(i,jp,istart[2]) + nkp/2) % nkp + KS(i,jp,istart[2]);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                  if (vv==U2 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #if (DO_RADIATION)
                  if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==1 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #endif
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;
        #endif /* GEOM==SPHERICAL && NDIM==3 */

        case PROB:
          ISLOOP(i) {
            jstart = JS(i,istart[1]);
            for (j=jstart-NG; j<jstart; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;
            
        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].lo[1] = %d\n", vv, bc[vv].lo[1]);
          //exit(1);
      }
    }
  }

  // Outer j-boundary
  if (istop[1] == global_grid_dims[1]) {

    VLOOP {

      switch (bc[vv].hi[1]) {

        case OUTFLOW:
          ISLOOP(i) {
            jstop = JS(i,istop[1]);
            for (j=jstop; j<jstop+NG; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jstop-1,k,vv);
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case REFLECT:
          ISLOOP(i) {
            jstop = JS(i,istop[1]);
            for (j=jstop; j<jstop+NG; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                jp = 2*jstop - 1 - j;
                kp = k*DKS(i,j)/DKS(i,jp);
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                if (vv==U2) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case PERIODIC:
          if (istart[1] == 0) {  // This proc has ALL zones in 1-direction
            ISLOOP(i) {
              jstop = JS(i,istop[1]);
              for (j=jstop; j<jstop+NG; j++) {
                njp = my_grid_dims[1]/DJS(i);
                jp = j - njp;
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  kp = KS(i,jp,k);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;

        #if (GEOM==SPHERICAL)
        case SPHERICAL_ORIGIN:
        if (istart[2] == 0 && istop[2] == global_grid_dims[2]) {  // only do this if self-periodic
            ISLOOP(i) {
              jstop = JS(i,istop[1]);
              for (j=jstop; j<jstop+NG; j++) {
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  jp = 2*jstop - 1 - j;
                  nkp = my_grid_dims[2]/DKS(i,jp);
                  kp = (k - KS(i,jp,istart[2]) + nkp/2) % nkp + KS(i,jp,istart[2]);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                  if (vv==U2 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #if (DO_RADIATION)
                  if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==1 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #endif
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;
        #endif /* GEOM==SPHERICAL && NDIM==3 */

        case PROB:
          ISLOOP(i) {
            jstop = JS(i,istop[1]);
            for (j=jstop; j<jstop+NG; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].hi[1] = %d\n", vv, bc[vv].hi[1]);
          //exit(1);
      }
    }
  }
  #endif


  #if (NDIM==3)
  // Inner k-boundary
  if (istart[2] == 0) {

    VLOOP {

      switch (bc[vv].lo[2]) {

        case OUTFLOW:
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstart = KS(i,j,istart[2]);
              for (k=kstart-NG; k<kstart; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstart,vv);
              }
            }
          }
          break;

        case REFLECT:
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstart = KS(i,j,istart[2]); 
              for(k=kstart-NG; k<kstart; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstart,vv);
                if (vv == U3) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= (-1.);
                }
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==2) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                }
                #endif
              }
            }
          }
          break;

        case PERIODIC:
          if (istop[2] == global_grid_dims[2]) {
            // Only do this if proc is periodic with itself; otherwise it's an MPI boundary
            ISLOOP(i) {
              JSLOOP(i,j) {
                kstart = KS(i,j,istart[2]);
                for (k=kstart-NG; k<kstart; k++) {
                  nkp = my_grid_dims[2]/DKS(i,j);
                  kp = k + nkp;
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kp,vv);
                }
              }
            }
          }
          break;

        case PROB:
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstart = KS(i,j,istart[2]);
              for (k=kstart-NG; k<kstart; k++) {
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              }
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].lo[2] = %d\n", vv, bc[vv].lo[2]);
          //exit(1);
      }
    }
  }

  // Outer k-boundary
  if (istop[2] == global_grid_dims[2]) {
    
    VLOOP {

      switch (bc[vv].hi[2]) {

        case OUTFLOW:
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstop = KS(i,j,istop[2]);
              for (k=kstop; k<kstop+NG; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstop-1,vv);
              }
            }
          }
          break;

        case REFLECT:
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstop = KS(i,j,istop[2]);
              for(k=kstop; k<kstop+NG; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstop-1,vv);
                if (vv == U3) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= (-1.);
                }
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==2) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                }
                #endif
              }
            }
          }
          break;

        case PERIODIC:
          if (istart[2] == 0) {
            // Only do this if proc is periodic with itself; otherwise it's an MPI boundary
            ISLOOP(i) {
              JSLOOP(i,j) {
                kstop = KS(i,j,istop[2]);
                for (k=kstop; k<kstop+NG; k++) {
                  nkp = my_grid_dims[2]/DKS(i,j);
                  kp = k - nkp;
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kp,vv);
                }
              }
            }
          }
          break;

        case PROB:
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstop = KS(i,j,istop[2]);
              for (k=kstop; k<kstop+NG; k++) {
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              }
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].hi[2] = %d\n", vv, bc[vv].hi[2]);
          //exit(1);
      }
    }
  }
  #endif

  if (half_step_sync) t -= dt;

  TIMER_STOP;
  return;
}

#if (USE_LINEAR_ALLOCATION==TRUE)
void reset_boundaries_gpu(double ** p,int on_gpu)
#else
void reset_boundaries_gpu(double NDP_PTR p,int on_gpu)
#endif
{
  int i=0,j=0,k=0,vv,dd,ip,jp,kp,njp,nkp,jj,kk;
  int jstart,jstop,kstart,kstop;
  double r,r0,vol_sum,var_sum;

  TIMER_START("reset_boundaries");

  // communicate boundary data between MPI domains
  if (p == sim_p) half_step_sync = 0;
  else            half_step_sync = 1;


  if (half_step_sync) t += dt;

  // set the physical boundary conditions
  // Inner i-boundary
  if (istart[0] == 0) {

    GPU_PRAGMA(omp target teams distribute parallel for)
    VLOOP {

      switch (bc[vv].lo[0]) {

        case OUTFLOW:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istart[0],j,k,vv);
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case REFLECT:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                ip = istart[0] - 1 - i;
                jp = j*DJS(i)/DJS(ip);
                kp = k*DKS(i,j)/DKS(ip,jp);
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,ip,jp,kp,vv);
                if (vv==U1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==0) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        #if (GEOM==SPHERICAL)
        case SPHERICAL_ORIGIN:
          if (istart[2] == 0 && istop[2] == global_grid_dims[2]) {  // only do this if self-periodic
            for (i=istart[0]-NG; i<istart[0]; i++) {
              ip  = istart[0] - 1 - i;  // i reflected across origin
              #if (NDIM>1)
              JSLOOP(i,j) {
                njp = my_grid_dims[1]/DJS(ip);
                jp  = njp - 1 - (j - JS(ip,istart[1])) + JS(ip,istart[1]);  // j reflected across origin
              #endif
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                  nkp = my_grid_dims[2]/DKS(ip,jp);
                  kp  = (k - KS(ip,jp,istart[2]) + nkp/2) % nkp + KS(ip,jp,istart[2]);  // k periodic across origin
                #endif
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,ip,jp,kp,vv);
                  if (vv==U1 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #if (DO_RADIATION)
                  if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==0 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #endif
                #if (NDIM==3)
                }
                #endif
              #if (NDIM>1)
              }
              #endif
            }
          }
          break;
        #endif /* GEOM==SPHERICAL && NDIM>1 */

        case EXTRAP:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istart[0],j,k,vv) + (i-istart[0])*(NDP_ELEM_LINEAR(p,istart[0]+1,j,k,vv)-NDP_ELEM_LINEAR(p,istart[0],j,k,vv));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;


        case PROB:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case PERIODIC:
          if (istop[0] == global_grid_dims[0]) {
            for (i=istart[0]-NG; i<istart[0]; i++) {
              #if (NDIM>1)
              JSLOOP(i,j) {
              #endif
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i+my_grid_dims[0],j,k,vv);
                #if (NDIM==3)
                }
                #endif
              #if (NDIM>1)
              }
              #endif
            }
          }
          break;

        case DISK:
          for (i=istart[0]-NG; i<istart[0]; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) =  NDP_ELEM_LINEAR(p,0,j,k,vv) + i*(NDP_ELEM_LINEAR(p,1,j,k,vv)-NDP_ELEM_LINEAR(p,0,j,k,vv));
                if (vv == U1 && NDP_ELEM_LINEAR(p,i,j,k,vv) > 0.0) NDP_ELEM_LINEAR(p,i,j,k,vv) = 0.0;
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].lo[0] = %d\n", vv, bc[vv].lo[0]);
          //exit(1);
      }
    }
  }

  // Outer i-boundary:  must assume there is no j- or k-refinement with i here
  if (istop[0] == global_grid_dims[0]) {

    GPU_PRAGMA(omp target teams distribute parallel for)
    VLOOP {

      switch (bc[vv].hi[0]) {

        case OUTFLOW:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv);
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case REFLECT:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,2*istop[0]-i-1,j,k,vv);
                if (vv==U1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==0) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case EXTRAP:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv) + (i-istop[0]+1)*(NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv)-NDP_ELEM_LINEAR(p,istop[0]-2,j,k,vv));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        #if (GEOM==SPHERICAL)
        case RADEXTRAP:
          r0 = r_of_x(rx_info,(istop[0]-0.5)*dx[0]+startx[0]);
          for (i=istop[0]; i<istop[0]+NG; i++) {
            r = r_of_x(rx_info,startx[0] + (i+0.5)*dx[0]);
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                if (vv>=irad1 && vv < ifrad1) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = pow(r0/r,2.0)*NDP_ELEM_LINEAR(p,istop[0]-1,j,k,vv);
                } else if (vv >= ifrad1 && vv < nvars) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = CLIGHT*NDP_ELEM_LINEAR(p,i,j,k,irad1 + (vv-ifrad1)/NDIM)/ND_ELEM_LINEAR(geom,i,j,k).scale[0][(vv-ifrad1)%NDIM];
                } else {
                  printf("SET RADEXTRAP boundary condition for non-rad variable\n");
                  //exit(1);
                }
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;
        #endif

        case PROB:
          for (i=istop[0]; i<istop[0]+NG; i++) {
            #if (NDIM>1)
            JSLOOP(i,j) {
            #endif
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            #if (NDIM>1)
            }
            #endif
          }
          break;

        case PERIODIC:
          if (istart[0] == 0) {
            for (i=istop[0]; i<istop[0]+NG; i++) {
              #if (NDIM>1)
              JSLOOP(i,j) {
              #endif
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i-my_grid_dims[0],j,k,vv);
                #if (NDIM==3)
                }
                #endif
              #if (NDIM>1)
              }
              #endif
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].hi[0] = %d\n", vv, bc[vv].hi[0]);
          //exit(1);
      }
    }
  }


  #if (NDIM>1)
  // Inner j-boundary
  if (istart[1] == 0) {

    VLOOP {

      switch (bc[vv].lo[1]) {

        case OUTFLOW:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            jstart = JS(i,istart[1]);
            for (j=jstart-NG; j<jstart; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jstart,k,vv);
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case REFLECT:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            jstart = JS(i,istart[1]);
            for (j=jstart-NG; j<jstart; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                jp = jstart - 1 - j;
                kp = k*DKS(i,j)/DKS(i,jp);  // probably unnecessary; ought to be that DKS(i,jp)=DKS(i,j)
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                if (vv==U2) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case PERIODIC:
          if (istop[1] == global_grid_dims[1]) {
    GPU_PRAGMA(omp target teams distribute parallel for)
            ISLOOP(i) {
              jstart = JS(i,istart[1]);
              for (j=jstart-NG; j<jstart; j++) {
                njp = my_grid_dims[1]/DJS(i);
                jp = j + njp;
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  kp = KS(i,jp,k);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;

        #if (GEOM==SPHERICAL)
        case SPHERICAL_ORIGIN:
          if (istart[2] == 0 && istop[2] == global_grid_dims[2]) {  // only do this if self-periodic
    GPU_PRAGMA(omp target teams distribute parallel for)
            ISLOOP(i) {
              jstart = JS(i,istart[1]);
              for (j=jstart-NG; j<jstart; j++) {
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  jp = 2*jstart - 1 - j;
                  nkp = my_grid_dims[2]/DKS(i,jp);
                  kp = (k - KS(i,jp,istart[2]) + nkp/2) % nkp + KS(i,jp,istart[2]);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                  if (vv==U2 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #if (DO_RADIATION)
                  if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==1 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #endif
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;
        #endif /* GEOM==SPHERICAL && NDIM==3 */

        case PROB:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            jstart = JS(i,istart[1]);
            for (j=jstart-NG; j<jstart; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;
            
        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].lo[1] = %d\n", vv, bc[vv].lo[1]);
          //exit(1);
      }
    }
  }

  // Outer j-boundary
  if (istop[1] == global_grid_dims[1]) {

    VLOOP {

      switch (bc[vv].hi[1]) {

        case OUTFLOW:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            jstop = JS(i,istop[1]);
            for (j=jstop; j<jstop+NG; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jstop-1,k,vv);
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case REFLECT:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            jstop = JS(i,istop[1]);
            for (j=jstop; j<jstop+NG; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                jp = 2*jstop - 1 - j;
                kp = k*DKS(i,j)/DKS(i,jp);
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                if (vv==U2) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==1) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                #endif
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        case PERIODIC:
          if (istart[1] == 0) {  // This proc has ALL zones in 1-direction
    GPU_PRAGMA(omp target teams distribute parallel for)
            ISLOOP(i) {
              jstop = JS(i,istop[1]);
              for (j=jstop; j<jstop+NG; j++) {
                njp = my_grid_dims[1]/DJS(i);
                jp = j - njp;
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  kp = KS(i,jp,k);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;

        #if (GEOM==SPHERICAL)
        case SPHERICAL_ORIGIN:
        if (istart[2] == 0 && istop[2] == global_grid_dims[2]) {  // only do this if self-periodic
    GPU_PRAGMA(omp target teams distribute parallel for)
            ISLOOP(i) {
              jstop = JS(i,istop[1]);
              for (j=jstop; j<jstop+NG; j++) {
                #if (NDIM==3)
                KSLOOP(i,j,k) {
                #endif
                  jp = 2*jstop - 1 - j;
                  nkp = my_grid_dims[2]/DKS(i,jp);
                  kp = (k - KS(i,jp,istart[2]) + nkp/2) % nkp + KS(i,jp,istart[2]);
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,jp,kp,vv);
                  if (vv==U2 || vv==U3) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #if (DO_RADIATION)
                  if (vv>=ifrad1 && ((vv-ifrad1)%NDIM==1 || (vv-ifrad1)%NDIM==2)) NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                  #endif
                #if (NDIM==3)
                }
                #endif
              }
            }
          }
          break;
        #endif /* GEOM==SPHERICAL && NDIM==3 */

        case PROB:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            jstop = JS(i,istop[1]);
            for (j=jstop; j<jstop+NG; j++) {
              #if (NDIM==3)
              KSLOOP(i,j,k) {
              #endif
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              #if (NDIM==3)
              }
              #endif
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].hi[1] = %d\n", vv, bc[vv].hi[1]);
          //exit(1);
      }
    }
  }
  #endif


  #if (NDIM==3)
  // Inner k-boundary
  if (istart[2] == 0) {

    VLOOP {

      switch (bc[vv].lo[2]) {

        case OUTFLOW:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstart = KS(i,j,istart[2]);
              for (k=kstart-NG; k<kstart; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstart,vv);
              }
            }
          }
          break;

        case REFLECT:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstart = KS(i,j,istart[2]); 
              for(k=kstart-NG; k<kstart; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstart,vv);
                if (vv == U3) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= (-1.);
                }
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==2) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                }
                #endif
              }
            }
          }
          break;

        case PERIODIC:
          if (istop[2] == global_grid_dims[2]) {
            // Only do this if proc is periodic with itself; otherwise it's an MPI boundary
    GPU_PRAGMA(omp target teams distribute parallel for)
            ISLOOP(i) {
              JSLOOP(i,j) {
                kstart = KS(i,j,istart[2]);
                for (k=kstart-NG; k<kstart; k++) {
                  nkp = my_grid_dims[2]/DKS(i,j);
                  kp = k + nkp;
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kp,vv);
                }
              }
            }
          }
          break;

        case PROB:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstart = KS(i,j,istart[2]);
              for (k=kstart-NG; k<kstart; k++) {
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              }
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].lo[2] = %d\n", vv, bc[vv].lo[2]);
          //exit(1);
      }
    }
  }

  // Outer k-boundary
  if (istop[2] == global_grid_dims[2]) {
    
    VLOOP {

      switch (bc[vv].hi[2]) {

        case OUTFLOW:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstop = KS(i,j,istop[2]);
              for (k=kstop; k<kstop+NG; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstop-1,vv);
              }
            }
          }
          break;

        case REFLECT:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstop = KS(i,j,istop[2]);
              for(k=kstop; k<kstop+NG; k++) {
                NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kstop-1,vv);
                if (vv == U3) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= (-1.);
                }
                #if (DO_RADIATION)
                if (vv>=ifrad1 && (vv-ifrad1)%NDIM==2) {
                  NDP_ELEM_LINEAR(p,i,j,k,vv) *= -1.0;
                }
                #endif
              }
            }
          }
          break;

        case PERIODIC:
          if (istart[2] == 0) {
            // Only do this if proc is periodic with itself; otherwise it's an MPI boundary
    GPU_PRAGMA(omp target teams distribute parallel for)
            ISLOOP(i) {
              JSLOOP(i,j) {
                kstop = KS(i,j,istop[2]);
                for (k=kstop; k<kstop+NG; k++) {
                  nkp = my_grid_dims[2]/DKS(i,j);
                  kp = k - nkp;
                  NDP_ELEM_LINEAR(p,i,j,k,vv) = NDP_ELEM_LINEAR(p,i,j,kp,vv);
                }
              }
            }
          }
          break;

        case PROB:
    GPU_PRAGMA(omp target teams distribute parallel for)
          ISLOOP(i) {
            JSLOOP(i,j) {
              kstop = KS(i,j,istop[2]);
              for (k=kstop; k<kstop+NG; k++) {
                prob_bounds(i,j,k,&NDP_ELEM_LINEAR(p,i,j,k,0));
              }
            }
          }
          break;

        default:
          printf("INVALID BOUNDARY CONDITION FOR bc[%d].hi[2] = %d\n", vv, bc[vv].hi[2]);
          //exit(1);
      }
    }
  }
  #endif

  if (half_step_sync) t -= dt;

  TIMER_STOP;
  return;
}

