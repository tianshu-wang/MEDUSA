
#include "decs.h"

static int limit;

static double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2)
{
  return m[4*r0+c0] * (m[4*r1+c1] * m[4*r2+c2] - m[4*r2+c1] * m[4*r1+c2]) -
         m[4*r0+c1] * (m[4*r1+c0] * m[4*r2+c2] - m[4*r2+c0] * m[4*r1+c2]) +
         m[4*r0+c2] * (m[4*r1+c0] * m[4*r2+c1] - m[4*r2+c0] * m[4*r1+c1]);
}
 
 
static void adjoint(double m[16], double adjOut[16])
{
  adjOut[ 0] =  MINOR(m,1,2,3,1,2,3); adjOut[ 1] = -MINOR(m,0,2,3,1,2,3); adjOut[ 2] =  MINOR(m,0,1,3,1,2,3); adjOut[ 3] = -MINOR(m,0,1,2,1,2,3);
  adjOut[ 4] = -MINOR(m,1,2,3,0,2,3); adjOut[ 5] =  MINOR(m,0,2,3,0,2,3); adjOut[ 6] = -MINOR(m,0,1,3,0,2,3); adjOut[ 7] =  MINOR(m,0,1,2,0,2,3);
  adjOut[ 8] =  MINOR(m,1,2,3,0,1,3); adjOut[ 9] = -MINOR(m,0,2,3,0,1,3); adjOut[10] =  MINOR(m,0,1,3,0,1,3); adjOut[11] = -MINOR(m,0,1,2,0,1,3);
  adjOut[12] = -MINOR(m,1,2,3,0,1,2); adjOut[13] =  MINOR(m,0,2,3,0,1,2); adjOut[14] = -MINOR(m,0,1,3,0,1,2); adjOut[15] =  MINOR(m,0,1,2,0,1,2);
}
 
static double det(double m[16], double invOut[16])
{
  return m[0] * invOut[0] + //MINOR(m, 1, 2, 3, 1, 2, 3) -
         m[1] * invOut[4] + //MINOR(m, 1, 2, 3, 0, 2, 3) +
         m[2] * invOut[8] + //MINOR(m, 1, 2, 3, 0, 1, 3) -
         m[3] * invOut[12]; //MINOR(m, 1, 2, 3, 0, 1, 2);
}
 
 
static void invertRowMajor(double m[16], double invOut[16])
{
  int i;
  double inv_det;
  
  adjoint(m, invOut);

  inv_det = 1.0f / det(m,invOut);
  for (i=0; i<16; ++i) 
    invOut[i] = invOut[i] * inv_det;
}

void init_interp_order()
{
	int i;

	for (i=0; i<NPRIM+ncomp; i++) {
		interp_order[i] = hydro_interp_order;
	}

	#if (DO_RADIATION==TRUE)
	for (i=irad1; i<nvars; i++) {
		interp_order[i] = rad_interp_order;
	}
	#endif

	interp_order[nvars  ] = hydro_interp_order;	// pressure
	interp_order[nvars+1] = hydro_interp_order;	// Gamma

	return;
}


static double phihat(const double theta)
{
	return MAX(0, MIN((2.0+theta)/3.0, MAX(-0.5*theta, MIN(MIN(2.0*theta, (2.0 + theta)/3.0), 1.6))));
}

static double slope_limiter(const double d1, const double d2, const double r, const double dx)
{
	#if (1)
	const double theta = d1/(d2+1.0e-14);
	const double eta = (SQR(d1) + SQR(d2))/SQR(r*dx);

	if (eta <= 1-1.0e-14) {
		return (2.0 + theta)/3.0;
	} else if (eta >= 1.0 + 1.0e-14){
		return phihat(theta);
	} else {
		return 0.5*((1.0 - (eta-1.0)/1.0e-14) * (2.0+theta)/3.0 + (1.0 + (eta - 1.0)/1.0e-14) * phihat(theta));
	}
	#endif

	#if (0)
	return MAX(0,MIN(1,d2/d1));
	#endif
}

void lin_recon(const double* restrict a, double* restrict al, double* restrict ar, 
               const double* restrict const beta, const double* restrict const Gamma,
               const double* restrict const x, const int istart, const int istop)
{
	int i;
	double sp,sc,sm;
	double A,B,xm,xp,x0;
// double sg;

	//GPU_PRAGMA(omp target teams distribute parallel for)
	for (i=istart+NG; i<=istop+NG; i++) {
		xm = beta[i-1]/Gamma[i-1];
                xp = beta[i+1]/Gamma[i+1];
		x0 = beta[i  ]/Gamma[i  ];
                sm = (a[i  ] - a[i-1])/(x0-0.5*(x[i]+xm));     // maximum allowable derivative
                sp = (a[i+1] - a[i  ])/(0.5*(x[i+1]+xp)-x0);
                // sm = (a[i  ] - a[i-1])/(x0-xm);     // maximum allowable derivative
                // sp = (a[i+1] - a[i  ])/(xp-x0);
		if (sm*sp <= 0.0) {
			al[i] = ar[i] = a[i];  // local extremum, therefore flatten
			continue;
		}
                // else {
                // sg = 2.0*sm*sp/(sm+sp);
                // }

		sc = (a[i+1]-a[i-1])/(xp-xm);
                A = SIGN(sc)*MIN(MIN(fabs(sm),fabs(sp)),fabs(sc));
                // A = SIGN(sc)*MIN(2.0*MIN(fabs(sm),fabs(sp)),MIN(fabs(sc),fabs(sg)));
		B = a[i] - A*beta[i]/Gamma[i];
		al[i] = A*x[i  ] + B;
		ar[i] = A*x[i+1] + B;	
	}

	return;
}


void para_recon(const double* restrict a, double* restrict al, double* restrict ar,
                       const double* restrict const alpha, const double* restrict const beta, const double* restrict const Gamma, 
                       const double* restrict const x, const int istart, const int istop)
{
  // parabolic reconstruction
  // this isn't exactly TVD yet, but seems at least to be TVB

	int i;
	double dqm,dqp,dqc,den,dadx,xd;
	double A, B, C, Al, Bl, Cl, Ar, Br, Cr;
	double xzero,xmid,xlo,xhi,mask_sign,sign,mask_zero,mean;
	double a0b1g2,a0b2g1,a1b0g2,a1b2g0,a2b0g1,a2b1g0;
	static int count = 0;

	//GPU_PRAGMA(omp target teams distribute parallel for simd shared(alpha,beta,Gamma))
	for (i=istart+NG; i<=istop+NG; i++) {
		dqm = a[i  ] - a[i-1];
		dqp = a[i+1] - a[i  ];
		if (dqm*dqp <= 0.0) {
                al[i] = ar[i] = a[i];
			continue;
		}

		dqc = a[i+1] - a[i-1];
		a0b1g2 = alpha[i-1]*beta[i  ]*Gamma[i+1];
		a0b2g1 = alpha[i-1]*beta[i+1]*Gamma[i  ];
		a1b0g2 = alpha[i  ]*beta[i-1]*Gamma[i+1];
		a1b2g0 = alpha[i  ]*beta[i+1]*Gamma[i-1];
		a2b0g1 = alpha[i+1]*beta[i-1]*Gamma[i  ];
		a2b1g0 = alpha[i+1]*beta[i  ]*Gamma[i-1];
		den = 1.0/(a0b1g2 + a1b2g0 + a2b0g1 - a0b2g1 - a1b0g2 - a2b1g0);
		A =  (dqm* beta[i+1]*Gamma[i-1]*Gamma[i] + Gamma[i+1]*(dqp* beta[i-1]*Gamma[i] - dqc* beta[i]*Gamma[i-1]))*den;
		B = -(dqm*alpha[i+1]*Gamma[i-1]*Gamma[i] + Gamma[i+1]*(dqp*alpha[i-1]*Gamma[i] - dqc*alpha[i]*Gamma[i-1]))*den;
		C = a[i] - (alpha[i]*A + beta[i]*B)/Gamma[i];

		al[i] = x[i  ]*(A*x[i  ] + B) + C;
		ar[i] = x[i+1]*(A*x[i+1] + B) + C;

                mean = 0.75*a[i-1] + 0.25*a[i]; 
		mask_sign = (a[i] - al[i])*(al[i] - mean);
		al[i] = (mask_sign < 0.0) ? mean : al[i];

                mean = 0.75*a[i+1] + 0.25*a[i];
		mask_sign = (mean - ar[i])*(ar[i] - a[i]);
		ar[i] = (mask_sign < 0.0) ? mean : ar[i];

		// reconstruct the parabola (again)
		den = (x[i] - x[i+1])*(alpha[i] + Gamma[i]*x[i]*x[i+1] - beta[i]*(x[i]+x[i+1]));
		den = 1.0/den;
		A = ( beta[i]*(ar[i] - al[i]) + Gamma[i]*( a[i]*(x[i] - x[i+1]) - ar[i]*x[i] + al[i]*x[i+1]))*den;
		B = (alpha[i]*(al[i] - ar[i]) + Gamma[i]*(ar[i] - a[i])*SQR(x[i]) + Gamma[i]*(a[i]-al[i])*SQR(x[i+1]))*den;
		C = a[i] - (alpha[i]*A + beta[i]*B)/Gamma[i];

		xzero = -B/(2*A+1.0e-16);
		if (isinf(xzero) || isnan(xzero)) {
			xzero = -1.0e-200;
		}
		xlo = x[i];
		xhi = x[i+1];
		xmid = 0.5*(x[i]+x[i+1]);

		if (xzero > xlo && xzero < xhi) {
			//lin_recon(a, al, ar, flatten, i, i);
			if (xzero > xmid) {	// move zero to right
				A = Gamma[i]*(a[i] - ar[i])/(alpha[i] + xhi*(Gamma[i]*xhi - 2.0*beta[i]));
				B = -2.0*A*xhi;
				C = a[i] - (alpha[i]*A + beta[i]*B)/Gamma[i];
				al[i] = x[i]*(A*x[i] + B) + C;
			} else {	// move zero to left
				A = Gamma[i]*(a[i] - al[i])/(alpha[i] + xlo*(Gamma[i]*xlo - 2.0*beta[i]));
				B = -2.0*A*xlo;
				C = a[i] - (alpha[i]*A + beta[i]*B)/Gamma[i];
				ar[i] = x[i+1]*(A*x[i+1] + B) + C;
			}
		}
/*
		double flat_fact = flatten[i];//MAX(MAX(flatten[MAX(i-1,-1)],flatten[i]),flatten[MIN(i+1,istop)]);

		flat_fact = 1.;///(pow(flat_fact,8) + 1.);
		if (flat_fact < 0. || flat_fact > 1.) fprintf(stderr,"WHAT THE HELL HAPPENED\n");
		al[i] = flat_fact*al[i] + (1.-flat_fact)*a[i];
		ar[i] = flat_fact*ar[i] + (1.-flat_fact)*a[i];
*/
	}

	return;
}


//GPU_PRAGMA(omp declare target)
//void lin_recon(const double* restrict a, double* restrict al, double* restrict ar, 
//               const double* restrict const beta, const double* restrict const Gamma,
//               const double* restrict const x, const int istart, const int istop);
//void para_recon(const double* restrict a, double* restrict al, double* restrict ar,
//                       const double* restrict const alpha, const double* restrict const beta, const double* restrict const Gamma, 
//                       const double* restrict const x, const int istart, const int istop);
//GPU_PRAGMA(omp end declare target)

void interp(int ninterp,int *interp_order,double p[NSIZE][2*NG], double pl[NSIZE][2*NG], double pr[NSIZE][2*NG], const double* restrict alpha, const double* restrict beta, const double* restrict Gamma, const double* restrict x, const int size)
{
  // reconstruct along the provided pencil using the variable specific interpolation order
	int i,vv;
	double dp,pmin,max_flat,pavg;
	#if (DUMP_RECON==TRUE)
	void dump_recon(double **p, double **pl, double **pr, int start, int stop);
	#endif

    //TIMER_START("interp");

  // tag cells near large pressure jumps to flatten the reconstruction
	max_flat = 0.0;
	//GPU_PRAGMA(omp target teams distribute parallel for)
	//for (i=-1; i<=size; i++) {
	//	pavg = 3.0/(1.0/p[nvars][i-1+NG] + 1.0/p[nvars][i+NG] + 1.0/p[nvars][i+1+NG]);
	//	dp = p[nvars][i+NG];
	//	flatten[i+NG] = dp/pavg - 1.0;
	//}

	ILOOP {
		switch(interp_order[vv]) {
			//case 5:
			//	printf(stderr,"quartic reconstruction not implemented in this version\n");
			//	exit(1);
			//	//quar_recon(p[vv], pl[vv], pr[vv], alpha, beta, Gamma, x, -1, size);
			//	break;
			//case 4:
			//	fprintf(stderr,"cubic reconstruction not implemented in this version\n");
			//	exit(1);
			//	//cubic_recon(p[vv], pl[vv], pr[vv], alpha, beta, Gamma, x, -1, size);
			//	break;
			case 3:
				para_recon(p[vv], pl[vv], pr[vv], alpha, beta, Gamma, x, -1, size);
				break;
			case 2:
				lin_recon(p[vv], pl[vv], pr[vv], beta, Gamma, x, -1, size);
				break;
			default:
				for (i=-1; i<=size; i++) {
					pl[vv][i+NG] = pr[vv][i+NG] = p[vv][i+NG];
				}
				break;
		}
	}

	//#if (DUMP_RECON==TRUE)
	//if (interp_order[0] > 1 && size > 10) {
	//	dump_recon(p, pl, pr, -1, 10);
	//}
	//#endif

    //TIMER_STOP;
	return;
}


void dump_recon(double **p, double **pl, double **pr, int start, int stop)
{
	int i,vv;
	char name[128];
	FILE *fp;
	static int last_step = -1;
	static int tag;

	if (istep == last_step) {
		tag++;
	} else {
		tag = 0;
	}

	sprintf(name,"dumps/recon_p_%05d_%02d.dat", istep, tag);
	fp = fopen(name,"w");
	for (i=start; i<=stop; i++) {
		fprintf(fp,"%g", startx[0] + (i+0.5)*dx[0]);
		ILOOP {
			fprintf(fp," %g", p[vv][i+NG]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	sprintf(name,"dumps/recon_lr_%05d_%02d.dat", istep, tag);
	fp = fopen(name,"w");
	for (i=start; i<=stop; i++) {
		fprintf(fp,"%g", startx[0] + i*dx[0]);
		ILOOP {
			fprintf(fp," %g", pl[vv][i+NG]);
		}
		fprintf(fp,"\n");
		fprintf(fp,"%g", startx[0] + (i+1)*dx[0]);
		ILOOP {
			fprintf(fp," %g", pr[vv][i+NG]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	last_step = istep;

	return;
}
