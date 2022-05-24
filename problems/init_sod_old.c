
#include "../decs.h"

/* Sod's shock tube */

/* Last checked: 4/5/11 */

double init_model() {

	int i;
	double rhol,uul,ul;
	double rhor,uur,ur;
	double xl,xr;
	double dt;

	gam = 1.4;
	cour = 0.4;

	bc_choice[1] = OUTFLOW;
	bc_choice[2] = OUTFLOW;
	bc_choice[3] = OUTFLOW;

	rhol = 1.;
	uul = 1./(gam-1.);
	ul = 0.;

	rhor = 0.125;
	uur = 0.1/(gam-1.);
	ur = 0.;

	tmax = 0.25;
	dtdump = tmax/20.;
	dx[1] = 1.e100;
	dx[2] = 1.e100;
	dx[3] = 1./N3;

	startx[1] = startx[2] = startx[3] = 0.;

	for(i=0;i<N3;i++) {
		xl = i*dx[3];
		xr = xl + dx[3];

		if(xr < 0.5) {
			p[0][0][i][RHO] = rhol;
			p[0][0][i][UU] = uul;
			p[0][0][i][U3] = ul;
		} else if(xl < 0.5) {
			p[0][0][i][RHO] = rhol*(0.5-xl)/dx[3] + rhor*(xr-0.5)/dx[3];
			p[0][0][i][UU] = uul*(0.5-xl)/dx[3] + uur*(xr-0.5)/dx[3];
			p[0][0][i][U3] = ul*(0.5-xl)/dx[3] + ur*(xr-0.5)/dx[3];
		} else {
			p[0][0][i][RHO] = rhor;
			p[0][0][i][UU] = uur;
			p[0][0][i][U3] = ur;
		}

		p[0][0][i][U1] = p[0][0][i][U2] = 0.;

	}

	bc(p);

	dt = 0.5*stepsize(CONST_GRIDFUNC p);

	return MIN(dt, dtdump);
}
		
	


