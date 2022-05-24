
#include "../../decs.h"
#include <hdf5.h>

#define TABLE_TOL	1.e-10

#define EOS_ELEM(irho,iT,iY)	(Nrho*((iY)*NT + (iT)) + (irho))

static int Nrho,NT,NYe;
static double *tab_lrho;
static double *tab_lT;
static double *tab_Ye;
static double *tab_lP;
static double *tab_ent;
static double *tab_cs2;
static double *tab_le;
static double *tab_xn;
static double *tab_xp;

static double tab_lrho_min,tab_lrho_max;
static double tab_lT_min,tab_lT_max;
static double tab_Ye_min,tab_Ye_max;
static double tab_dlrho,tab_dlT,tab_dYe;

static double energy_shift;

void eos_init(char *name) {

	hsize_t file_grid_dims[3], file_start[3], file_count[3];
	hsize_t mem_grid_dims[3], mem_start[3];

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	hid_t file_id = H5Fopen(name, H5F_ACC_RDONLY, plist_id);
	H5Pclose(plist_id);

	hsize_t one = 1;
	file_grid_dims[0] = 1;
	hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = 1;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	hid_t memspace  = H5Screate_simple(1, &one, NULL);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	hid_t dset_id = H5Dopen(file_id, "pointsrho", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	int status = H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, plist_id, &Nrho);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "pointstemp", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, plist_id, &NT);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "pointsye", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_STD_I32LE, memspace, filespace, plist_id, &NYe);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "energy_shift", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, &energy_shift);
	H5Dclose(dset_id);
	H5Pclose(plist_id);


	H5Sclose(filespace);
	H5Sclose(memspace);

	int tab_size = NYe*NT*Nrho;

	tab_lrho = malloc_rank1(Nrho, sizeof(double));
	tab_lT = malloc_rank1(NT, sizeof(double));
	tab_Ye = malloc_rank1(NYe, sizeof(double));

	tab_lP = malloc_rank1(tab_size, sizeof(double));
	tab_ent = malloc_rank1(tab_size, sizeof(double));
	tab_cs2 = malloc_rank1(tab_size, sizeof(double));
	tab_le = malloc_rank1(tab_size, sizeof(double));
	tab_xn = malloc_rank1(tab_size, sizeof(double));
	tab_xp = malloc_rank1(tab_size, sizeof(double));

	hsize_t Narr = Nrho;
	file_grid_dims[0] = Narr;
	filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = Narr;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	memspace  = H5Screate_simple(1, &Narr, NULL);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "logrho", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_lrho);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
	H5Sclose(filespace);
	H5Sclose(memspace);

	Narr = NT;
	file_grid_dims[0] = Narr;
	filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = Narr;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	memspace  = H5Screate_simple(1, &Narr, NULL);
	
	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "logtemp", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_lT);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
	H5Sclose(filespace);
	H5Sclose(memspace);

	Narr = NYe;
	file_grid_dims[0] = Narr;
	filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = Narr;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	memspace  = H5Screate_simple(1, &Narr, NULL);
	
	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "ye", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_Ye);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
	H5Sclose(filespace);
	H5Sclose(memspace);

	file_grid_dims[0] = NYe;
	file_grid_dims[1] = NT;
	file_grid_dims[2] = Nrho;
	file_start[0] = 0;
	file_start[1] = 0;
	file_start[2] = 0;
	file_count[0] = NYe;
	file_count[1] = NT;
	file_count[2] = Nrho;

	filespace = H5Screate_simple(3, file_grid_dims, NULL);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

	mem_grid_dims[0] = NYe;
	mem_grid_dims[1] = NT;
	mem_grid_dims[2] = Nrho;
	memspace = H5Screate_simple(3, mem_grid_dims, NULL);
	mem_start[0] = 0;
	mem_start[1] = 0;
	mem_start[2] = 0;
	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "logpress", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_lP);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "logenergy", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_le);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "cs2", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_cs2);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "entropy", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_ent);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "Xn", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_xn);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	dset_id = H5Dopen(file_id, "Xp", plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	status = H5Dread(dset_id, H5T_IEEE_F64LE, memspace, filespace, plist_id, tab_xp);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	H5Sclose(filespace);
	H5Sclose(memspace);

	H5Fclose(file_id);

	tab_lrho_min = tab_lrho[0];
	tab_lrho_max = tab_lrho[Nrho-1];
	tab_dlrho = (tab_lrho_max - tab_lrho_min)/(Nrho-1);

	tab_lT_min = tab_lT[0];
	tab_lT_max = tab_lT[NT-1];
	tab_dlT = (tab_lT_max - tab_lT_min)/(NT-1);

	tab_Ye_min = tab_Ye[0];
	tab_Ye_max = tab_Ye[NYe-1];
	tab_dYe = (tab_Ye_max - tab_Ye_min)/(NYe-1);

}

double eos_interp(double lrho, double lT, double Ye, double* restrict tab) {

	int irho = (lrho - tab_lrho_min)/tab_dlrho;
	int iT = (lT - tab_lT_min)/tab_dlT;
	int iY = (Ye - tab_Ye_min)/tab_dYe;

	if(irho < 0) {
		irho = 0;
		lrho = tab_lrho_min;
	} else if(irho >= Nrho) {
		irho = Nrho-1;
		lrho = tab_lrho_max;
	}

	if(iT < 0) {
		iT = 0;
		lT = tab_lT_min;
	} else if (iT >= NT) {
		iT = NT-1;
		lT = tab_lT_max;
	}

	if(iY < 0) {
		iY = 0;
		Ye = tab_Ye_min;
	} else if(iY >= NYe) {
		iY = NYe-1;
		Ye = tab_Ye_max;
	}

	const double delrho = (lrho - tab_lrho[irho])/(tab_lrho[irho+1]-tab_lrho[irho]);
	const double delT = (lT - tab_lT[iT])/(tab_lT[iT+1]-tab_lT[iT]);
	const double delY = (Ye - tab_Ye[iY])/(tab_Ye[iY+1] - tab_Ye[iY]);

	return	  (1-delrho)	* (1-delT)	* (1-delY)	* tab[EOS_ELEM(irho,iT,iY)]
			+ delrho		* (1-delT)	* (1-delY)	* tab[EOS_ELEM(irho+1,iT,iY)]
			+ (1-delrho)	* delT		* (1-delY)	* tab[EOS_ELEM(irho,iT+1,iY)]
			+ delrho		* delT		* (1-delY)	* tab[EOS_ELEM(irho+1,iT+1,iY)]
			+ (1-delrho)	* (1-delT)	* delY		* tab[EOS_ELEM(irho,iT,iY+1)]
			+ delrho		* (1-delT)	* delY		* tab[EOS_ELEM(irho+1,iT,iY+1)]
			+ (1-delrho)	* delT		* delY		* tab[EOS_ELEM(irho,iT+1,iY+1)]
			+ delrho		* delT		* delY		* tab[EOS_ELEM(irho+1,iT+1,iY+1)];
}

int eos_bisect(const double lrho, double lT, const double Ye, double* restrict tab, double val, double *leosTemp) {

	int status = 1;
	double lTl,lTr,fl,fr;

	double grow = 0.01;
	do {
		double dTl = fabs(grow*lT) + 1.e-10;
		lTl = lT - dTl;
		fl = eos_interp(lrho, lTl, Ye, tab) - val;
		lTr = lT + dTl;
		fr = eos_interp(lrho, lTr, Ye, tab) - val;
		grow *= 1.1;
	} while(fl*fr > 0 && lTl > tab_lT_min && lTr < tab_lT_max);

	if(fl*fr > 0) {
		lTl = tab_lT_min;
		fl = eos_interp(lrho, lTl, Ye, tab) - val;
		lTr = tab_lT_max;
		fr = eos_interp(lrho, lTr, Ye, tab) - val;
		if(fl*fr > 0) return 0;
	}

	do {
		double lTm = 0.5*(lTl + lTr);
		double fm = eos_interp(lrho, lTm, Ye, tab) - val;
		if(fl*fm <= 0.) {
			lTr = lTm;
			fr = fm;
		} else {
			lTl = lTm;
			fl = fm;
		}
	} while(lTr-lTl > TABLE_TOL);

	*leosTemp = 0.5*(lTl + lTr);

	return status;
}

int eos_root(const double lrho, double lT, const double Ye, double* restrict tab, double val, double *leosTemp) {

	int status = 1;
	double lT_last;
	double lT_save = lT;

	// first try Newton iteration
	int iter = 0;
	do {
		double f = eos_interp(lrho, lT, Ye, tab) - val;
		double dlT = fabs(1.e-6*lT)+1.e-10;
		double df = (eos_interp(lrho, lT+dlT, Ye, tab) - eos_interp(lrho, lT-dlT, Ye, tab))/(2.*dlT);
		lT_last = lT;
		lT -= f/df;
		iter++;
	} while(iter < 10 && fabs(lT-lT_last)/(fabs(lT)+1.e-10) > TABLE_TOL);

	*leosTemp = lT;

	if((iter==10 && fabs(lT-lT_last) > TABLE_TOL) || (isnan(lT) || isinf(lT))) {
		// Newton failed, try bisection
		fprintf(stderr,"Newton failed, trying bisection: %g %g %g\n", lrho, val, Ye);
		double min_val = 1.e100;
		double max_val = -1.e100;
		for(iter=0;iter<=NT;iter++) {
			lT = tab_lT_min + iter*tab_dlT;
			double f = eos_interp(lrho, lT, Ye, tab);
			if(f < min_val) min_val = f;
			if(f > max_val) max_val = f;
		}
		fprintf(stderr,"min/max = %g %g\n", min_val, max_val);
		status = eos_bisect(lrho, lT_save, Ye, tab, val, leosTemp);
	}

	return status;
}
		
int eos_fill(double* restrict p, double* restrict eos) {

	double lTguess,leosTemp;

	double le = log10(p[UU]/p[RHO]);
	if(eos[TEMP] < 1.e-10) lTguess = 0.;
	else lTguess = log10(eos[TEMP]);

	double lrho = log10(p[RHO]);

	int success = eos_root(lrho, lTguess, p[YE], tab_le, le, &leosTemp);
	if(!success) {
		fprintf(stderr,"Failed to root find table %g %g %g %g   %g %d\n", lrho, lTguess, p[YE], leosTemp, p[RHO], mpi_io_proc());
		return 0;
	}

	eos[PRESS] = eos_interp(lrho, leosTemp, p[YE], tab_lP);
	eos[PRESS] = pow(10.,eos[PRESS]);

	eos[CS] = eos_interp(lrho, leosTemp, p[YE], tab_cs2);
	eos[GAMMA] = eos[CS]*p[RHO]/eos[PRESS];
	eos[GAMMA]  = MAX(eos[GAMMA],1.1);
	eos[CS] = sqrt(eos[CS]);

	eos[TEMP] = pow(10.,leosTemp);
	eos[ENT] = eos_interp(lrho, leosTemp, p[YE], tab_ent);

	return 1;
}

double u_given_rtx(double rho, double temp, double ye) {

	double e = eos_interp(log10(rho), log10(temp), ye, tab_le);
	e = pow(10.,e);
	return rho*e;
}

void eos_get_xnxp(double rho, double temp, double ye, double *xn, double *xp) {

	double lrho = log10(rho);
	double lT = log10(temp);

	*xn = eos_interp(lrho, lT, ye, tab_xn);
	*xp = eos_interp(lrho, lT, ye, tab_xp);
}
