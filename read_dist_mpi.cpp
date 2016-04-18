
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string.h> // for memset
#include "timer.h"
#include "./special_funcs/specialfunctions.h" // for erf function
#include <mpi.h>

using namespace std;

// flags
#define READ_PLASMA_FILE          1
#define READ_NEUTRAL_FILE_MOMENTS 1
#define COMPUTE_SOURCE_TERMS      0  // requires moments as well

#define PRINT_PLASMA_VTK          1
#define PRINT_NEUTRAL_MOMENTS_VTK 1
#define PRINT_SOURCE_TERMS        0

// MPI constants
#define ROOT 0
MPI_Datatype mpi_vel_dist_header;

// all other constants
const double EPS     = 1.0E-6;
const double AU_to_M = 149597870700.0;
const double J_to_EV = 6.242E18;
const double MP      = 1.6726219E-27; // mass of proton in kg
const double KB      = 1.3806488E-23; // Boltzmann constant in J/K
const double third_root_16_over_pi = pow(16.0/M_PI,1.0/3.0);

struct vel_dist_header {
	int    nxr;    int    nvr;
	int    nxt;    int    nvt;
	int    nxp;    int    nvp;
	double xalpha; double valpha;
	double xmin;   double vmin;
	double xmax;   double vmax;
};

struct grid_cell_sph {
	double *gxr, *cxr, *gxt, *cxt, *gxp, *cxp;
	double *gvr, *cvr, *gvt, *cvt, *gvp, *cvp;
};

// i = rho, j = theta, k = phi
// nn = [nr, nt, np]
int arrayIdx3(int i, int j, int k, int* nn) 
{
	return k + nn[2]*(j + nn[1]*i);
}

// Analytic expression for charge exchange rate from maxwellian protons
//  Example vars:
// dens = local density of protons
// temp  = temperature of protons
// sig_ex = pointer to function of the charge exchange (requires units of 1keV)
// bulkv = local bulk flow of protons as [x, y, z] double array
// m     = mass of protons
// v     = velocity of comparing distribution (neutrals) as [x, y, z] double array
double beta_erf(double dens, double temp, double (*sig_ex)(double),
	double* bulkv, double m, double* v)
{
	double vth = sqrt(2.0 * KB * temp / m); // m/s
	double sqrel_v = (v[0]-bulkv[0]) * (v[0]-bulkv[0])
		           + (v[1]-bulkv[1]) * (v[1]-bulkv[1])
		           + (v[2]-bulkv[2]) * (v[2]-bulkv[2]);
	double rel_v = sqrt(sqrel_v) / vth;
	double rel_energy = J_to_EV * 0.5 * m * sqrel_v / 1000.0; // keV
	
	return dens * vth * sig_ex(rel_energy)
		* ((rel_v + 1.0/(2.0 * rel_v))* alglib::errorfunction(rel_v)  
		+ exp(-rel_v * rel_v) / sqrt(M_PI));
}

// A more simplified version of the above function, which has a maximum error
// of 3% according to Pauls et. al. 1995
double beta_sqrt(double dens, double temp, double (*sig_ex)(double),
	double* bulkv, double m, double* v)
{
	double vth = sqrt(2.0 * KB * temp / m); // m/s
	double sqrel_v = (v[0]-bulkv[0]) * (v[0]-bulkv[0])
		           + (v[1]-bulkv[1]) * (v[1]-bulkv[1])
		           + (v[2]-bulkv[2]) * (v[2]-bulkv[2]);
	//double rel_v = sqrt(sqrel_v) / vth;
	double rel_energy = J_to_EV * 0.5 * m * sqrel_v / 1000.0; // keV
	
	return dens * vth * sig_ex(rel_energy)
		* sqrt(4.0/M_PI + sqrel_v /(vth*vth));
}

// calculate charge exchange rate, and total source term of momentum and energy
//  using Maxwellian protons; exact solution
void source_maxwellian_erf(double dens, double temp, double (*sig_ex)(double),
	double* bulkv, double m, double *v,
	double &beta_ij, double *source_m, double &source_E)
{
	double vth = sqrt(2.0 * KB * temp / m); // m/s
	double sqrel_v = (v[0]-bulkv[0]) * (v[0]-bulkv[0])
		           + (v[1]-bulkv[1]) * (v[1]-bulkv[1])
		           + (v[2]-bulkv[2]) * (v[2]-bulkv[2]);
	double rel_v = sqrt(sqrel_v) / vth;
	double rel_energy = J_to_EV * 0.5 * m * sqrel_v / 1000.0; // keV
	double vj_dot_ui  = (v[0]*(v[0] - bulkv[0]) + v[1]*(v[1] - bulkv[1]) + v[2]*(v[2] - bulkv[2]) ) / (vth * vth);
	
	double erf_x  = alglib::errorfunction(rel_v);
	double exp_x2_pi = exp(-rel_v * rel_v) / sqrt(M_PI);
	double ch_ex  = sig_ex(rel_energy);
	
	// update output variables
	// charge exchange rate
	beta_ij = dens * vth * ch_ex * ( (rel_v + 1.0/(2.0 * rel_v)) * erf_x
	        + exp_x2_pi );
	
	// total momentum source term
	double sqrt_2nd_moment = (rel_v + 1.0/rel_v - 1.0 / (4.0 * rel_v*rel_v*rel_v)) * erf_x
		+ (1.0 + 1.0 / (2.0 * rel_v*rel_v)) * exp_x2_pi;
	double source_m_coeff = m * dens * vth * ch_ex * sqrt_2nd_moment;
	source_m[0] = (v[0]-bulkv[0]) * source_m_coeff;
	source_m[1] = (v[1]-bulkv[1]) * source_m_coeff;
	source_m[2] = (v[2]-bulkv[2]) * source_m_coeff;
	
	// total energy source term
	source_E = 0.5 * m * vth*vth*vth * dens * ch_ex
		* ( -(rel_v * (rel_v*rel_v + 3.0) + 3.0 / (4.0 * rel_v) ) * erf_x
		- (rel_v*rel_v + 5.0/2.0 ) * exp_x2_pi + 2.0 * vj_dot_ui * sqrt_2nd_moment);
}

// using square root approximations to avoid exp and erf calls
void source_maxwellian_sqrt(double dens, double temp, double (*sig_ex)(double),
	double* bulkv, double m, double *v,
	double &beta_ij, double *source_m, double &source_E)
{
	double vth = sqrt(2.0 * KB * temp / m); // m/s
	double sqrel_v = (v[0]-bulkv[0]) * (v[0]-bulkv[0])
		           + (v[1]-bulkv[1]) * (v[1]-bulkv[1])
		           + (v[2]-bulkv[2]) * (v[2]-bulkv[2]);
	double rel_v = sqrt(sqrel_v) / vth;
	double rel_energy = J_to_EV * 0.5 * m * sqrel_v / 1000.0; // keV
	double vj_dot_ui  = (v[0]*(v[0] - bulkv[0]) + v[1]*(v[1] - bulkv[1]) + v[2]*(v[2] - bulkv[2]) ) / (vth * vth);
	
	//double erf_x  = alglib::errorfunction(rel_v);
	//double exp_x2_pi = exp(-rel_v * rel_v) / sqrt(M_PI);
	double ch_ex  = sig_ex(rel_energy);
	
	// update output variables
	// charge exchange rate
	beta_ij = dens * vth * ch_ex * sqrt(4.0/M_PI + rel_v*rel_v);
	
	// total momentum source term
	double sqrt_2nd_moment = sqrt(64.0/(9.0*M_PI) + rel_v*rel_v);
	double source_m_coeff = m * dens * vth * ch_ex * sqrt_2nd_moment;
	source_m[0] = (v[0]-bulkv[0]) * source_m_coeff;
	source_m[1] = (v[1]-bulkv[1]) * source_m_coeff;
	source_m[2] = (v[2]-bulkv[2]) * source_m_coeff;
	
	// total energy source term
	source_E = 0.5 * m * vth*vth*vth * dens * ch_ex
		* ( -sqrt((third_root_16_over_pi +  rel_v*rel_v)*(third_root_16_over_pi
		+rel_v*rel_v)*(third_root_16_over_pi +  rel_v*rel_v))//pow(third_root_16_over_pi +  rel_v*rel_v,1.5) 
		+ 2.0 * vj_dot_ui * sqrt_2nd_moment);
}

// Assuming a Kappa distribution, in the same limit as beta_sqrt, ie x >> 1,
// from Heerikhuisen et. al. 2008. Note, kappa < 171.6 for the gamma function
double beta_sqrt_kappa(double dens, double temp, double (*sig_ex)(double),
	double* bulkv, double m, double* v, double kappa)
{
	double vth = sqrt(2.0 * KB * temp / m) / sqrt(kappa/(kappa-1.5)); // m/s
	double sqrel_v = (v[0]-bulkv[0]) * (v[0]-bulkv[0])
		           + (v[1]-bulkv[1]) * (v[1]-bulkv[1])
		           + (v[2]-bulkv[2]) * (v[2]-bulkv[2]);
	double g_k1  = alglib::gammafunction(kappa + 1.0);
	double g_k12 = alglib::gammafunction(kappa - 0.5);
	double rel_energy = J_to_EV * 0.5 * m * sqrel_v / 1000.0; // keV
	
	return dens * vth * sig_ex(rel_energy)
		* sqrt( 4.0 *g_k1*g_k1 / (M_PI * kappa * (kappa+1.0) * g_k12*g_k12)
		+ sqrel_v /(vth*vth) );
}

// This function will use the kappa approx if the region is the inner heliosheath
// for the protons (for example), found heristically by the temperature
double beta_sqrt_kappa_IHS(double dens, double temp, double (*sig_ex)(double),
	double* bulkv, double m, double* v, double kappa, double LISM_temp)
{
	if(temp > 10.0 * LISM_temp) // log(T/T_LISM) > 1.6, 39.811
		return beta_sqrt_kappa(dens, temp, sig_ex, bulkv, m, v, kappa);
	else
		return beta_sqrt(dens, temp, sig_ex, bulkv, m, v);
}

// From Lindsay and Stebbings 2005 for H * H+ interactions
double charge_ex(double rel_eng) { // rel_eng in units of keV
	const double a1 = 4.15, a2 = 0.531, a3 = 67.3;
	// only calculate the hard tail if high enough energy
	double tail = (rel_eng < 35.0  ? 1.0 : pow(1.0 - exp(-a3/rel_eng), 4.5));
	// there needs to be enough energy to exchange the electron
	//double core = (rel_eng < 0.000005 ? 0.0 : (a1 - a2*log(rel_eng))); // (a1 - a2*log(0.005))(a1 - a2*log(rel_eng))
	double core = (a1 - a2*log(rel_eng));
	return core * core * tail;
}

// insertion sort is slow as hell, expecially with millions of vals to sort
void insert_sort_key(int* val, int* key, int size)
{
	//cout << "size = " << size << endl;
	int x, k, i, j;
	for(i=1; i < size; i++) {
		x = val[i];
		k = key[i];
		j = i;
		while(j > 0 && (key[j-1] > k || (key[j-1] == k && val[j-1] > x)))
		{
			val[j] = val[j-1];
			key[j] = key[j-1];
			j--;
		}
		val[j] = x;
		key[j] = k;
		//if(i % (size/100) == 0) cout << i << endl;
	}
}

// To use: mergesort(array, temparray, keyarray, tempkeyarray, minindex, maxindex)
// if the key is not unique, and the array is already sorted, the merge sort
// is stable, ie it respects the order of array with the same keys
void merge(int* A, int* Ak, int* B, int* Bk, int p, int q, int r)
{
	int i = p, k = p;
	int j = q + 1;
	
	while(i <= q && j <= r) { // while both subarrays are non empty
		if(Ak[i] <= Ak[j]) {
			Bk[k]  = Ak[i];
			B[k++] = A[i++]; // copy from left sub array
		} else {
			Bk[k]  = Ak[j];
			B[k++] = A[j++]; // copy from right sub array
		}
	}
	while(i <= q) {
		Bk[k]  = Ak[i];
		B[k++] = A[i++]; // copy any left over to B
	}
	while(j <= r) {
		Bk[k]  = Ak[j];
		B[k++] = A[j++];
	}
	memcpy(A+p,  B+p,  (1 + r - p)*sizeof(int)); // copy B back to A
	memcpy(Ak+p, Bk+p, (1 + r - p)*sizeof(int));
}

void mergesort(int* A, int* Ak, int* B, int* Bk, int p, int r) {
	if(p < r) {
		int q = (p + r) / 2;
		mergesort(A, Ak, B, Bk, p, q);
		mergesort(A, Ak, B, Bk, q + 1, r);
		merge(A, Ak, B, Bk, p, q, r);
	}
}

void display_header(vel_dist_header& h) {
	cout <<"nxr = "<< h.nxr <<"  nxt = "<< h.nxt <<"  nxp = "<< h.nxp << endl;
	cout <<"nvr = "<< h.nvr <<"  nvt = "<< h.nvt <<"  nvp = "<< h.nvp << endl;
	cout <<"xalpha = " << h.xalpha <<"  xmin = " << h.xmin <<"  xmax  = " << h.xmax << endl;
	cout <<"valpha = " << h.valpha <<"  vmin = " << h.vmin <<"  vmax  = " << h.vmax << endl;
}

void init_grid_cell_sph(grid_cell_sph& gc, vel_dist_header& h) {
	int i;
	// make a contiguous chunk of mem and manually position the pointers
	int size = 2*(h.nxr + h.nxt + h.nxp + h.nvr + h.nvt + h.nvp) + 6;
	gc.gxr = new double[size];
	gc.cxr = gc.gxr + h.nxr + 1;
	gc.gxt = gc.cxr + h.nxr;
	gc.cxt = gc.gxt + h.nxt + 1;
	gc.gxp = gc.cxt + h.nxt;
	gc.cxp = gc.gxp + h.nxp + 1;
	
	gc.gvr = gc.cxp + h.nxp;
	gc.cvr = gc.gvr + h.nvr + 1;
	gc.gvt = gc.cvr + h.nvr;
	gc.cvt = gc.gvt + h.nvt + 1;
	gc.gvp = gc.cvt + h.nvt;
	gc.cvp = gc.gvp + h.nvp + 1;
	
	// radial position grid and cell locations
	for(i=0; i <= h.nxr; i++)
		gc.gxr[i] = h.xmin + (h.xmax - h.xmin)
			* (exp(h.xalpha*i/double(h.nxr)) - 1.0)
			/ (exp(h.xalpha) - 1.0);
	for(i=0; i < h.nxr; i++)
		gc.cxr[i] = 0.5 * (gc.gxr[i] + gc.gxr[i+1]);
	// theta position grid and cell locations
	for(i=0; i <= h.nxt; i++)
		gc.gxt[i] = M_PI * i / double(h.nxt);
	for(i=0; i < h.nxt; i++)
		gc.cxt[i] = 0.5 * (gc.gxt[i] + gc.gxt[i+1]);
	// phi position grid and cell locations
	for(i=0; i < h.nxp; i++)
		gc.cxp[i] = 2.0 * M_PI * i / double(h.nxp);
	gc.gxp[0] = -M_PI / double(h.nxp);
	for(i=0; i < h.nxp; i++)
		gc.gxp[i+1] = gc.cxp[i] + M_PI / double(h.nxp);
	
	// radial velocity grid and cell locations
	for(i=0; i <= h.nvr; i++)
		gc.gvr[i] = h.vmin + (h.vmax - h.vmin)
			* (exp(h.valpha*i/double(h.nvr)) - 1.0)
			/ (exp(h.valpha) - 1.0);
	for(i=0; i < h.nvr; i++)
		gc.cvr[i] = 0.5 * (gc.gvr[i] + gc.gvr[i+1]);
	// theta velocity grid and cell locations
	for(i=0; i <= h.nvt; i++)
		gc.gvt[i] = M_PI * i / double(h.nvt);
	for(i=0; i < h.nvt; i++)
		gc.cvt[i] = 0.5 * (gc.gvt[i] + gc.gvt[i+1]);
	// phi velocity grid and cell locations
	for(i=0; i < h.nvp; i++)
		gc.cvp[i] = 2.0 * M_PI * i / double(h.nvp);
	gc.gvp[0] = -M_PI / double(h.nvp);
	for(i=0; i < h.nvp; i++)
		gc.gvp[i+1] = gc.cvp[i] + M_PI / double(h.nvp);
	
}

void delete_grid_cell_sph(grid_cell_sph& gc) {
	delete[] gc.gxr;
}

void display_grid_cell_sph(grid_cell_sph& gc, vel_dist_header& h) {
	int i;
	// display grid
	cout << "\n  gxr[i]: \n";
	for(i=0; i <= h.nxr; i++)
		cout << gc.gxr[i] << endl;
	cout << "\n  gxt[i]: \n";
	for(i=0; i <= h.nxt; i++)
	 	cout << gc.gxt[i] << endl;
	cout << "\n  gxp[i]: \n";
	for(i=0; i <= h.nxp; i++)
	 	cout << gc.gxp[i] << endl;
	
	cout << "\n  gvr[i]: \n";
	for(i=0; i <= h.nvr; i++)
		cout << gc.gvr[i] << endl;
	cout << "\n  gvt[i]: \n";
	for(i=0; i <= h.nvt; i++)
	 	cout << gc.gvt[i] << endl;
	cout << "\n  gvp[i]: \n";
	for(i=0; i <= h.nvp; i++)
	 	cout << gc.gvp[i] << endl;	
	 	
	// display cells
	cout << "\n  cxr[i]: \n";
	for(i=0; i < h.nxr; i++)
		cout << gc.cxr[i] << endl;
	cout << "\n  cxt[i]: \n";
	for(i=0; i < h.nxt; i++)
	 	cout << gc.cxt[i] << endl;
	cout << "\n  cxp[i]: \n";
	for(i=0; i < h.nxp; i++)
	 	cout << gc.cxp[i] << endl;
	
	cout << "\n  cvr[i]: \n";
	for(i=0; i < h.nvr; i++)
		cout << gc.cvr[i] << endl;
	cout << "\n  cvt[i]: \n";
	for(i=0; i < h.nvt; i++)
	 	cout << gc.cvt[i] << endl;
	cout << "\n  cvp[i]: \n";
	for(i=0; i < h.nvp; i++)
	 	cout << gc.cvp[i] << endl;
}

void ff_bskipRead(ifstream& binfile, int type, int n) // n = size of array of type
{
	char* t_char;
	int* t_int;
	float* t_float;
	double* t_double;
	
	switch(type) {
		case 0: // char
			t_char = new char[n];
			binfile.read (reinterpret_cast<char*>(t_char), n);
			delete[] t_char;
			break;
		case 1: // int
			t_int = new int[n];
			binfile.read (reinterpret_cast<char*>(t_int), 4*n);
			delete[] t_int;
			break;
		case 2: // float
			t_float = new float[n];
			binfile.read (reinterpret_cast<char*>(t_float), 4*n);
			delete[] t_float;
			break;
		case 3: // double
			t_double = new double[n];
			binfile.read (reinterpret_cast<char*>(t_double), 8*n);
			delete[] t_double;
			break;
		default:
			cout << "ERROR: Undefined datatype in reading binfile.\n";
	}
}

template<class T>
bool ff_bdataRead(ifstream& binfile, int type, int n, int loopidxnum,
	T* data_array, int i, int j, int k, int* nn)
{
	int countstart, idx;
	bool good_read = 1;
	if(loopidxnum == 0)      countstart = i;
	else if(loopidxnum == 1) countstart = j;
	else if(loopidxnum == 2) countstart = k;
	else {
		cout << "\nERROR: loopidxnum = {0, 1, 2}\n\n";
		return 0; // error in read
	}
	
	T t_read[1];
	
	for(int count = countstart; count < n+countstart; count++) {
		if(loopidxnum == 0)
			idx = arrayIdx3(count, j, k, nn);
		else if(loopidxnum == 1)
			idx = arrayIdx3(i, count, k, nn);
		else if(loopidxnum == 2)
			idx = arrayIdx3(i, j, count, nn);
			
		binfile.read (reinterpret_cast<char*>(t_read), sizeof(T));
		data_array[idx] = t_read[0];
	}
	
	return good_read;
}

void print2vtk(const char* vtkfname, vel_dist_header& h, grid_cell_sph& gc,
	double *dens, double *velx, double *vely, double *velz,
	double *vxvar, double *vyvar, double *vzvar, double *enrg, int PARAM)
{
	if(PARAM == 0) {
		cout << " No valid parameters, exit print2vtk\n\n";
		return;
	}
	int ixr, ixt, ixp, dk, NUM_CELL_VALS, i, j, k, idx;
	int cell_type, cell_iter = 0;
	float x, y, z, xy;
	int POINTS = (h.nxr+1)*(h.nxt+1)*(h.nxp+1);
	int CELLS  = h.nxr * h.nxt * h.nxp;
	int *cell_type_array;
	cell_type_array = new int[CELLS];
	int nxp[3] = {h.nxr+1, h.nxt+1, h.nxp+1}; // these 1's are very important!!!
	int nxc[3] = {h.nxr, h.nxt, h.nxp};
	// set up .vtk header
	ofstream fvtk(vtkfname);
	fvtk << "# vtk DataFile Version 3.1\n";
	fvtk << "convert 6D neutral data from jh helio to .vtk format\n";
	fvtk << "ASCII\n";
	// set up POINTS header section
	fvtk << "DATASET UNSTRUCTURED_GRID\n";
	fvtk << "POINTS " << POINTS << " FLOAT\n";
	
	// Note: Many of the degenerate points are not used (origin and/or z-axis)
	//	but for indexing purposes they are still here
	for(ixr = 0; ixr <= h.nxr; ixr++)  // rho
		for(ixt = 0; ixt <= h.nxt; ixt++)  // theta
			for(ixp = 0; ixp <= h.nxp; ixp++) { // phi
				//sph2cart(rgrid[k], tgrid[j], pgrid[i], x, y, z);
				z  = gc.gxr[ixr] * cos(gc.gxt[ixt]);
		  		xy = gc.gxr[ixr] * sin(gc.gxt[ixt]);
		  		x  = xy * cos(gc.gxp[ixp]);
		  		y  = xy * sin(gc.gxp[ixp]);
				fvtk << x << ' ' << y << ' ' << z << endl;
			}
	cout << " POINTS finished\n";
	
	// get number of cell values
	NUM_CELL_VALS = 0;
	for(ixr = 0; ixr < h.nxr; ixr++) { // rho
		for(ixt = 0; ixt < h.nxt; ixt++) { // theta
			for(ixp = 0; ixp < h.nxp; ixp++) { // phi
				// At origin and + z-axis
				if( gc.gxr[ixr] < EPS && abs(gc.gxt[ixt]) < EPS ) {
					NUM_CELL_VALS += 5; // tetra
				// At origin and - z-axis
				} else if( gc.gxr[ixr] < EPS && abs(gc.gxt[ixt+1] - M_PI) < EPS ) {
					NUM_CELL_VALS += 5; // tetra
				// At origin, but not z-axis
				} else if( gc.gxr[ixr] < EPS ) {
					NUM_CELL_VALS += 6; // pyramid
				// At + z-axis, but not origin
				} else if( abs(gc.gxt[ixt]) < EPS ) {
					NUM_CELL_VALS += 7; // wedge
				// At - z-axis, but not origin
				} else if( abs(gc.gxt[ixt+1] - M_PI) < EPS ) {
					NUM_CELL_VALS += 7; // wedge
				// everything else
				} else {
					NUM_CELL_VALS += 9; // hexahedron
				}
			}
		}
	}
	cout << " GET NUMBER OF CELL VALUES finished\n";
	// assuming rho = [rho_min, rho_max], and full domain of theta and phi
	// set up CELLS header section and create CELL data
	fvtk << "\nCELLS " << CELLS << ' ' << NUM_CELL_VALS << endl;
	for(i = 0; i < h.nxr; i++) { // rho
		for(j = 0; j < h.nxt; j++) { // theta
			for(k = 0; k < h.nxp; k++) { // phi
				// we want phi to be connected at the last cell
				dk = (k == h.nxp-1 ? 0 : k+1);
				
				// At origin and + z-axis
				if( gc.gxr[i] < EPS && abs(gc.gxt[j]) < EPS ) {
					cell_type = 10; // tetra
					fvtk << 4;
					fvtk << ' ' << arrayIdx3(i, 0, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);
				// At origin and - z-axis
				} else if( gc.gxr[i] < EPS && abs(gc.gxt[j+1] - M_PI) < EPS ) {
					cell_type = 10; // tetra
					fvtk << 4;
					fvtk << ' ' << arrayIdx3(i, 0, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, dk, nxp);
				// At origin, but not z-axis
				} else if( gc.gxr[i] < EPS ) {
					cell_type = 14; // pyramid
					fvtk << 5;
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, dk, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);
					fvtk << ' ' << arrayIdx3(i, 0, 0, nxp);
				// At + z-axis, but not origin
				} else if( abs(gc.gxt[j]) < EPS ) {
					cell_type = 13; // wedge
					fvtk << 6;
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, 0, nxp);
					fvtk << ' ' << arrayIdx3(i, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i, j+1, dk, nxp);
					fvtk << ' ' << arrayIdx3(i, j, 0, nxp);
				// At - z-axis, but not origin
				} else if( abs(gc.gxt[j+1] - M_PI) < EPS ) {
					cell_type = 13; // wedge
					fvtk << 6;
					fvtk << ' ' << arrayIdx3(i+1, j+1, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, dk, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i, j+1, 0, nxp);
					fvtk << ' ' << arrayIdx3(i, j, dk, nxp);
					fvtk << ' ' << arrayIdx3(i, j, k, nxp);
				// everything else
				} else {
					cell_type = 12; // hexahedron
					fvtk << 8;
					fvtk << ' ' << arrayIdx3(i, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i, j, dk, nxp);          
					fvtk << ' ' << arrayIdx3(i+1, j, dk, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);          
					fvtk << ' ' << arrayIdx3(i, j+1, dk, nxp);
				}

				fvtk << endl;
				cell_type_array[cell_iter] = cell_type;
				cell_iter++;
			}
		}
	}
	
	fvtk << "\nCELL_TYPES " << CELLS << endl;
	for(k=0; k < CELLS; k++)
		fvtk << cell_type_array[k] << endl; 
	cout << " CELLS finished\n";
	delete[] cell_type_array;
	// set up the data
	fvtk << "\nCELL_DATA " << CELLS << endl;
	if(PARAM & 1) {
		// output density
		fvtk << "SCALARS density FLOAT\n";
		fvtk << "LOOKUP_TABLE default\n";
		for(ixr = 0; ixr < h.nxr; ixr++) { // rho
			for(ixt = 0; ixt < h.nxt; ixt++) { // theta
				for(ixp = 0; ixp < h.nxp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << dens[idx] << endl;
				}
			}
		}
		cout << " density data finished\n";
	}
	if((PARAM & 2) && (PARAM & 4) && (PARAM & 8)) {
		// output bulk_flow data
		fvtk << "\nVECTORS bulk_flow FLOAT\n";
		for(ixr = 0; ixr < h.nxr; ixr++) { // rho
			for(ixt = 0; ixt < h.nxt; ixt++) { // theta
				for(ixp = 0; ixp < h.nxp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << velx[idx] << ' '
						 << vely[idx] << ' '
						 << velz[idx] << endl;
				}
			}
		}
		cout << " bulk_flow vector data finished\n";
	}
	if(PARAM & 16) {
		// output temperature data
		fvtk << "SCALARS temperature FLOAT\n";
		fvtk << "LOOKUP_TABLE default\n";
		for(ixr = 0; ixr < h.nxr; ixr++) { // rho
			for(ixt = 0; ixt < h.nxt; ixt++) { // theta
				for(ixp = 0; ixp < h.nxp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					// 40.4085783 = mass of hydrogen / boltzman constant / 3
					// = [K / (km/s)^2]
					if((PARAM & 16) && (PARAM & 32) && (PARAM & 64)) // input are the variances
						fvtk << 40.4085783 * (vxvar[idx] + vyvar[idx] + vzvar[idx])
						     << endl;
					else if(PARAM & 16) // the temperature was precalculated
						fvtk << vxvar[idx] << endl;
				}
			}
		}
		cout << " temperature data finished\n";
	}
	if(PARAM & 128) {
		// output energy data
		fvtk << "SCALARS energy FLOAT\n";
		fvtk << "LOOKUP_TABLE default\n";
		for(ixr = 0; ixr < h.nxr; ixr++) { // rho
			for(ixt = 0; ixt < h.nxt; ixt++) { // theta
				for(ixp = 0; ixp < h.nxp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << enrg[idx] << endl;
				}
			}
		}
		cout << " energy data finished\n";
	}
	cout << vtkfname << " printing success!\n\n";
	fvtk.close();
}

int read_plasma_binary(const char* pfname, vel_dist_header& h, double* pdens,
	double* pvelx, double* pvely, double* pvelz, double* ptemp)
{
	int i, j, k, idx;
	// hard coded numbers from vel_dist_moments.m
	double gamma     = 5.0/3.0;
	double dens_LISM = 0.082; // cm^-3
	double v_LISM    = 26.0;  // km/s
	double temp_LISM = 8000.0;// K
	double AMB2      = v_LISM*v_LISM * 72.69 / (2.0*temp_LISM);
	ifstream pfile(pfname, ios::in|ios::binary|ios::ate);
	// check if file is opened correctly, exit if not
	cout << endl << pfname << endl;
	if(pfile.good()) {
		cout << " File open successfully\n";
	} else {
		cout << " File open failed\n";
		return -1;
	}
	streampos file_size = pfile.tellg(); // get size of file in bytes
	cout << "Size of file = " << file_size/1024.0/1024.0 << " MB\n\n";
	pfile.seekg(0, ios::beg); // start at beginning of file
	
	int nx[3] = {h.nxr, h.nxt, h.nxp}; 
	
	// skip ghost cells, following "vel_dist_moments.m"
	for(i=0; i < 2; i++) {
	  for(j=0; j < h.nxt; j++) {
		for(k=0; k < 8; k++) {
		  ff_bskipRead(pfile, 1, 1);
		  ff_bskipRead(pfile, 3, 2 + h.nxr);
		  ff_bskipRead(pfile, 1, 1);
		  ff_bskipRead(pfile, 3, 2);
		}
		ff_bskipRead(pfile, 1, 6 + h.nxr);
	  }
	}
	// note this is opposite of the way we access the neutral file,
	// which is phi in the inner most loop, not rho
	// so, we will force it to match the same access pattern as the neutral file
	for(k=0; k < h.nxp; k++) { // phi
	    for(j=0; j < h.nxt; j++) { // theta
	    	ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bdataRead(pfile, 3, h.nxr, 0, pdens, 0, j, k, nx);
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bdataRead(pfile, 3, h.nxr, 0, ptemp, 0, j, k, nx);
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bdataRead(pfile, 3, h.nxr, 0, pvelx, 0, j, k, nx);
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bdataRead(pfile, 3, h.nxr, 0, pvely, 0, j, k, nx);
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bdataRead(pfile, 3, h.nxr, 0, pvelz, 0, j, k, nx);
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 3, h.nxr); // skip Bx
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 3, h.nxr); // skip By
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 3, h.nxr); // skip Bz
		    
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 3, 2);
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 1, 2);
		    ff_bskipRead(pfile, 1, h.nxr); // skip region
		    ff_bskipRead(pfile, 1, 1); ff_bskipRead(pfile, 1, 2);
		    
		    // convert newly read data to physical units
		    for(i=0; i < h.nxr; i++) {
		      idx = arrayIdx3(i, j, k, nx);
		      ptemp[idx] *= AMB2 * gamma / pdens[idx] * temp_LISM;
		      pdens[idx] *= dens_LISM;
		      pvelx[idx] *= v_LISM;
		      pvely[idx] *= v_LISM;
		      pvelz[idx] *= v_LISM;
		    }
	    }
	}
	
	pfile.close();
}

void init_hydrogen_header(const char* fname, vel_dist_header& h)
{
	int myrank;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

	if(!myrank) {
		int    tidata[3];
		double tddata[3];
		int    nx[3] = {0, 0, 0}; // {nxr, nxt, nxp}
		int    nv[3] = {0, 0, 0}; // to read the header, these vals don't matter yet

		ifstream vfile(fname, ios::in|ios::binary|ios::ate);
		// check if file is opened correctly, exit if not
		cout << endl << fname << endl;
		if(vfile.good()) {
			cout << " File open successfully\n";
		} else {
			cout << " File open failed\n";
			return; // will cause the program to hang, aka crash
		}
		streampos file_size = vfile.tellg(); // get size of file in bytes
		cout << "Size of file = " << file_size/1024.0/1024.0 << " MB\n\n";
		vfile.seekg(0, ios::beg); // start at beginning of file

		// Read header from binary file
		ff_bskipRead(vfile, 1, 1);
		ff_bdataRead(vfile, 1, 3, 2, tidata, 0, 0, 0, nx);
		h.nxr = tidata[0];
		h.nxt = tidata[1];
		h.nxp = tidata[2];
		ff_bskipRead(vfile, 1, 2);
		ff_bdataRead(vfile, 3, 3, 2, tddata, 0, 0, 0, nx);
		h.xalpha = tddata[0];
		h.xmin   = tddata[1] / 1.5E11;//AU_to_M;
		h.xmax   = tddata[2] / 1.5E11;//AU_to_M;
		ff_bskipRead(vfile, 1, 2);
		ff_bdataRead(vfile, 1, 3, 2, tidata, 0, 0, 0, nv);
		h.nvr = tidata[0];
		h.nvt = tidata[1];
		h.nvp = tidata[2];
		ff_bskipRead(vfile, 1, 2);
		ff_bdataRead(vfile, 3, 3, 2, tddata, 0, 0, 0, nv);
		h.valpha = tddata[0];
		h.vmin   = tddata[1]; // cm/s
		h.vmax   = tddata[2]; // cm/s
		ff_bskipRead(vfile, 1, 2);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 1);
		vfile.close();
	}
	MPI_Bcast(&h, 1, mpi_vel_dist_header, ROOT, MPI_COMM_WORLD);
}

// hard coded in header !
void init_proton_header(vel_dist_header& ph)
{
	// init header for 3D plasma dist
	ph.nxr    = 380;
	ph.nxt    = 120;
	ph.nxp    = 118;
	ph.xalpha = 3.0;
	ph.xmin   = 10.0; // in AU
	ph.xmax   = 1000.0; // in AU
	// set velocity space to nil
	ph.nvr = ph.nvt = ph.nvp = 0;
	ph.valpha = 0.0;
	ph.vmin = ph.vmax = 0.0;
}

void compute_moments_hydrogen(const char* vel_dist_fname, vel_dist_header& h,
 	double *dens, double *velx, double *vely, double *velz,
 	double *vxvar, double *vyvar, double *vzvar,
 	double *vxn, double *vyn, double *vzn, double *d3vn)
{
	int myrank, wrank, nrank, m_msg, w_msg;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &nrank);
	nrank--; // number of workers total
	if(nrank < 1) {
		cout << "\nERROR: NEED AT LEAST 2 THREADS TO RUN\n\n";
		return;
	}

	int ixr, ixt, ixp, ivr, ivt, ivp, idx, idv, w_idx;
	int nv[3] = {h.nvr, h.nvt, h.nvp};
	int nx[3] = {h.nxr, h.nxt, h.nxp};
	int vcells = nv[0] * nv[1] * nv[2];

	double vtemp[3], rmean[3], wd3v, tdens;
	float *dist_chunk;
	double *moment_packet;
	moment_packet = new double[8];
	dist_chunk    = new float[vcells];

	m_msg = 1; // there is still work
	if(!myrank) { // master thread
		timer loop_timer;
		wrank = 1;
		// open file and skip header
		ifstream vfile(vel_dist_fname, ios::in|ios::binary);
		ff_bskipRead(vfile, 1, 6);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 7);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 2);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 1);
		// loop over all positions as we read through the file
		start_timer(loop_timer);
		for(ixr = 0; ixr < h.nxr; ixr++) {
		  for(ixt = 0; ixt < h.nxt; ixt++) {
		    for(ixp = 0; ixp < h.nxp; ixp++) {
		    	//cout << wrank << endl;
		      idx = arrayIdx3(ixr, ixt, ixp, nx); // index into neutral grid
		        //cout << "index: " << idx << endl;
		      // Loop over velocity space
		      for(ivr = 0; ivr < h.nvr; ivr++) {
		      	for(ivt = 0; ivt < h.nvt; ivt++) {
		      	  // Read next chunk from bin file (a strip of phi)
			      ff_bskipRead(vfile, 1, 1);
			      ff_bdataRead(vfile, 2, h.nvp, 2, dist_chunk, ivr, ivt, 0, nv);
			      ff_bskipRead(vfile, 1, 1);
		      	}
		      }
		      // Finished with reading vel space chunk, now let's send it to next worker
		      // First, revc the msg from the worker so we know if there is data to get
		      //   or not, then let them know we are going to send work
		      MPI_Recv(&w_msg, 1, MPI_INT, wrank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      MPI_Send(&m_msg, 1, MPI_INT, wrank, 2, MPI_COMM_WORLD);

		      if(w_msg == 1) {
		      	//MPI_Recv(&w_idx, 1, MPI_INT, wrank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      	MPI_Recv(moment_packet, 8, MPI_DOUBLE, wrank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      	// store recved data at w_idx
		      	w_idx = moment_packet[7];
		      	dens[w_idx]  = moment_packet[0];
		      	velx[w_idx]  = moment_packet[1];
		      	vely[w_idx]  = moment_packet[2];
		      	velz[w_idx]  = moment_packet[3];
		      	vxvar[w_idx] = moment_packet[4];
		      	vyvar[w_idx] = moment_packet[5];
		      	vzvar[w_idx] = moment_packet[6];
		      }
		      // send the 3D vel space data and the index
		      w_idx = idx;
		      MPI_Send(&w_idx, 1, MPI_INT, wrank, 3, MPI_COMM_WORLD);
		      MPI_Send(dist_chunk, vcells, MPI_FLOAT, wrank, 0, MPI_COMM_WORLD);

		      // round robin increment worker
		      wrank = (wrank >= nrank ? 1 : wrank+1);
		    }
		  }
		  current_time(loop_timer);
		  cout << 100.0 * (ixr + 1.0) / double(h.nxr) << "% done\n";
		  display_remaining_time(loop_timer, ixr+1.0, h.nxr, 2);
		}
		vfile.close();
		m_msg = 0;
		// loop through all workers once to let them all know there is no more work
		//  but, make sure to recv the last bits of data
		for(int w_exit = 0; w_exit < nrank; w_exit++) {
			printf("send finished message to %i\n", wrank);
			MPI_Recv(&w_msg, 1, MPI_INT, wrank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    MPI_Send(&m_msg, 1, MPI_INT, wrank, 2, MPI_COMM_WORLD);

		    if(w_msg == 1) {
		      	//MPI_Recv(&w_idx, 1, MPI_INT, wrank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      	MPI_Recv(moment_packet, 8, MPI_DOUBLE, wrank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		      	// store recved data at w_idx
		      	w_idx = moment_packet[7];
		      	dens[w_idx]  = moment_packet[0];
		      	velx[w_idx]  = moment_packet[1];
		      	vely[w_idx]  = moment_packet[2];
		      	velz[w_idx]  = moment_packet[3];
		      	vxvar[w_idx] = moment_packet[4];
		      	vyvar[w_idx] = moment_packet[5];
		      	vzvar[w_idx] = moment_packet[6];
		    }
		    wrank = (wrank >= nrank ? 1 : wrank+1);
		}
		cout << "master thread finished hydrogen momemts\n";
	} else { // worker thread
		w_msg = 0; // no work finished yet
		while(m_msg == 1) { // still work to recv
			MPI_Send(&w_msg, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			MPI_Recv(&m_msg, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if(w_msg == 1) { // finished computation to send from prev loop
				//MPI_Send(&w_idx, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
				moment_packet[7] = w_idx;
				MPI_Send(moment_packet, 8, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
			}
			if(m_msg == 1) { // new data to work on
				MPI_Recv(&w_idx, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(dist_chunk, vcells, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				memset(moment_packet, 0, 8*sizeof(double));
				// loop through vel space, computing moments
				for(ivr = 0; ivr < h.nvr; ivr++) {
				  for(ivt = 0; ivt < h.nvt; ivt++) {
				  	for(ivp = 0; ivp < h.nvp; ivp++) {
				  	  idv  = arrayIdx3(ivr, ivt, ivp, nv);
				  	  wd3v = dist_chunk[idv] * d3vn[idv];

				  	  if(wd3v > 1.0E-32) { // avoid div by zero
				  	  	// prestore first bulk speed values
						if((ivp + ivt + ivr)== 0) {
						  moment_packet[1] = vxn[0];
						  moment_packet[2] = vyn[0];
						  moment_packet[3] = vzn[0];
						}
						// online variance and mean algorithm (West 1979)
				  		vtemp[0] = vxn[idv] - moment_packet[1];
				  		vtemp[1] = vyn[idv] - moment_packet[2];
				  		vtemp[2] = vzn[idv] - moment_packet[3];
				  		
				  		tdens = moment_packet[0] + wd3v;
				  		
				  		rmean[0] = vtemp[0] * wd3v / tdens;
				  		rmean[1] = vtemp[1] * wd3v / tdens;
				  		rmean[2] = vtemp[2] * wd3v / tdens;
				  		
				  		moment_packet[1] += rmean[0];
				  		moment_packet[2] += rmean[1];
				  		moment_packet[3] += rmean[2];
				  		
				  		moment_packet[4] += rmean[0] * moment_packet[0] * vtemp[0];
				  		moment_packet[5] += rmean[1] * moment_packet[0] * vtemp[1];
				  		moment_packet[6] += rmean[2] * moment_packet[0] * vtemp[2];
				  		
				  		moment_packet[0] = tdens;
				  	  }
				  	}
				  }
				}// Finish vel space loop
				// Normalize to density (cancels out d3x)
				// [velocity] = [km/s]
				moment_packet[1] *= 1.0E-3;
			    moment_packet[2] *= 1.0E-3;
			    moment_packet[3] *= 1.0E-3;
			    moment_packet[4] *= 1.0E-6* vcells / ((vcells-1.0)*moment_packet[0]);
			    moment_packet[5] *= 1.0E-6* vcells / ((vcells-1.0)*moment_packet[0]);
			    moment_packet[6] *= 1.0E-6* vcells / ((vcells-1.0)*moment_packet[0]);
			    // [density] = [1/cm^3]
			    moment_packet[0]  *= 1.0E-6;
				w_msg = 1;
			}
		} // END WHILE (of still work to-do)
		printf("worker %i/%i finished hydrogen moments\n", myrank, nrank);
	}
	delete[] moment_packet, dist_chunk;
}

void interp_h2p_grid(double *data, double *data_2p, vel_dist_header &h, vel_dist_header &ph,
	grid_cell_sph &gc, grid_cell_sph &pgc)
{
	int ixr, ixt, ixp, idx, idx2, myrank, nrank, irank, interp_case;
	int ir0, it0, ip0, idx8[8];
	int pnx[3] = {ph.nxr, ph.nxt, ph.nxp};
	int nx[3] = {h.nxr, h.nxt, h.nxp}, nxcells = h.nxr * h.nxt * h.nxp;
	timer loop_timer;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nrank);
	int offset[nrank+1];

	printf("Enter interp: thread %i/%i\n", myrank, nrank);

	for(irank = 0; irank < nrank+1; irank++) {
		offset[irank] = (irank  * ph.nxr) / nrank;
		//if(!myrank) printf("%i offset: %i/%i\n", irank, offset[irank],ph.nxr);
	}

	
	for(ixr = offset[myrank]; ixr < offset[myrank+1]; ixr++) {
		ir0 = 0;
	  if(!myrank) start_timer(loop_timer);
	  // Find left most neutral index next to the current plasma index position for rho
	  if(pgc.cxr[ixr] < h.xmin) {        // no z-axis
	  	interp_case = 1;
	  	ir0 = 0;
	  } else if(pgc.gxr[ixr] > gc.gxr[nx[0]-2]){  // outside domain
	  	interp_case = 2;
	  	ir0 = h.nxr - 2;
	  } else if(gc.gxr[0] == 0.0 &&       // near z-axis 
	  			pgc.gxr[ixr] < gc.gxr[1] )
	  {
	  	interp_case = 0;
	  	ir0 = 0;
	  } else {                            // inside domain
	  	interp_case = 3;

	  	while(!(pgc.gxr[ixr] > gc.gxr[ir0] && pgc.gxr[ixr] <= gc.gxr[ir0+1] )) ir0++;
	  	//if(!myrank)
	  	//	printf("%f < %f <= %f\n",gc.gxr[ir0], pgc.gxr[ixr], gc.gxr[ir0+1]  );
	  }
	  
	  for(ixt = 0; ixt < ph.nxt; ixt++) {
	  	it0 = 0;
	  	// Find left most neutral index next to the current plasma index position for theta
	  	while(!(pgc.gxt[ixt] > gc.cxt[it0] && pgc.gxt[ixt] <= gc.cxt[it0+1] )) it0++;
	  	//if(!myrank && ixr == 0)
	  	//	printf("%f < %f <= %f\n", gc.gxt[it0], pgc.gxt[ixt], gc.gxt[it0+1]);
	  	
	    for(ixp = 0; ixp < ph.nxp; ixp++) {
	    	ip0 = 0;
	      // Find left most neutral index next to the current plasma index position for phi
	  	  while(!(pgc.cxp[ixp] > gc.cxp[ip0] && pgc.cxp[ixp] <= gc.cxp[ip0+1] )) ip0++;
		  //if(!myrank && ixr == 0 && ixt == 0)
	  	  //printf("%f < %f <= %f\n", gc.gxp[ip0], pgc.gxp[ixp], gc.gxp[ip0+1]);
	  	
	      idx = arrayIdx3(ixr, ixt, ixp, pnx);

	      // if(!myrank) {
	      // 	printf("Proton idx: %f %f %f\n Neutral idx: %f %f %f\n",
	      // 	  pgc.cxr[ixr], pgc.cxt[ixt], pgc.cxp[ixp], gc.cxr[ir0], gc.cxt[it0], gc.cxp[ip0]);
	      // }

	      switch(interp_case) {
	      	// r0 < r < r1 (near z axis), petrahedron (aka wedge) interp
	        // we will approx the tetra and pyramid cases as a wedge for simplicity
	      	case 0:
	        //  break;
	      	// r < r0 (no z axis), rectangle interp
	      	case 1: // fall through to case 2
	      	// r > rn (outside of domain), rectangle interp
	      	case 2:
	      	  idx8[0] = arrayIdx3(ir0  , it0  , ip0  , nx);
	      	  idx8[1] = arrayIdx3(ir0  , it0  , ip0+1, nx);
	      	  idx8[2] = arrayIdx3(ir0  , it0+1, ip0  , nx);
	      	  idx8[3] = arrayIdx3(ir0  , it0+1, ip0+1, nx);

	      	  if(idx8[3] > nxcells)
	      	  	printf("ERROR: idx4 %i/%i out of bounds in thread %i \n", idx8[3], nxcells, myrank);

	      	  data_2p[idx] =  int(idx8[3] > nxcells);// myrank + 0.5;
	      	   // ((gc.gxt[it0+1] - pgc.gxt[ixt])
	      	   //  *  ((gc.gxp[ip0+1] - pgc.gxp[ixp])*data[idx8[0]] + (pgc.gxp[ixp] - gc.gxp[ip0])*data[idx8[1]] )
	      	   //  +  (pgc.gxt[ixt] - gc.gxt[it0])

	      	   //  *  ((gc.gxp[ip0+1] - pgc.gxp[ixp])*data[idx8[2]] + (pgc.gxp[ixp] - gc.gxp[ip0])*data[idx8[3]] ))
	      	   //  / ((gc.gxt[it0+1]-gc.gxt[it0]) * (gc.gxp[ip0+1]-gc.gxp[ip0]));
	      	  break;

	      	case 3:
	      	  idx8[0] = arrayIdx3(ir0  , it0  , ip0  , nx);
	      	  idx8[1] = arrayIdx3(ir0  , it0  , ip0+1, nx);
	      	  idx8[2] = arrayIdx3(ir0  , it0+1, ip0  , nx);
	      	  idx8[3] = arrayIdx3(ir0  , it0+1, ip0+1, nx);
	      	  idx8[4] = arrayIdx3(ir0+1, it0  , ip0  , nx);
	      	  idx8[5] = arrayIdx3(ir0+1, it0  , ip0+1, nx);
	      	  idx8[6] = arrayIdx3(ir0+1, it0+1, ip0  , nx);
	      	  idx8[7] = arrayIdx3(ir0+1, it0+1, ip0+1, nx);

	      	  if(idx8[7] > nxcells)
	      	  	printf("ERROR: idx8 %i/%i out of bounds in thread %i \n", idx8[7], nxcells, myrank);

	      	  data_2p[idx] = 2*int(idx8[7] > nxcells);
	      	   // ((gc.gxr[ir0+1] - pgc.gxr[ixr])
	      	   //  *  ((gc.gxt[it0+1] - pgc.gxt[ixt])
	      	   //  *   ((gc.gxp[ip0+1] - pgc.gxp[ixp])*data[idx8[0]] + (pgc.gxp[ixp] - gc.gxp[ip0])*data[idx8[1]])
	      	   //  +   (pgc.gxt[ixt] - gc.gxt[it0])
	      	   //  *   ((gc.gxp[ip0+1] - pgc.gxp[ixp])*data[idx8[2]] + (pgc.gxp[ixp] - gc.gxp[ip0])*data[idx8[3]]))

	      	   //  +  (pgc.gxr[ixr] - gc.gxr[ir0])
	      	   //  *  ((gc.gxt[it0+1] - pgc.gxt[ixt])
	      	   //  *   ((gc.gxp[ip0+1] - pgc.gxp[ixp])*data[idx8[4]] + (pgc.gxp[ixp] - gc.gxp[ip0])*data[idx8[5]])
	      	   //  +   (pgc.gxt[ixt] - gc.gxt[it0])
	      	   //  *   ((gc.gxp[ip0+1] - pgc.gxp[ixp])*data[idx8[6]] + (pgc.gxp[ixp] - gc.gxp[ip0])*data[idx8[7]])))
	      	   //  / ((gc.gxr[ir0+1]-gc.gxr[ir0]) * (gc.gxt[it0+1]-gc.gxt[it0]) * (gc.gxp[ip0+1]-gc.gxp[ip0]));
	      	  break;
	      	default: break;
	      }
	      
	    }
	  }
	  if(!myrank && (ixr % 10 == 0)) {
	  	current_time(loop_timer);
	  	cout << 100.0 * (ixr + 1.0) / double(offset[1]) << "% done\n";
		display_remaining_time(loop_timer, ixr+1.0, offset[1], 2);
	  }
	}
	printf("Thread %i finished interp\n", myrank);

	// MPI_Allreduce(
 //    void* send_data,
 //    void* recv_data,
 //    int count,
 //    MPI_Datatype datatype,
 //    MPI_Op op,
 //    MPI_Comm communicator)

	//send all missing chucks to all other processes
	for(irank = 0; irank < nrank; irank++) {
		idx  = arrayIdx3(offset[irank], 0, 0, pnx);
		idx2 = arrayIdx3(offset[irank+1], 0, 0, pnx);
		//MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(data_2p + idx, idx2-idx, MPI_DOUBLE, irank, MPI_COMM_WORLD);
		if(irank == myrank) printf("Thread %i sending broacast of size %i\n", myrank, idx2-idx);
	}
}


// mpirun -np 8 ./rdmpi
//***********************************************************************************//
int main (int argc, char* argv[])
{
	double *VPTR = NULL;
	timer main_timer, vtk_timer, plasma_timer, source_timer;
	// MPI startup
	int rank, size;
	MPI_Init (&argc, &argv);      /* starts MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);        /* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);        /* get number of processes */

	// define mpi_struct type for sending the header data from master to slaves
	int blocks[12]={1,1,1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype types[12]={MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint intex, doubex;
	MPI_Type_extent(MPI_INT, &intex);
	MPI_Type_extent(MPI_DOUBLE, &doubex);
	MPI_Aint displacements[12] = {0, intex, 2*intex, 3*intex, 4*intex, 5*intex,
		6*intex, 6*intex+doubex, 6*intex+2*doubex, 6*intex+3*doubex, 6*intex+4*doubex,
		6*intex+5*doubex};
	MPI_Type_struct(12, blocks, displacements, types, &mpi_vel_dist_header);
	MPI_Type_commit(&mpi_vel_dist_header);
	MPI_Barrier(MPI_COMM_WORLD);

	// indices
	int ixr, ixt, ixp, ivr, ivt, ivp, idx, idv;
	// temp vars for intermediate calculations
	double dvr, dvt, dvp, vx, vy, vz, d3v, dv, vxy;

	const char vel_dist_fname[] = "./circ2015_3muG_HP120B/vel_dist_circ2015_3muG_HP120B_FULL.grid";
	const char plasma_prop_fname[] = "./circ2015_3muG_HP120B/plasma_circ2015_3muG_HP120B.plasma";
	// alloc headers and grid for hydrogen and protons
	vel_dist_header h, ph;
	grid_cell_sph gc, pgc;

	start_timer(main_timer);

	init_hydrogen_header(vel_dist_fname, h);
	if(!rank) {
		cout << "Neutral distribution header:\n";
		display_header(h);
	}
	init_grid_cell_sph(gc, h);
	// if(!rank)
	// 	display_grid_cell_sph(gc, h);

	init_proton_header(ph);
	if(!rank) {
		cout << "\nPlasma distribution header:\n";
		display_header(ph);
	}
	init_grid_cell_sph(pgc, ph);
	// if(!rank)
	// 	display_grid_cell_sph(pgc, ph);

	//**********************************************************************************//
	// read plasma data from binary file
	int pxcells = ph.nxr * ph.nxt * ph.nxp;
	double *pdens, *ptemp, *pvelx, *pvely, *pvelz;
	if(READ_PLASMA_FILE) {
		pdens  = new double[pxcells]; memset(pdens, 0, pxcells*sizeof(double));
		pvelx  = new double[pxcells]; memset(pvelx, 0, pxcells*sizeof(double));
		pvely  = new double[pxcells]; memset(pvely, 0, pxcells*sizeof(double));
		pvelz  = new double[pxcells]; memset(pvelz, 0, pxcells*sizeof(double));
		ptemp  = new double[pxcells]; memset(ptemp, 0, pxcells*sizeof(double));
	}
	
	start_timer(plasma_timer);
	if(READ_PLASMA_FILE) {
		if(!rank) {
			cout << "READ_PLASMA_FILE\n\n";
			read_plasma_binary(plasma_prop_fname, ph,
				pdens, pvelx, pvely, pvelz, ptemp);
		}
		MPI_Bcast(pdens, pxcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(pvelx, pxcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(pvely, pxcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(pvelz, pxcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		MPI_Bcast(ptemp, pxcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	}
	current_time(plasma_timer);

	//**********************************************************************************//
	// Precalculate the velocities and d3v in neutral vel space
	int xcells = h.nxr * h.nxt * h.nxp;
	int vcells = h.nvr * h.nvt * h.nvp;
	int nx[3] = {h.nxr, h.nxt, h.nxp};
	int nv[3] = {h.nvr, h.nvt, h.nvp};
	double *vxn, *vyn, *vzn, *d3vn;
	d3vn = new double[vcells];
	vxn  = new double[vcells]; memset(vxn , 0, vcells*sizeof(double));
	vyn  = new double[vcells]; memset(vyn , 0, vcells*sizeof(double));
	vzn  = new double[vcells]; memset(vzn , 0, vcells*sizeof(double));

	dvt = gc.gvt[1] - gc.gvt[0];
	dvp = gc.gvp[1] - gc.gvp[0];
	for(ivr = 0; ivr < h.nvr; ivr++) {
	  dvr = gc.gvr[ivr+1] - gc.gvr[ivr];
	  dv  = dvr * dvp * dvt;
	  for(ivt = 0; ivt < h.nvt; ivt++) {
	    d3v = gc.cvr[ivr] * gc.cvr[ivr] * sin(gc.cvt[ivt]) * dv;
	    for(ivp = 0; ivp < h.nvp; ivp++) {
	      idv = arrayIdx3(ivr, ivt, ivp, nv);
	      // in units of m/s
	      vz  = gc.cvr[ivr] * cos(gc.cvt[ivt]);
		  vxy = gc.cvr[ivr] * sin(gc.cvt[ivt]);
		  vx  = vxy * cos(gc.cvp[ivp]);
		  vy  = vxy * sin(gc.cvp[ivp]);
		  
		  vxn[idv]  = vx;
		  vyn[idv]  = vy;
		  vzn[idv]  = vz;
		  d3vn[idv] = d3v;
	    }
	  }
	}

	double *dens, *velx, *vely, *velz, *vxvar, *vyvar, *vzvar;
	if(READ_NEUTRAL_FILE_MOMENTS) {
		dens  = new double[xcells]; memset(dens, 0, xcells*sizeof(double));
		if(!rank) {
			velx  = new double[xcells]; memset(velx, 0, xcells*sizeof(double));
			vely  = new double[xcells]; memset(vely, 0, xcells*sizeof(double));
			velz  = new double[xcells]; memset(velz, 0, xcells*sizeof(double));
			vxvar = new double[xcells]; memset(vxvar, 0, xcells*sizeof(double));
			vyvar = new double[xcells]; memset(vyvar, 0, xcells*sizeof(double));
			vzvar = new double[xcells]; memset(vzvar, 0, xcells*sizeof(double));
		}
		compute_moments_hydrogen(vel_dist_fname, h,
			dens, velx, vely, velz, vxvar, vyvar, vzvar, vxn, vyn, vzn, d3vn);
		MPI_Bcast(dens, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD); // for interp grid later
		//MPI_Bcast(velx, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		//MPI_Bcast(vely, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		//MPI_Bcast(velz, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		//MPI_Bcast(vxvar, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		//MPI_Bcast(vyvar, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
		//MPI_Bcast(vzvar, xcells, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
	}
		

	//*******************************************************************************
	// Calculate source & loss terms (flux to protons from hydrgoen)
	// For this algorithm, we will first interpolate the neutral density data onto the
	//  plasma grid. There are many ways to do this, so first we will try the easiest,
	//  but maybe not the fastest, which will be related to the inverse of the distances
	// We will also need the density to be precalculated, so for now require the moments
	//  to be calculated.
	double *dens_2p;
	if(!READ_NEUTRAL_FILE_MOMENTS && COMPUTE_SOURCE_TERMS) {
		if(!rank)
			cout << "To compute source terms set READ_NEUTRAL_FILE_MOMENTS = 1\n";
	} else if(READ_NEUTRAL_FILE_MOMENTS && COMPUTE_SOURCE_TERMS){
		dens_2p = new double[pxcells];
		memset(dens_2p, 0, pxcells*sizeof(double));
		printf("Before interp: thread %i\n", rank);
		interp_h2p_grid(dens, dens_2p, h, ph, gc, pgc);
		//compute_source_terms(h, ph, dens_2p, pdens, pvelx, pvely, pvelz, ptemp
		//	vxn, vyn, vzn, d3vn, fdens, fvelx, fvely, fvelz, fenrg);
	}
		

	//*******************************************************************************
	// Print plasma data and neutral moments to vtk file
	start_timer(vtk_timer);
	if(!rank && PRINT_PLASMA_VTK) {
		print2vtk(
		  "./vtk/plasma_circ2015_3muG_HP120B.vtk",
		  ph, pgc, pdens, pvelx, pvely, pvelz, ptemp, VPTR, VPTR, VPTR, 31);
	}
	if(!rank && PRINT_NEUTRAL_MOMENTS_VTK) {
		print2vtk(
		  "./vtk/vel_dist_circ2015_3muG_HP120B_FULL.vtk",
		  h, gc, dens, velx, vely, velz, vxvar, vyvar, vzvar, VPTR, 127);
	}
	if(!rank && PRINT_SOURCE_TERMS) {
		print2vtk(
		  "./vtk/vel_dist_circ2015_3muG_HP120B_FULL_interp_dens_only.vtk",
		  ph, pgc, dens_2p, VPTR, VPTR, VPTR, VPTR, VPTR, VPTR, VPTR, 1);
		// print2vtk(
		//   "./vtk/plasma_circ2015_3muG_HP120B_source_terms_maxwellian_total_sqrt.vtk",
		//   ph, pgc, fdens, fvelx, fvely, fvelz, VPTR, VPTR, VPTR, fenrg, 143);
	}
		
	current_time(vtk_timer);


	if(!rank) {
		current_time(main_timer);
		cout << "\n Main timer:\n";
		display_elapse_time(main_timer, 3);
		cout << " Plasma reader timer:\n";
		display_elapse_time(plasma_timer);
		if(PRINT_PLASMA_VTK) {
			cout << " VTK printer timer:\n";
			display_elapse_time(vtk_timer);
		}
	}

	delete_grid_cell_sph(gc);
	delete_grid_cell_sph(pgc);
	MPI_Type_free(&mpi_vel_dist_header);
	delete[] d3vn, vxn, vyn, vzn;
	if(READ_NEUTRAL_FILE_MOMENTS && COMPUTE_SOURCE_TERMS)
		delete[] dens_2p;
	if(READ_PLASMA_FILE)
		delete[] pdens, pvelx, pvely, pvelz, ptemp;
	if(READ_NEUTRAL_FILE_MOMENTS) {
		delete[] dens;
		if(!rank)
			delete[] velx, vely, velz, vxvar, vyvar, vzvar;
	}
		
	MPI_Finalize();
	return 0;
}