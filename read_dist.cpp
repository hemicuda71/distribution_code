// to compile:
// g++ -O3 read_dist.cpp timer.cpp ./special_funcs/specialfunctions.cpp ./special_funcs/alglibinternal.cpp ./special_funcs/ap.cpp 


#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string.h> // for memset
#include "timer.h"
#include "./special_funcs/specialfunctions.h" // for erf function

using namespace std;

#define READ_NEUTRAL_FILE_MOMENTS 1
#define READ_NEUTRAL_FILE_LOSS_TERMS 1
#define READ_PLASMA_FILE 1
#define READ_OUTPUT_SOURCE 1

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
	// only calculate the hard tail if high enough energy (2nd term at 10%)
	double tail = (rel_eng < 17.8698  ? 1.0 : pow(1.0 - exp(-a3/rel_eng), 4.5));
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


void print2vtk_v2_vspace(const char* vtkfname, vel_dist_header& h, grid_cell_sph& gc,
	double *scalar_field[], double *vector_field[], const char *field_name[],
	int NSF, int NVF)
{
	int ixr, ixt, ixp, dk, NUM_CELL_VALS, i, j, k, idx;
	int cell_type, cell_iter = 0;
	float x, y, z, xy;
	int POINTS = (h.nvr+1)*(h.nvt+1)*(h.nvp+1);
	int CELLS  = h.nvr * h.nvt * h.nvp;
	int *cell_type_array;
	cell_type_array = new int[CELLS];
	int nxp[3] = {h.nvr+1, h.nvt+1, h.nvp+1}; // these 1's are very important!!!
	int nxc[3] = {h.nvr, h.nvt, h.nvp};
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
	for(ixr = 0; ixr <= h.nvr; ixr++)  // rho
		for(ixt = 0; ixt <= h.nvt; ixt++)  // theta
			for(ixp = 0; ixp <= h.nvp; ixp++) { // phi
				//sph2cart(rgrid[k], tgrid[j], pgrid[i], x, y, z);
				z  = gc.gvr[ixr] * cos(gc.gvt[ixt]);
		  		xy = gc.gvr[ixr] * sin(gc.gvt[ixt]);
		  		x  = xy * cos(gc.gvp[ixp]);
		  		y  = xy * sin(gc.gvp[ixp]);
				fvtk << x << ' ' << y << ' ' << z << endl;
			}
	cout << " POINTS finished\n";
	
	// get number of cell values
	NUM_CELL_VALS = 0;
	for(ixr = 0; ixr < h.nvr; ixr++) { // rho
		for(ixt = 0; ixt < h.nvt; ixt++) { // theta
			for(ixp = 0; ixp < h.nvp; ixp++) { // phi
				// At origin and + z-axis
				if( gc.gvr[ixr] < EPS && abs(gc.gvt[ixt]) < EPS ) {
					NUM_CELL_VALS += 5; // tetra
				// At origin and - z-axis
				} else if( gc.gvr[ixr] < EPS && abs(gc.gvt[ixt+1] - M_PI) < EPS ) {
					NUM_CELL_VALS += 5; // tetra
				// At origin, but not z-axis
				} else if( gc.gvr[ixr] < EPS ) {
					NUM_CELL_VALS += 6; // pyramid
				// At + z-axis, but not origin
				} else if( abs(gc.gvt[ixt]) < EPS ) {
					NUM_CELL_VALS += 7; // wedge
				// At - z-axis, but not origin
				} else if( abs(gc.gvt[ixt+1] - M_PI) < EPS ) {
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
	for(i = 0; i < h.nvr; i++) { // rho
		for(j = 0; j < h.nvt; j++) { // theta
			for(k = 0; k < h.nvp; k++) { // phi
				// we want phi to be connected at the last cell
				dk = (k == h.nvp-1 ? 0 : k+1);
				
				// At origin and + z-axis
				if( gc.gvr[i] < EPS && abs(gc.gvt[j]) < EPS ) {
					cell_type = 10; // tetra
					fvtk << 4;
					fvtk << ' ' << arrayIdx3(i, 0, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);
				// At origin and - z-axis
				} else if( gc.gvr[i] < EPS && abs(gc.gvt[j+1] - M_PI) < EPS ) {
					cell_type = 10; // tetra
					fvtk << 4;
					fvtk << ' ' << arrayIdx3(i, 0, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, 0, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, dk, nxp);
				// At origin, but not z-axis
				} else if( gc.gvr[i] < EPS ) {
					cell_type = 14; // pyramid
					fvtk << 5;
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, dk, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);
					fvtk << ' ' << arrayIdx3(i, 0, 0, nxp);
				// At + z-axis, but not origin
				} else if( abs(gc.gvt[j]) < EPS ) {
					cell_type = 13; // wedge
					fvtk << 6;
					fvtk << ' ' << arrayIdx3(i+1, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j+1, dk, nxp);
					fvtk << ' ' << arrayIdx3(i+1, j, 0, nxp);
					fvtk << ' ' << arrayIdx3(i, j+1, k, nxp);
					fvtk << ' ' << arrayIdx3(i, j+1, dk, nxp);
					fvtk << ' ' << arrayIdx3(i, j, 0, nxp);
				// At - z-axis, but not origin
				} else if( abs(gc.gvt[j+1] - M_PI) < EPS ) {
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
	// Scalar data
	fvtk << "\nCELL_DATA " << CELLS << endl;
	for(i = 0; i < NSF; i++) {
		fvtk << "SCALARS " << field_name[i] << " FLOAT\n";
		fvtk << "LOOKUP_TABLE default\n";
		for(ixr = 0; ixr < h.nvr; ixr++) { // rho
			for(ixt = 0; ixt < h.nvt; ixt++) { // theta
				for(ixp = 0; ixp < h.nvp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << scalar_field[i][idx] << endl;
				}
			}
		}
		cout << ' ' <<  field_name[i] << " data finished\n";
	}
	// Vector data
	for(i = 0; i < NVF; i++) {
		fvtk << "VECTORS " << field_name[NSF+i] << " FLOAT\n";
		for(ixr = 0; ixr < h.nvr; ixr++) { // rho
			for(ixt = 0; ixt < h.nvt; ixt++) { // theta
				for(ixp = 0; ixp < h.nvp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << vector_field[3*i  ][idx] << ' '
						 << vector_field[3*i+1][idx] << ' '
						 << vector_field[3*i+2][idx] << endl;
				}
			}
		}
		cout << ' ' <<  field_name[NSF+i] << " data finished\n";
	}
	cout << vtkfname << " printing success!\n\n";
	fvtk.close();
}

void print2vtk_v2(const char* vtkfname, vel_dist_header& h, grid_cell_sph& gc,
	double *scalar_field[], double *vector_field[], const char *field_name[],
	int NSF, int NVF)
{
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
		  		x  = xy * cos(gc.gvp[ixp]);
		  		y  = xy * sin(gc.gvp[ixp]);
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
	// Scalar data
	fvtk << "\nCELL_DATA " << CELLS << endl;
	for(i = 0; i < NSF; i++) {
		fvtk << "SCALARS " << field_name[i] << " FLOAT\n";
		fvtk << "LOOKUP_TABLE default\n";
		for(ixr = 0; ixr < h.nxr; ixr++) { // rho
			for(ixt = 0; ixt < h.nxt; ixt++) { // theta
				for(ixp = 0; ixp < h.nxp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << scalar_field[i][idx] << endl;
				}
			}
		}
		cout << ' ' <<  field_name[i] << " data finished\n";
	}
	// Vector data
	for(i = 0; i < NVF; i++) {
		fvtk << "VECTORS " << field_name[NSF+i] << " FLOAT\n";
		for(ixr = 0; ixr < h.nxr; ixr++) { // rho
			for(ixt = 0; ixt < h.nxt; ixt++) { // theta
				for(ixp = 0; ixp < h.nxp; ixp++) { // phi
					idx = arrayIdx3(ixr, ixt, ixp, nxc);
					fvtk << vector_field[3*i  ][idx] << ' '
						 << vector_field[3*i+1][idx] << ' '
						 << vector_field[3*i+2][idx] << endl;
				}
			}
		}
		cout << ' ' <<  field_name[NSF+i] << " data finished\n";
	}
	cout << vtkfname << " printing success!\n\n";
	fvtk.close();
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

// void interp_h2p_grid(double *data, double *data_2p, vel_dist_header &h, vel_dist_header &ph,
// 	grid_cell_sph &gc, grid_cell_sph &pgc)
// {
// 	int ixr, ixt, ixp, idx, idx2, interp_case;
// 	int ir0, it0, ip0, idx8[8];
// 	int pnx[3] = {ph.nxr, ph.nxt, ph.nxp};
// 	int nx[3] = {h.nxr, h.nxt, h.nxp}, nxcells = h.nxr * h.nxt * h.nxp;
// 	timer loop_timer;
// 	start_timer(loop_timer);

// 	for(ixr = 0; ixr < ph.nxr; ixr++) {
// 	  ir0 = 0;
// 	  interp_case = 0;
// 	  // Find left most neutral index next to the current plasma index position for rho
// 	    //  interp_case = xx00 r in domain
// 	    //              = xx01 ixr == 0 w/ origin
// 	    //              = xx10 ixr == 0 w/o origin
// 	    //              = xx11 ixr == nxr
// 	    //              = 00xx theta in domain
// 	    //              = 01xx ixt == 0
// 	    //              = 10xx ixt == nxt

// 	  if(gc.gxr[0] <= EPS && pgc.cxr[ixr] <= gc.cxr[0]) { // origin, inner boundary
// 	  	interp_case += 1;
// 	  } else if(gc.gxr[0] > EPS && pgc.cxr[ixr] <= gc.cxr[0]) { // no origin, inner boundary
// 	  	interp_case += 2;
// 	  } else if(pgc.cxr[ixr] >= gc.cxr[h.nxr-1]) { // outer boundary
// 	  	interp_case += 3;
// 	  } // else, rho inside domain (already 0)

// 	  if(pgc.cxr[ixr] > gc.cxr[0])
// 	  while((ir0 < h.nxr-1)
// 	  	&& !(pgc.cxr[ixr] > gc.cxr[ir0] && pgc.cxr[ixr] <= gc.cxr[ir0+1]))
// 	  	ir0++;
// 	  	//if(!myrank)
// 	  	//	printf("%f < %f <= %f\n",gc.gxr[ir0], pgc.gxr[ixr], gc.gxr[ir0+1]  );
	  
// 	  for(ixt = 0; ixt < ph.nxt; ixt++) {
// 	  	it0 = 0;
// 	  	interp_case &= 3; // erase all bits above the first 2
// 	  	if(pgc.cxt[ixt] <= gc.cxt[0]) { // min edge of theta
// 	  		interp_case += 4;
// 	  	} else if(pgc.cxt[ixt] >= gc.cxt[h.nxt-1]) { // max edge of theta
// 			interp_case += 8;
// 	  	} // else, theta inside domain (already 0-ed out)

// 	  	// Find left most neutral index next to the current plasma index position for theta
// 	  	if(pgc.cxt[ixt] > gc.cxt[0])
// 	  	while((it0 < h.nxt-1)
// 	  		&& !(pgc.cxt[ixt] > gc.cxt[it0] && pgc.cxt[ixt] <= gc.cxt[it0+1] ))
// 	  		it0++;
// 	  	//if(!myrank && ixr == 0)
// 	  	//	printf("%f < %f <= %f\n", gc.gxt[it0], pgc.gxt[ixt], gc.gxt[it0+1]);
	  	
// 	    for(ixp = 0; ixp < ph.nxp; ixp++) {
// 	    	ip0 = 0;
// 	      // Find left most neutral index next to the current plasma index position for phi
// 	  	  if(pgc.cxp[ixr] > gc.cxp[0])
// 	  	  while((ip0 < h.nxp-1)
// 	  	  	&& !(pgc.cxp[ixp] > gc.cxp[ip0] && pgc.cxp[ixp] <= gc.cxp[ip0+1]))
// 	  	  	ip0++;
// 		  //if(!myrank && ixr == 0 && ixt == 0)
// 	  	  //printf("%f < %f <= %f\n", gc.gxp[ip0], pgc.gxp[ixp], gc.gxp[ip0+1]);
	  	
// 	      idx = arrayIdx3(ixr, ixt, ixp, pnx);

// 	      // if(!myrank) {
// 	      // 	printf("Proton idx: %f %f %f\n Neutral idx: %f %f %f\n",
// 	      // 	  pgc.cxr[ixr], pgc.cxt[ixt], pgc.cxp[ixp], gc.cxr[ir0], gc.cxt[it0], gc.cxp[ip0]);
// 	      // }

// 	      // *******
// 	      // NOTE: instead of having the crazy below logic, just put everything in
// 	      //       the LINT functions and let them figure it out, which in the end
// 	      //       will get rid of the interp_case variable
// 	      // ******
// 	      data_2p[idx]
// 	        = LINT3D_Spherical(data, h, gc, pgc.cxr[ixr], pgc.cxt[ixt], pgc.cxp[ixp]);

// 	      // if((interp_case & 1) && !(interp_case & 2)) { // at origin

// 	      // } else {
// 	      // 	if((interp_case & 4) || (interp_case & 8)) { // theta at inner or outer edge
// 	      // 	  if(interp_case & 2) { // rho at inner or outer edge

// 	      // 	  } else { // rho is in domain

// 	      // 	  }
// 	      // 	} else { // theta is in domain
// 	      // 	  if(interp_case & 2) { // rho at inner or outer edge

// 	      // 	  } else { // everything is in domain

// 	      // 	  }
// 	      // 	}
// 	      // }



// 	    } // END PHI
// 	  } // END THETA
// 	} // END RHO
// }

// void LINT3D_Spherical(double *data, vel_dist_header &h, grid_cell_sph &gc,
// 	double rho, double theta, double phi) 
// {


// if(gc.gxr[0] <= EPS && pgc.cxr[ixr] <= gc.cxr[0]) { // origin, inner boundary
// 	  	interp_case += 1;
// 	  } else if(gc.gxr[0] > EPS && pgc.cxr[ixr] <= gc.cxr[0]) { // no origin, inner boundary
// 	  	interp_case += 2;
// 	  } else if(pgc.cxr[ixr] >= gc.cxr[h.nxr-1]) { // outer boundary
// 	  	interp_case += 3;
// 	  } // else, rho inside domain (already 0)


// 	// interp rho
// 	if(gc.gxr[0] <= EPS && pgc.cxr[ixr] <= gc.cxr[0]) { // at origin, inner boundary
// 		data_2p[idx] =
// 		  LINTr(data, h, gc,
// 		  LINTu_iyz(data, h, gc, 0, M_PI-theta, (phi+M_PI >= 2.0*M_PI ? phi-M_PI : phi+M_PI)),
// 		  LINTu_iyz(data, h, gc, 0, theta, phi));
// 	} else if(gc.gxr[0] > EPS && pgc.cxr[ixr] <= gc.cxr[0]) {// no origin, inner boundary
// 		data_2p[idx] =

// 	} else if(pgc.cxr[ixr] >= gc.cxr[h.nxr-1]) { // outer boundary

// 	} else { // rho in domain

// 	}
// }

//*** MAIN CODE ***///
int main() {
	timer main_timer, vtk_timer, plasma_timer, source_timer;
	//const char vel_dist_fname[] = "vel_dist_circ2015_4muG_HP120A.grid";
	//const char vel_dist_fname[] = "./circ2015_3muG_HP120B/vel_dist_circ2015_3muG_HP120B_FULL.grid";
	const char vel_dist_fname[] = "./circ2015_3muG_HP120B/vel_dist_circ2015_3muG_HP120B_comp0123.grid";
	//const char plasma_prop_fname[] = "plasma_circ2015_4muG_HP120A.plasma";
	const char plasma_prop_fname[] = "./circ2015_3muG_HP120B/plasma_circ2015_3muG_HP120B.plasma";
	vel_dist_header h, ph;
	// Open bin file and init header
	// Read in plasma file , ios::ate sets file cursor to EOF
	ifstream vfile(vel_dist_fname, ios::in|ios::binary|ios::ate);
	// check if file is opened correctly, exit if not
	cout << endl << vel_dist_fname << endl;
	if(vfile.good()) {
		cout << " File open successfully\n";
	} else {
		cout << " File open failed\n";
		return -1;
	}
	streampos file_size = vfile.tellg(); // get size of file in bytes
	cout << "Size of file = " << file_size/1024.0/1024.0 << " MB\n\n";
	vfile.seekg(0, ios::beg); // start at beginning of file
	
	int    nx[3] = {0, 0, 0}; // {nxr, nxt, nxp}
	int    nv[3] = {0, 0, 0}; // to read the header, these vals don't matter yet
	int    npx[3] = {0, 0, 0};
	int    tidata[3];
	double tddata[3];
	 
	// Read header from binary file
	ff_bskipRead(vfile, 1, 1);
	ff_bdataRead(vfile, 1, 3, 2, tidata, 0, 0, 0, nx);
	h.nxr = tidata[0]; nx[0] = h.nxr;
	h.nxt = tidata[1]; nx[1] = h.nxt;
	h.nxp = tidata[2]; nx[2] = h.nxp;
	ff_bskipRead(vfile, 1, 2);
	ff_bdataRead(vfile, 3, 3, 2, tddata, 0, 0, 0, nx);
	h.xalpha = tddata[0];
	h.xmin   = tddata[1] / 1.5E11;//AU_to_M;
	h.xmax   = tddata[2] / 1.5E11;//AU_to_M;
	ff_bskipRead(vfile, 1, 2);
	ff_bdataRead(vfile, 1, 3, 2, tidata, 0, 0, 0, nv);
	h.nvr = tidata[0]; nv[0] = h.nvr;
	h.nvt = tidata[1]; nv[1] = h.nvt;
	h.nvp = tidata[2]; nv[2] = h.nvp;
	ff_bskipRead(vfile, 1, 2);
	ff_bdataRead(vfile, 3, 3, 2, tddata, 0, 0, 0, nv);
	h.valpha = tddata[0];
	h.vmin   = tddata[1]; // cm/s
	h.vmax   = tddata[2]; // cm/s
	ff_bskipRead(vfile, 1, 2);
	ff_bskipRead(vfile, 3, 3); // LISM_vx, vy, vz
	ff_bskipRead(vfile, 1, 1);
	cout << "Neutral distribution header:\n";
	display_header(h);
	
	// Create position and velocity grid for neutrals
	grid_cell_sph gc;
	init_grid_cell_sph(gc, h);
	//display_grid_cell_sph(gc, h);
	
	// init header for 3D plasma dist
	ph.nxr    = 380; npx[0] = ph.nxr;
	ph.nxt    = 120; npx[1] = ph.nxt;
	ph.nxp    = 118; npx[2] = ph.nxp;
	ph.xalpha = 3.0;
	ph.xmin   = 10.0; // in AU
	ph.xmax   = 1000.0; // in AU
	// set velocity space to nil
	ph.nvr = ph.nvt = ph.nvp = 0;
	ph.valpha = 0.0;
	ph.vmin = ph.vmax = 0.0;
	
	cout << "\nPlasma distribution header:\n";
	display_header(ph);
	
	// Create position and velocity grid for plasma
	grid_cell_sph pgc;
	init_grid_cell_sph(pgc, ph);
	//display_grid_cell_sph(pgc, ph);
	
	// Read plasma data from binary file
	double *VPTR = NULL;
	int pxcells = ph.nxr * ph.nxt * ph.nxp;
	int ipxr, ipxt, ipxp;
	double *pdens, *ptemp, *pvelx, *pvely, *pvelz;
	pdens  = new double[pxcells]; memset(pdens, 0, pxcells*sizeof(double));
	pvelx  = new double[pxcells]; memset(pvelx, 0, pxcells*sizeof(double));
	pvely  = new double[pxcells]; memset(pvely, 0, pxcells*sizeof(double));
	pvelz  = new double[pxcells]; memset(pvelz, 0, pxcells*sizeof(double));
	ptemp  = new double[pxcells]; memset(ptemp, 0, pxcells*sizeof(double));
	double pbulkv[3], curvel[3];
	
	start_timer(plasma_timer);
	
	if(READ_PLASMA_FILE) {
		cout << "READ_PLASMA_FILE\n\n";
		read_plasma_binary(plasma_prop_fname, ph,
			pdens, pvelx, pvely, pvelz, ptemp);
	}
	current_time(plasma_timer);
	
	// Calculate moments
	int xcells = h.nxr * h.nxt * h.nxp;
	int vcells = h.nvr * h.nvt * h.nvp;
	int ixr, ixt, ixp, ivr, ivt, ivp, idx, idv;
	double dvr, dvt, dvp, vx, vy, vz, d3v, wd3v, dv, vxy;
	double dxr, dxt, dxp, d3x, dx, tdens;
	double *dens, *velx, *vely, *velz, *vxvar, *vyvar, *vzvar;
	dens  = new double[xcells]; memset(dens, 0, xcells*sizeof(double));
	velx  = new double[xcells]; memset(velx, 0, xcells*sizeof(double));
	vely  = new double[xcells]; memset(vely, 0, xcells*sizeof(double));
	velz  = new double[xcells]; memset(velz, 0, xcells*sizeof(double));
	vxvar = new double[xcells]; memset(vxvar, 0, xcells*sizeof(double));
	vyvar = new double[xcells]; memset(vyvar, 0, xcells*sizeof(double));
	vzvar = new double[xcells]; memset(vzvar, 0, xcells*sizeof(double));
	
	double vtemp[3], rmean[3];
	int vcell_count, i;
	
	// init vars for flux moments (charge exchange, etc), follows plasma grid
	double *fdens, *fvelx, *fvely, *fvelz, *fenrg;
	fdens  = new double[pxcells]; memset(fdens, 0, pxcells*sizeof(double));
	fvelx  = new double[pxcells]; memset(fvelx, 0, pxcells*sizeof(double));
	fvely  = new double[pxcells]; memset(fvely, 0, pxcells*sizeof(double));
	fvelz  = new double[pxcells]; memset(fvelz, 0, pxcells*sizeof(double));
	fenrg  = new double[pxcells]; memset(fenrg, 0, pxcells*sizeof(double));
	
	float *dist_chunk;
	dist_chunk = new float[h.nvp];
	dvt = gc.gvt[1] - gc.gvt[0];
	dvp = gc.gvp[1] - gc.gvp[0];
	dxt = gc.gxt[1] - gc.gxt[0];
	dxp = gc.gxp[1] - gc.gxp[0];

	// get the vel space of a few positions
	double *IHS_vspace, *OHS_vspace, *HP_vspace;
	IHS_vspace = new double[vcells]; memset(IHS_vspace, 0, vcells*sizeof(double));
	OHS_vspace = new double[vcells]; memset(OHS_vspace, 0, vcells*sizeof(double));
	HP_vspace  = new double[vcells]; memset(HP_vspace , 0, vcells*sizeof(double));
	// this line must be after the new operator is called, NOT BEFORE
	double *vspace_array[3] = {IHS_vspace, HP_vspace, OHS_vspace}; 

	double vspace_locs[3] = {93.5, 107, 122};
	bool store_vspace = 0;
	int vspace_count = 0;

	start_timer(main_timer);
	// Loop over all positions
	
	if(READ_NEUTRAL_FILE_MOMENTS) {
	    cout << "READ_NEUTRAL_FILE_MOMENTS\n\n";
		for(ixr = 0; ixr < h.nxr; ixr++) {
		  for(ixt = 0; ixt < h.nxt; ixt++) {
			for(ixp = 0; ixp < h.nxp; ixp++) {
			  idx = arrayIdx3(ixr, ixt, ixp, nx); // index into neutral arrays
			  vcell_count = 0;

			  // // check if need to store vspace
			  if(vspace_count < 3){
				  if((gc.gxr[ixr] <= vspace_locs[vspace_count]
				  && gc.gxr[ixr+1] >= vspace_locs[vspace_count])
				  &&(gc.gxt[ixt] <= M_PI/2.0 && gc.gxt[ixt+1] >= M_PI/2.0 )
				  &&(gc.gxp[ixp] <= 0.0 && gc.gxp[ixp+1] >= 0.0 ) ) {
				  	store_vspace = 1;
				  	cout << "starting to store vspace " << vspace_count << '|'
				  	  << gc.gxr[ixr] << ' ' << vspace_locs[vspace_count] << ' ' << gc.gxr[ixr+1]
				  	  << endl;
				  } else {
				  	store_vspace = 0;
				  }
			   }

			  // Loop over velocity space
			  for(ivr = 0; ivr < h.nvr; ivr++) {
				dvr = gc.gvr[ivr+1] - gc.gvr[ivr];
				dv  = dvr * dvp * dvt;
			    for(ivt = 0; ivt < h.nvt; ivt++) {
			      // Read next chunk from bin file
			      ff_bskipRead(vfile, 1, 1);
			      ff_bdataRead(vfile, 2, h.nvp, 2, dist_chunk, 0, 0, 0, nv);
			      ff_bskipRead(vfile, 1, 1);
			      
			      d3v = gc.cvr[ivr] * gc.cvr[ivr] * sin(gc.cvt[ivt]) * dv;
			      for(ivp = 0; ivp < h.nvp; ivp++) { 
			        // if(ixr+ixt+ixp == 0 && ivp == 0)
			        // 	cout << ivr << ' ' << ivt << ' ' << ivp << " | " << d3v << endl;
					if(dist_chunk[ivp] > 0.0) { // avoid div by zero
						wd3v = dist_chunk[ivp] * d3v;
						// store vel space
						if(store_vspace) {
							//cout << vspace_count << '|'; 
							//cout << arrayIdx3(ivr, ivt, ivp, nv) << ' ' << vcells << endl;
							vspace_array[vspace_count][arrayIdx3(ivr, ivt, ivp, nv)] = dist_chunk[ivp];
						}

						// in units of m/s
				  		vz  = gc.cvr[ivr] * cos(gc.cvt[ivt]);
				  		vxy = gc.cvr[ivr] * sin(gc.cvt[ivt]);
				  		vx  = vxy * cos(gc.cvp[ivp]);
				  		vy  = vxy * sin(gc.cvp[ivp]);
						// prestore first bulk speed values
						if((ivp + ivt + ivr)== 0) {
						  velx[idx] = vx;
						  vely[idx] = vy;
						  velz[idx] = vz;
						}
				        // online variance and mean algorithm (West 1979)
				  		vtemp[0] = vx - velx[idx];
				  		vtemp[1] = vy - vely[idx];
				  		vtemp[2] = vz - velz[idx];
				  		
				  		tdens = dens[idx] + wd3v;
				  		
				  		rmean[0] = vtemp[0] * wd3v / tdens;
				  		rmean[1] = vtemp[1] * wd3v / tdens;
				  		rmean[2] = vtemp[2] * wd3v / tdens;
				  		
				  		velx[idx] += rmean[0];
				  		vely[idx] += rmean[1];
				  		velz[idx] += rmean[2];
				  		
				  		vxvar[idx] += rmean[0] * dens[idx] * vtemp[0];
				  		vyvar[idx] += rmean[1] * dens[idx] * vtemp[1];
				  		vzvar[idx] += rmean[2] * dens[idx] * vtemp[2];
				  		
				  		dens[idx] = tdens;			  		
				  		vcell_count++;
			  		}
			  	  }
			  	}
			  } // Finish vel space loop
			  if(store_vspace) {
			  	cout << "Finished storing vspace at " << vspace_locs[vspace_count]
			  	  << "AU\n"; 
			  	vspace_count++;
			  	store_vspace = 0;
			  }

			  // Normalize to density (cancels out d3x)
			  // [velocity] = [km/s]
			  vxvar[idx] *= 1.0E-6* vcell_count / ((vcell_count-1.0)*dens[idx]);
			  vyvar[idx] *= 1.0E-6* vcell_count / ((vcell_count-1.0)*dens[idx]);
			  vzvar[idx] *= 1.0E-6* vcell_count / ((vcell_count-1.0)*dens[idx]);
			  velx[idx]  *= 1.0E-3;
			  vely[idx]  *= 1.0E-3;
			  velz[idx]  *= 1.0E-3;
			  // [density] = [1/cm^3]
			  dens[idx]  *= 1.0E-6;
			}
		  }
		  current_time(main_timer);
		  cout << 100.0 * (ixr + 1.0) / double(h.nxr) << "% done\n";
		  display_remaining_time(main_timer, ixr+1.0, h.nxr, 2);
		}
	} // END NEUTRAL CALCS
	vfile.close();
	
	// Print moments to vtk file
	start_timer(vtk_timer); //"./vtk/vel_dist_circ2015_4muG_HP120A.vtk"
	if(READ_NEUTRAL_FILE_MOMENTS)
		print2vtk(
		  "./vtk/vel_dist_circ2015_3muG_HP120B_comp0123.vtk",
			//"./vtk/vel_dist_circ2015_3muG_HP120B_FULL.vtk",
		  h, gc, dens, velx, vely, velz, vxvar, vyvar, vzvar, VPTR, 127);

		// print vspace
		const char *field_names[] = {"IHS_vspace", "HP_vspace", "OHS_vspace"};
		print2vtk_v2_vspace("./vtk/vel_dist_circ2015_3muG_HP120B_FULL_Vspace_93-5_107_122AU.vtk", h, gc,
			vspace_array, NULL, field_names, 3, 0);

	if(READ_PLASMA_FILE)
	{
		// print2vtk(
		//   //"./vtk/plasma_circ2015_4muG_HP120A.vtk",
		//   "./vtk/plasma_circ2015_3muG_HP120B.vtk",
		//   ph, pgc, pdens, pvelx, pvely, pvelz, ptemp, VPTR, VPTR, VPTR, 31);

		double *scalar_fields[] = {pdens, ptemp};
		double *vector_fields[] = {pvelx, pvely, pvelz};
		const char *field_names[] = {"density", "temperature", "bulk_flow"};

		print2vtk_v2("./vtk/plasma_circ2015_3muG_HP120B.vtk", ph, pgc,
			scalar_fields, vector_fields, field_names, 2, 1);
	}


	current_time(vtk_timer);
	
	
	if(READ_NEUTRAL_FILE_LOSS_TERMS && !READ_PLASMA_FILE)
		cout << " Need to read in the plasma file to calculate source terms!\n\n";
	if(READ_NEUTRAL_FILE_LOSS_TERMS && READ_PLASMA_FILE) {
	    cout << "READ_NEUTRAL_FILE_LOSS_TERMS\n\n";
	// calculate source & loss terms (flux to protons from hydrgoen)
	//  To be able to loop over plasma position, while reading in the
	//  neutral file, we need to figure out which plasma indices
	//  correspond to which neutral indices. We are assuming that the
	//  plasma grid is at a higher resolution in all rho, theta, phi directions.
	
	// Create two arrays, one of the linear plasma index (the value) and a second array
	// of the corresponding linear neutral index (which will be the key).
	// Now sort these arrays using the key, and because the merge sort is stable,
	// if the plasma indices are presorted, their ordering will be preserved for
	// similar key values (which is what we want)
		//ofstream tfile("test_sort.txt");
		int read_count = 0, pivr, pivt, pivp;
		bool integrate_dens = 0;
		int *index_p, *index_n, *temp_B, *temp_Bk;
		float *vel_chunk;
		double *vxn, *vyn, *vzn, *d3vn, *dens_2p;
		double nextr, nextt, nextp, prevr, prevt, prevp, tbeta_wd3v, source[5];
		index_p   = new int[pxcells];
		index_n   = new int[pxcells];
		temp_B    = new int[pxcells];
		temp_Bk   = new int[pxcells];
		vel_chunk = new float[vcells];
		d3vn      = new double[vcells];
		vxn       = new double[vcells];
		vyn       = new double[vcells];
		vzn       = new double[vcells];
		dens_2p   = new double[pxcells];
		
		memset(dens, 0, xcells*sizeof(double));
		memset(vxn , 0, vcells*sizeof(double));
		memset(vyn , 0, vcells*sizeof(double));
		memset(vzn , 0, vcells*sizeof(double));
		
		// interp neutral density onto plasma grid
		//interp_h2p_grid(dens, dens_2p, h, ph, gc, pgc);


		// Calc conversion of plasma idx to neutral idx
		ivr = 1;
		for(ixr = 0; ixr < ph.nxr; ixr++) {
		  // rho: find closest cell in neutral grid to plasma cell
		  prevr = h.xmax;
		  for(ivr = 0; ivr < h.nxr; ivr++) {
		    nextr = abs(pgc.cxr[ixr] - gc.cxr[ivr]);
		    if(nextr < prevr) {
		      pivr = ivr;
		      prevr = nextr;
		    }
		  }
		  ivr = pivr;
		  for(ixt = 0; ixt < ph.nxt; ixt++) {
		    // theta: find closest cell in neutral grid to plasma cell
		    prevt = M_PI;
		    for(ivt = 0; ivt < h.nxt; ivt++) {
		      nextt = abs(pgc.cxt[ixt] - gc.cxt[ivt]);
		      if(nextt < prevt) {
		        pivt = ivt;
		        prevt = nextt;
		      }
		    }
		    ivt = pivt;
		    for(ixp = 0; ixp < ph.nxp; ixp++) {
		      // phi: find closest cell in neutral grid to plasma cell
		      prevp = 2.0*M_PI;
		      for(ivp = 0; ivp < h.nxp; ivp++) {
		        nextp = abs(pgc.cxp[ixp] - gc.cxp[ivp]);
		        if(nextp < prevp) {
		          pivp = ivp;
		          prevp = nextp;
		        }
		      }
		      ivp = pivp;
		      // plasma cell index
		      idx = arrayIdx3(ixr, ixt, ixp, npx);
		      index_p[idx] = idx;
		      // neutral cell index
		      index_n[idx] = arrayIdx3(ivr, ivt, ivp, nx);
		    }
		  }
		}
		/*
		// Calc conversion of plasma idx to neutral idx
		ivr = 1;
		for(ixr = 0; ixr < ph.nxr; ixr++) {
		  // search for the r index of the neutrals just to the left of the
		  // plasma grid center
		  ivr--;
		  // assuming the neutral grid fits inside the plasma grid
		  while(gc.gxr[ivr] < pgc.cxr[ixr]) ivr++; // gc.gxr[ivr] < pgc.cxr[ixr]
		  ivr--;
		  
		  for(ixt = 0; ixt < ph.nxt; ixt++) {
		    ivt = (h.nxt - 1.0) / (ph.nxt - 1.0) * ixt;
		    for(ixp = 0; ixp < ph.nxp; ixp++) {
		      // plasma grid index
		      idx = arrayIdx3(ixr, ixt, ixp, npx);
		      index_p[idx] = idx;
		      // calc neutral grid indices (reusing variables)
		      ivp = (h.nxp - 1.0) / (ph.nxp - 1.0) * ixp;
		      index_n[idx] = arrayIdx3(ivr, ivt, ivp, nx);
		    }
		  }
		}*/
		
		// sort by key and value if keys are similar
		cout << "Starting sorting of plasma indices wrt to neutral indices\n\n";
		//insert_sort_key(index_p, index_n, pxcells);
		mergesort(index_p, index_n, temp_B, temp_Bk, 0, pxcells-1);
		//for(i=0; i < pxcells; i++)
		//	tfile << index_n[i] << ' ' << index_p[i] << ' ' << endl;
		//tfile.close();
		delete[] temp_B, temp_Bk;
	    // Once the arrays are sorted, check all keys for out of bounds keys.
	    // For the out of bounds keys, make the corresponding value -1 as to 
	    // skip over this iteration.
		for(i = 0; i < pxcells; i++) {
			if(index_n[i] < 0 || index_n[i] >= xcells)
				index_n[i] = -1;
		}
		// precalculate the velocities and d3v in neutral vel space
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
	
	    start_timer(source_timer);
		// Read in neutral file
		vfile.open(vel_dist_fname, ios::in|ios::binary);
		// skip header info, this was already used to init the neutral grids
		ff_bskipRead(vfile, 1, 6);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 7);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 2);
		ff_bskipRead(vfile, 3, 3);
		ff_bskipRead(vfile, 1, 1);
		
		// loop over plasma position grid
		for(ixr = 0; ixr < ph.nxr; ixr++) {
		  for(ixt = 0; ixt < ph.nxt; ixt++) {
		    for(ixp = 0; ixp < ph.nxp; ixp++) {
		      idx = arrayIdx3(ixr, ixt, ixp, npx); // index into plasma grid
		      // Is the neutral index valid?, if not skip everything
		      if(index_n[idx] != -1) {
				  if(read_count <= index_n[idx]) {
				  	read_count++;
				  	integrate_dens = 1;
				  	// Read next chunk of velocity space from bin file
				  	for(ivr = 0; ivr < h.nvr; ivr++) {
					  for(ivt = 0; ivt < h.nvt; ivt++) {
					    // Read next chunk from bin file
					    ff_bskipRead(vfile, 1, 1);
					    ff_bdataRead(vfile, 2, h.nvp, 2, vel_chunk, ivr, ivt, 0, nv);
					    ff_bskipRead(vfile, 1, 1);
					  }
					}
				  }
				  // bulk vel in plasma at idxp position, convert from km/s to m/s
				  pbulkv[0] = 1000.0 * pvelx[index_p[idx]];
				  pbulkv[1] = 1000.0 * pvely[index_p[idx]];
				  pbulkv[2] = 1000.0 * pvelz[index_p[idx]];
				  
				  // loop over neutral velocity grid
				  for(ivr = 0; ivr < h.nvr; ivr++) {
				    for(ivt = 0; ivt < h.nvt; ivt++) {
				      //#pragma omp parallel for reduction(+:fdens[index_p[idx]], fvelx[index_p[idx]] \
				      //fvely[index_p[idx]], fvelz[index_p[idx]], fenrg[index_p[idx]]) \
				     // private(curvel[0], curvel[1], curvel[2], wd3v, idv)
				      for(ivp = 0; ivp < h.nvp; ivp++) {
				        idv = arrayIdx3(ivr, ivt, ivp, nv);
						
						if(vel_chunk[idv] > 0.0) { // avoid div by zero
						  wd3v = vel_chunk[idv] * d3vn[idv];
						  // current vel of vel space cell, in m/s
						  curvel[0] = vxn[idv];
					  	  curvel[1] = vyn[idv];
					  	  curvel[2] = vzn[idv];
					  	  
					  	  // only integrate the density for the first time going over
					  	  // the vel space
					  	  if(integrate_dens){
					  	  //#pragma omp critical
					  	   // {
					  	  	  dens[index_n[idx]] += wd3v;
					  	  //	}
					  	  }
					  	 // tbeta_wd3v = beta_sqrt_kappa_IHS(pdens[index_p[idx]],
						 // 	ptemp[index_p[idx]], &charge_ex, pbulkv ,MP ,curvel,
						  //	1.63, 8000.0) * wd3v; // kappa = 1.63, temp_LISM = 8000 K
					  	  //tbeta_wd3v = beta_sqrt(pdens[index_p[idx]],
						  //	ptemp[index_p[idx]], &charge_ex, pbulkv ,MP ,curvel) * wd3v;
						  
						  //source_maxwellian_erf(pdens[index_p[idx]], ptemp[index_p[idx]],
						  //&charge_ex, pbulkv ,MP ,curvel, source[0], source+1, source[4]);
						  
						  source_maxwellian_sqrt(pdens[index_p[idx]], ptemp[index_p[idx]],
						  &charge_ex, pbulkv ,MP ,curvel, source[0], source+1, source[4]);

						  fdens[index_p[idx]] += source[0] * wd3v;//tbeta_wd3v;
						  // MP is factored out here and done in normalizing
						  fvelx[index_p[idx]] += source[1] * wd3v;//vxn[idv] * tbeta_wd3v;
						  fvely[index_p[idx]] += source[2] * wd3v;//vyn[idv] * tbeta_wd3v;
						  fvelz[index_p[idx]] += source[3] * wd3v;//vzn[idv] * tbeta_wd3v;
						  fenrg[index_p[idx]] += source[4] * wd3v;//(vxn[idv]*vxn[idv] + vyn[idv]*vyn[idv]
						    //+ vzn[idv]*vzn[idv]) * tbeta_wd3v;
						}
				      }
				    }
				  } // END of NEUTRAL VEL GRID LOOP
				  // do not normalize density in order to cancel out units of wd3v,
				  //   which is also not normalized
				  //if(integrate_dens) // normalize the density
				  //	dens[index_n[idx]] *= 1.0E-6;
				  integrate_dens = 0;
				  // Normalize to density (cancels out d3x)
				  
				  // charge exchange rate [fdens] = [1/s / m^3]
				  fdens[index_p[idx]] *= 1.0E-20; // not sure about this number
				  // momentum exchange rate [fvel] = [kg*m/s /s / m^3]
				  fvelx[index_p[idx]] *= 1.0E-18;
				  fvely[index_p[idx]] *= 1.0E-18;
				  fvelz[index_p[idx]] *= 1.0E-18;
				  // energy exchange rate [fenrg] = [J / s / m^3]
				  fenrg[index_p[idx]] *= 1.0E-16;
				  
		      } // End of valid neutral index check
		    }
		  }
		  current_time(source_timer);
		  cout << 100.0 * (ixr + 1.0) / double(ph.nxr) << "% done\n";
		  display_remaining_time(source_timer, ixr+1.0, ph.nxr, 3);
		}
		vfile.close();
		print2vtk(
		  //"./vtk/plasma_circ2015_4muG_HP120A.vtk",
		  "./vtk/plasma_circ2015_3muG_HP120B_comp0123_source_terms_maxwellian_total_sqrt.vtk",
		  ph, pgc, fdens, fvelx, fvely, fvelz, VPTR, VPTR, VPTR, fenrg, 143);
		
		delete[] index_p, index_n, vel_chunk, vxn, vyn, vzn, d3vn, dens_2p;
	} // End Source Term Calculations
	
	// read source terms from Jacob's source file
	if(READ_OUTPUT_SOURCE) {
		cout << "READ_OUTPUT_SOURCE\n\n";
		ifstream source_file("./circ2015_3muG_HP120B/source_circ2015_3muG_HP120B.source");
		//npx[0]
		double *mass, *px, *py, *pz, *ee, *nchex, dud;
		mass  = new double[pxcells];
		px    = new double[pxcells];
		py    = new double[pxcells];
		pz    = new double[pxcells];
		ee    = new double[pxcells];
		nchex = new double[pxcells];
		
		// loop over plasma position grid. This assumes the source file
		// has the same dimenstions in rho, theta, and phi
		for(ixp = 0; ixp < ph.nxp; ixp++) {
		  for(ixt = 0; ixt < ph.nxt; ixt++) {
		    for(ixr = 0; ixr < ph.nxr; ixr++) {
		      // this index reorders the data so it's phi as the fastest changing dim, not rho
		      idx = arrayIdx3(ixr, ixt, ixp, npx); // index into plasma grid
		      
		      source_file >> dud >> dud >> dud; // read rho, theta, phi vals
		      source_file >> mass[idx]; // read mass, ie density
		      source_file >> px[idx] >> py[idx] >> pz[idx]; // momentum source terms
		      source_file >> ee[idx]; // energy source term
		      source_file >> nchex[idx]; // number of charge exchanges
		      //cout << nchex[idx] << endl;
		    }
		  }
		}
		
		double *scalar_fields[] = {mass, nchex, ee};
		double *vector_fields[] = {px, py, pz};
		const char *field_names[] =
			{"density", "charge_exchange_count", "energy_source", "momentum_source"};

		print2vtk_v2("./vtk/plasma_circ2015_3muG_HP120B_source_terms_jacob.vtk", ph, pgc,
			scalar_fields, vector_fields, field_names, 3, 1);

		// print2vtk(
		//   "./vtk/plasma_circ2015_3muG_HP120B_source_terms_jacob.vtk",
		//   ph, pgc, nchex, px, py, pz, VPTR, VPTR, VPTR, ee, 143);
		
		delete[] mass, px, py, pz, ee, nchex;
	}
	
	
	current_time(main_timer);
	cout << "\n Main timer:\n";
	display_elapse_time(main_timer, 3);
	cout << " VTK printer timer:\n";
	display_elapse_time(vtk_timer);
	cout << " Plasma reader timer:\n";
	display_elapse_time(plasma_timer);
	cout << " Source Calc timer:\n";
	display_elapse_time(source_timer, 3);
	
	
	delete[] dist_chunk, dens, velx, vely, velz;
	delete[] IHS_vspace, OHS_vspace, HP_vspace;
	delete[] vxvar, vyvar, vzvar;
	delete[] pdens, ptemp, pvelx, pvely, pvelz;
	delete[] fdens, fvelx, fvely, fvelz, fenrg;
	delete_grid_cell_sph(gc);
	delete_grid_cell_sph(pgc);
	return 0;
}

