/************************************************************************
NetConductV1
Program to compute conducted responses in microvascular networks
Reads Network.dat file, computes currents and writes NetworkNew.dat.
Note that only segments with segtyp = 4 or 5 in network.dat are included.
TWS 2023, based on work of Sara Djurich
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include <iostream>

using namespace std;

#if defined(__linux__)
	// Requires c++17 support, should be included in all current linux releases
	#include <experimental/filesystem> 
	namespace fs = std::experimental::filesystem::v1;
#elif defined(__APPLE__)
	// Requires removal of the -lstdc++fs flag from makefile
	#include <filesystem>
	namespace fs = std::filesystem;
#elif defined(_WIN32)    //Windows version
	#include <Windows.h>
#endif

void input();
void analyzenet();
void setuparrays1(int nseg, int nnod);
void current();
void cmgui(float *segvar);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void picturenetwork(float *nodvar, float *segvar, const char fname[]);
void analyzeresults();

int max=100,nseg,nnod,nnodbc,niter,nnodfl,nsegfl,mxx,myy,mzz,nodsegm;
int nitmax1,nitmax, type_stim;
int *segtyp,*segname,*bcnodname,*bcnod,*bctyp,*nodname,*nodtyp;
int *nodrank,*nodout,*nk,*ista,*iend;
int **nodnod,**nodseg,**segnodname;

float F, g_bg, g_d, g_gj, g_Kir0, g_KNV, K_e0, K_i, k_Kir, v_Kir, l_EC, d_EC, R, T, v_bg;
float alx, aly, alz, lb, maxl;
float tol, pot_tol, omega, optw, optlam;
float K_stim, r_stim, x_stim, y_stim, z_stim;
float *diam, *q, *bcpotqj, *lseg, *condq, *condqsum, *seg_pot;
float *segvar, *nodvar, *xsl0, *xsl1, *xsl2;
float *length_weight, *histogramdisplay, *histogramweight;
float **cnode;
float *n_EC_seg, *n_EC_nod, *n_EC_gj, *E_K, *K_e, *qjnodmem0, *qjnodmem1;
double *nod_pot, *nod_pot_old;
FILE *ofp;

int main(int argc, char *argv[])
{
	int iseg, inod, nnod_stim;
	float r2, r_stim2;

	//Create a Current subdirectory if needed. Copy data files to it.
#if defined(__unix__)
	if (!fs::exists("Current")) fs::create_directory("Current");
	fs::copy_file("ContourParams.dat", fs::path("Current/ContourParams.dat"), fs::copy_options::overwrite_existing);
	fs::copy_file("Network.dat", fs::path("Current/Network.dat"), fs::copy_options::overwrite_existing);
	fs::copy_file("ConductParams.dat", fs::path("Current/ConductParams.dat"), fs::copy_options::overwrite_existing);
#elif defined(_WIN32)
	BOOL NoOverwrite = FALSE;
	DWORD ftyp = GetFileAttributesA("Current\\");
	if (ftyp != FILE_ATTRIBUTE_DIRECTORY) system("mkdir Current");
	CopyFile("ContourParams.dat", "Current\\ContourParams.dat", NoOverwrite);
	CopyFile("Network.dat", "Current\\Network.dat", NoOverwrite);
	CopyFile("ConductParams.dat", "Current\\ConductParams.dat", NoOverwrite);
#endif

	input();	//read data files

	setuparrays1(nseg, nnod);

	analyzenet();

	for (iseg = 1; iseg <= nseg; iseg++) segvar[iseg] = iseg;
	for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] > 0) nodvar[inod] = inod;
	picturenetwork(nodvar, segvar, "Current/nod_seg.ps");	//2D projection of network

	if (type_stim == 0) {		//45-segment
		for (inod = 1; inod <= nnod; inod++) K_e[inod] = K_e0;
		K_e[18] = K_stim;
		K_e[19] = K_stim;
		K_e[20] = K_stim;
	}
	else if (type_stim == 1) {	//strict sphere
		nnod_stim = 0;
		for (inod = 1; inod <= nnod; inod++) {	//set extracellular potassium
			r2 = SQR(cnode[1][inod] - x_stim) + SQR(cnode[2][inod] - y_stim) + SQR(cnode[3][inod] - z_stim);
			if (r2 <= SQR(r_stim)) {
				K_e[inod] = K_stim;
				nnod_stim++;
			}
			else K_e[inod] = K_e0;
		}
		printf("nnod_stim = %i\n", nnod_stim);
	}
	else if (type_stim == 2) {	//Gaussian sphere
		r_stim2 = SQR(r_stim);
		for (inod = 1; inod <= nnod; inod++) {	//set extracellular potassium
			r2 = SQR(cnode[1][inod] - x_stim) + SQR(cnode[2][inod] - y_stim) + SQR(cnode[3][inod] - z_stim);
			K_e[inod] = K_e0 + (K_stim - K_e0) * exp(-r2 / r_stim2);
		}
	}
	else printf("*** Error: Invalid stimulus type ***\n");

	current();

	analyzeresults();	//statistics and histograms of resulting flows

	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) segvar[iseg] = seg_pot[iseg];
	for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] > 0) nodvar[inod] = nod_pot[inod];

	cmgui(segvar);		//files for 3D image of network

	picturenetwork(nodvar, segvar, "Current/Seg_pot.ps");	//2D projection of network

	cout << endl << "Press enter to exit program" << endl;
	std::cin.get();
		
	return 0;
}