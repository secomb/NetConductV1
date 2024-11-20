/************************************************************************
setuparrays1 - for NetConductV1
Set up arrays with dimensions nnod and nseg
TWS 2023
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void setuparrays1(int nseg, int nnod)
{
	extern int nsp,nodsegm;
	extern int *nk,*nodrank,*nodout,*nodtyp;
	extern int *ista, *iend;
	extern int **nodnod,**nodseg;
	extern float *lseg, *condq, *condqsum, *segvar, *nodvar, *seg_pot;
	extern float *length_weight, *histogramdisplay, *histogramweight;
	extern float *n_EC_seg, *n_EC_nod, *n_EC_gj, *E_K, *K_e, *qjnodmem0, *qjnodmem1;
	extern double *nod_pot, *nod_pot_old;

	ista = ivector(1,nseg); 
	iend = ivector(1,nseg); 
	nk = ivector(1,nnod); 
	nodout = ivector(1,nnod); 
	nodrank = ivector(1,nnod); 
	nodtyp = ivector(1,nnod);
	E_K = vector(1, nnod); 
	K_e = vector(1, nnod);

	nodnod = imatrix(1,nodsegm,1,nnod); 
	nodseg = imatrix(1,nodsegm,1,nnod); 

	condq = vector(1, nseg);
	condqsum = vector(1, nnod);
	qjnodmem0 = vector(1, nnod);
	qjnodmem1 = vector(1, nnod);
	seg_pot = vector(1, nseg);
	n_EC_seg = vector(1, nseg);
	n_EC_nod = vector(1, nnod);
	n_EC_gj = vector(1, nseg);
	lseg = vector(1,nseg); 
	segvar = vector(1,nseg);
	nodvar = vector(1, nnod);
	length_weight = vector(1, nnod);
	histogramdisplay = vector(1, LMAX(nseg, nnod));
	histogramweight = vector(1, LMAX(nseg, nnod));

	nod_pot = dvector(1, nnod);
	nod_pot_old = dvector(1, nnod);
}
