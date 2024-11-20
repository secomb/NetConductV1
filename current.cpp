/************************************************************************
current - for NetConductV1
Electrophysiology of conducted repsonses, including Kir channels
TWS 2023
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void solveq(double *nod_pot);

void current()
{
	extern int nseg, nnod, nitmax1;
	extern int *segtyp, *nodtyp, **nodseg, *ista, *iend;
	extern float F, g_bg, g_d, g_gj, g_Kir0, g_KNV, K_e0, K_i, k_Kir, v_Kir, l_EC, d_EC, R, T, v_bg;
	extern float pot_tol;
	extern float *diam, *q, *condq, *lseg, *seg_pot;
	extern float *n_EC_seg, *n_EC_nod, *n_EC_gj, *E_K, *K_e, *qjnodmem0, *qjnodmem1;
	extern double *nod_pot, *nod_pot_old;

	int i, iseg, inod, niter, errnodpot = 0;
	float maxpoterr = 0., potchange, g_Kir, flow_pot;

	for (iseg = 1; iseg <= nseg; iseg++) if (segtyp[iseg] == 4 || segtyp[iseg] == 5) {
		n_EC_seg[iseg] = diam[iseg] / d_EC * lseg[iseg] / l_EC;		//Number of endothelial cells per segment
		n_EC_gj[iseg] = diam[iseg] / d_EC * l_EC / lseg[iseg];		//Gap junction conductance factor due to endothelial cells
	}
	for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] > 0) {
		E_K[inod] = R * T / F * log(K_e[inod] / K_i) * 1.e3;		//Node K+ Nernst potential in mV
		n_EC_nod[inod] = 0.;
		for (i = 1; i <= nodtyp[inod]; i++) {		//number of endothelial cells associated with node
			iseg = nodseg[i][inod];
			if (segtyp[iseg] == 4 || segtyp[iseg] == 5) n_EC_nod[inod] += n_EC_seg[iseg] / 2.;
		}
		nod_pot[inod] = -30.5;		//initialize potentials
	}

	for(niter=1; niter<=nitmax1; niter++){
		for (inod = 1; inod <= nnod; inod++) if (nodtyp[inod] > 0) {
			nod_pot_old[inod] = nod_pot[inod];		//save old value for convergence testing
			g_Kir = g_Kir0 * sqrt(K_e[inod]) / (1 + exp((nod_pot[inod] - (E_K[inod] + v_Kir)) / k_Kir)); //Kir channel conductance per EC
			qjnodmem0[inod] = -n_EC_nod[inod] * (g_bg * v_bg + (g_Kir + g_KNV) * E_K[inod]);	//membrane current per node, const term
			qjnodmem1[inod] = n_EC_nod[inod] * (g_bg + g_Kir + g_KNV);	//membrane current per node, linear term
		}
		for(iseg = 1; iseg <= nseg; iseg++) if(segtyp[iseg] == 4 || segtyp[iseg] == 5){
			//if flow_pot > 0, current same direction as flow, no extra resistance to conducted response
			//if flow_pot < 0, current opposite direction as flow, extra resistance to conducted response
			flow_pot = (nod_pot[ista[iseg]] - nod_pot[iend[iseg]]) * q[iseg];
			//Segment gap junction conductance
			if (flow_pot >= 0.) condq[iseg] = n_EC_gj[iseg] * g_gj; 
			//else condq[iseg] = n_EC_gj[iseg] * g_d * g_gj / (g_d + g_gj); //Effect of diode conductance
			else condq[iseg] = n_EC_gj[iseg] * g_d; //Effect of diode conductance
		}

		solveq(nod_pot);

//compare nod_pot with previous values
		maxpoterr = 0.;
		errnodpot = 0;
		for(inod=1; inod<=nnod; inod++) if (nodtyp[inod] > 0) {
			potchange = nod_pot[inod] - nod_pot_old[inod];
			if(fabs(potchange) >= maxpoterr){
				maxpoterr = fabs(potchange);
				errnodpot = inod;
			}
		}
		if(maxpoterr < pot_tol) goto converged;
	}
	printf("*** Warning: Nonlinear iteration not converged\n");
	printf("*** current error = %f at node %i\n", maxpoterr, errnodpot);
	converged:;
	printf("current: %i iterations\n",niter);
	for (iseg = 1; iseg <= nseg; iseg++) seg_pot[iseg] = (nod_pot[ista[iseg]] + nod_pot[iend[iseg]]) / 2.;
}
