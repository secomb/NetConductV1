/************************************************************************
solveq - for NetConductV1
iterative solution for potentials in network
Need double precision potentials for networks with many segments
Lengths, diameters in microns, times in s
Potentials in mV
Conductances in nS
Currents in pA
Basic equation at node 0:
Node potentials: p_0, p_i, i = 1, 2,..
Segment conductances: g_i
Membrane current associated with the node is i_0 + g_0*p_0 
sum((p_i - p_0) * g_i) = i_0 + g_0*p_0
p_0 = [sum(p_i * g_i) - i_0]/(sum(g_i) + g_0)
At a boundary node with imposed (inward) current i_in:
(p_1 - p_0) * g_1 + i_in = i_0 + g_0*p_0
p_0 = [p_1 * g_1 + i_in - i_0]/(g_1 + g_0)
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void solveq(double *nod_pot)
{
	extern int nseg,nnod,nnodbc,nitmax,nodsegm;
	extern int *bcnod,*bctyp,*nodtyp,*segtyp;
	extern int **nodseg,**nodnod;
	extern float tol,omega;
	extern float *bcpotqj, *condq, *condqsum, *nodvar, *segvar, *qjnodmem0, *qjnodmem1;

	int iseg, inod, inod1, inodbc, i, j, errnode, niter;
	float maxerr;
	double **wk,pot1,pcondqsum;
	wk = dmatrix(1,nodsegm,1,nnod);

	//potential and current boundary nodes.  Temporarily set nodtyp of potential nodes to -1.
	for (inodbc = 1; inodbc <= nnodbc; inodbc++) {
		inod = bcnod[inodbc];
		if (bctyp[inodbc] == 0) {		//boundary condition on potential
			nod_pot[inod] = bcpotqj[inodbc];
			nodtyp[inod] = -1;
		}
		else wk[1][inod] = bcpotqj[inodbc];	//store current inflow to this boundary node
	}
//set up coefficients
	for(inod=1; inod<=nnod; inod++){
		if(nodtyp[inod] > 0){
			condqsum[inod] = qjnodmem1[inod];
			for(i=1; i<=nodtyp[inod]; i++){
				iseg = nodseg[i][inod];
				condqsum[inod] += condq[iseg];
				if (nodtyp[inod] > 1) wk[i][inod] = condq[iseg];	//store conductances associated with this node
			}
		}
	}
//iterative solution for potentials
	for(niter=1; niter<=nitmax; niter++){
		maxerr = 0.;
		for(inod=1; inod<=nnod; inod++){
			if (nodtyp[inod] == 1) {	//current type boundary condition
				iseg = nodseg[1][inod];
				inod1 = nodnod[1][inod];
				pot1 = omega * ((nod_pot[inod1] * condq[iseg] + wk[1][inod] - qjnodmem0[inod]) / condqsum[inod] - nod_pot[inod]);
			}
			if(nodtyp[inod] >= 2){
				pcondqsum = 0.;
				for(i=1; i<=nodtyp[inod]; i++) pcondqsum += wk[i][inod] * nod_pot[nodnod[i][inod]];
  				pot1 = omega * ((pcondqsum - qjnodmem0[inod]) / condqsum[inod] - nod_pot[inod]);
			}
 			if(nodtyp[inod] >= 1){
				nod_pot[inod] += pot1;
				if(fabs(pot1) >= maxerr){
					maxerr = fabs(pot1);
					errnode = inod;
				}
			}
		}
		if(maxerr < tol) goto converged;
	}
	printf("*** Warning: linear iteration not converged, maxerr = %g at node %i\n", maxerr,errnode);
	converged:;	
	for(inod=1; inod<=nnod; inod++)	if(nodtyp[inod] == -1) nodtyp[inod] = 1;//reset nodtyp of potential nodes
	free_dmatrix(wk,1,nodsegm,1,nnod);
}
