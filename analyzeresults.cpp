/************************************************************************
analyzeresults for NetConductV1
TWS 2023
*************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "nrutil.h"

void histogram(float *var, float *weight, int n, const char filename[]);

void analyzeresults()
{
	extern int nnod, nseg, varyviscosity, *ista, *iend;
	extern int *nodname, *nodtyp, *segname, *segtyp, **nodseg;
	extern float constvisc, pi1;
	extern float *histogramdisplay, *histogramweight, *diam, *seg_pot;
	extern float *length_weight, *lseg, *condq, **cnode;
	extern float x_stim, y_stim, z_stim;
	extern double *nod_pot;

	int iseg, inod, i;
	int min_pot_nod = 0, max_pot_nod = 0;
	float mean_nod_pot, nod_pot_deviation, totallength;
	float max_pot, min_pot, source_dist;
	float *cseg;
	FILE *ofp;

	cseg = vector(1, 3);
	totallength = 0.;
	mean_nod_pot = 0.;
	min_pot = 1.e6;
	max_pot = -1.e6;
	for (inod = 1; inod <= nnod; inod++) {
		//node weighting factors: 1/2 the sum of the lengths of the segments connected to the node
		length_weight[inod] = 0;
		for (iseg = 1; iseg <= nodtyp[inod]; iseg++) length_weight[inod] += 0.5*lseg[nodseg[iseg][inod]];
		totallength += length_weight[inod];
		histogramdisplay[inod] = -nod_pot[inod];
		histogramweight[inod] = length_weight[inod];
		// calculate average and max/min potentials
		mean_nod_pot += nod_pot[inod] * length_weight[inod];
		if (nod_pot[inod] < min_pot) {
			min_pot = nod_pot[inod];
			min_pot_nod = inod;
		}
		if (nod_pot[inod] > max_pot) {
			max_pot = nod_pot[inod];
			max_pot_nod = inod;
		}
	}
	mean_nod_pot = mean_nod_pot / totallength;
	nod_pot_deviation = 0.;
	for (inod = 1; inod <= nnod; inod++) nod_pot_deviation += length_weight[inod] * SQR(nod_pot[inod] - mean_nod_pot);
	nod_pot_deviation = sqrt(nod_pot_deviation / totallength);
	histogram(histogramdisplay, histogramweight, nnod, "Current/histogram-potentials.out");

	ofp = fopen("Current/Results.out", "w");
	fprintf(ofp, "All quantities weighted by segment length\n");
	fprintf(ofp, "Nodal potential mean +- s.d.: %6f +- %6f\n", mean_nod_pot, nod_pot_deviation);
	fprintf(ofp, "Maximum potential: %6f at node %4i\n", max_pot, nodname[max_pot_nod]);
	fprintf(ofp, "Minimum potential: %6f at node %4i\n", min_pot, nodname[min_pot_nod]);
	fprintf(ofp, "  seg. name length  conductance potent. coordinates x y z source distance\n");
	for (iseg = 1; iseg <= nseg; iseg++) {
		for (i = 1; i <= 3; i++) cseg[i] = (cnode[i][ista[iseg]] + cnode[i][iend[iseg]]) / 2.;
		source_dist = sqrt(SQR(x_stim - cseg[1]) + SQR(y_stim - cseg[2]) + SQR(z_stim - cseg[3]));
		fprintf(ofp, "%4i %4i %8.3f %8.3f %8.3f %6.1f %6.1f %6.1f %6.1f\n ",
			iseg, segname[iseg], lseg[iseg], condq[iseg], (nod_pot[ista[iseg]]+nod_pot[iend[iseg]])/2.,
			cseg[1], cseg[2], cseg[3], source_dist);
	}
	fclose(ofp);

	ofp = fopen("Current/Potentials.out", "w");
	for (inod = 1; inod <= nnod; inod++) fprintf(ofp, "%4i %8.3f \n ", inod, nod_pot[inod]);
	fclose(ofp);
}