/************************************************************************
input - for NetConductV1
TWS 2023
************************************************************************/
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

void input()
{
	extern int mxx,myy,mzz,nseg,nnod,nnodbc,nodsegm;
	extern int nitmax1,nitmax, type_stim;
	extern int *segname,*segtyp,*bcnodname,*bcnod,*bctyp,*nodname;
	extern int **segnodname; 

	extern float F, g_bg, g_d, g_gj, g_Kir0, g_KNV, K_e0, K_i, k_Kir, v_Kir, l_EC, d_EC, R, T, v_bg;
	extern float alx, aly, alz, lb, maxl;
	extern float tol,pot_tol,omega,optw,optlam;
	extern float K_stim, r_stim, x_stim, y_stim, z_stim;
	extern float *diam,*lseg,*q,*bcpotqj;
	extern float *xsl0, *xsl1, *xsl2;
	extern float **cnode;

	int i,inodbc,iseg,max=200;
	FILE *ifp;
	char bb[200];

	ifp = fopen("ConductParams.dat", "r");
	fscanf(ifp, "%f%*[^\n]", &F);		//Faraday's constant, C/mol
	fscanf(ifp, "%f%*[^\n]", &g_bg);	//Conductance of background current, nS
	fscanf(ifp, "%f%*[^\n]", &g_d);		//Diode conductance, random value
	fscanf(ifp, "%f%*[^\n]", &g_gj);	//Gap junction conductance, nS
	fscanf(ifp, "%f%*[^\n]", &g_Kir0);	//Kir channel maximum conductance in nS / (mM ^ 1 / 2)
	fscanf(ifp, "%f%*[^\n]", &g_KNV);	//Conductance of non - voltage gated K + channels
	fscanf(ifp, "%f%*[^\n]", &K_e0);	//Baseline extracellular potassium in mM
	fscanf(ifp, "%f%*[^\n]", &K_i);		//Intracellular potassium, mM
	fscanf(ifp, "%f%*[^\n]", &k_Kir);	//Kir slope factor, mV
	fscanf(ifp, "%f%*[^\n]", &v_Kir);	//Kir offset factor, mV
	fscanf(ifp, "%f%*[^\n]", &l_EC);	//Length of endothelial cell, um
	fscanf(ifp, "%f%*[^\n]", &d_EC);	//reference diameter for one endothelial cell in circumference, um
	fscanf(ifp, "%f%*[^\n]", &R);		//Gas constant, J*mol^-1 * K^-1
	fscanf(ifp, "%f%*[^\n]", &T);		//Temperature, Kelvin
	fscanf(ifp, "%f%*[^\n]", &v_bg);	//Background Nernst potential, mV
	fscanf(ifp, "%i %f %f%*[^\n]", &nitmax, &tol, &omega);
	fscanf(ifp, "%i %f%*[^\n]", &nitmax1, &pot_tol);
	fscanf(ifp, "%f%*[^\n]", &K_stim);	//potassium concentration in stimulus region
	fscanf(ifp, "%f%*[^\n]", &r_stim);	//radius of stimulus region
	fscanf(ifp, "%f %f %f%*[^\n]", &x_stim, &y_stim, &z_stim);	//coordinates of center of stimulus region
	fscanf(ifp, "%i%*[^\n]", &type_stim);	//type of stimulus. 0=45-segment. 1=strict sphere. 2=Gaussian sphere.
	fclose(ifp);

//network data file
	ifp = fopen("Network.dat", "r");
	fgets(bb,max,ifp);
	printf("%s\n",bb);
	fscanf(ifp, "%f %f %f%*[^\n]", &alx,&aly,&alz);
	fscanf(ifp, "%i %i %i%*[^\n]", &mxx,&myy,&mzz);
	fscanf(ifp, "%f%*[^\n]", &lb);
	fscanf(ifp, "%f%*[^\n]", &maxl);
	fscanf(ifp, "%i%*[^\n]", &nodsegm);
	fscanf(ifp, "%i%*[^\n]", &nseg);
	fgets(bb,max,ifp);
	fgets(bb,max,ifp);
	segname = ivector(1,nseg);
	segtyp = ivector(1,nseg);
	segnodname = imatrix(1,2,1,nseg);
	diam = vector(1,nseg);
	q = vector(1,nseg);
	for(iseg=1; iseg<=nseg; iseg++)	fscanf(ifp, "%i %i %i %i %f %f%*[^\n]", 
		&segname[iseg],&segtyp[iseg],&segnodname[1][iseg],&segnodname[2][iseg],&diam[iseg],&q[iseg]);
//number of nodes in vessel network
	fgets(bb,max,ifp);
	fscanf(ifp,"%i%*[^\n]", &nnod);
	fgets(bb,max,ifp);
	fgets(bb,max,ifp);
//coordinates of nodes
	nodname = ivector(1,nnod);
	cnode = matrix(1,3,1,nnod);
	for(i=1; i<=nnod; i++) fscanf(ifp, "%i %f %f %f%*[^\n]", &nodname[i],&cnode[1][i],&cnode[2][i],&cnode[3][i]);
//boundary nodes
	fscanf(ifp,"%i%*[^\n]", &nnodbc);
	fgets(bb,max,ifp);
	fgets(bb,max,ifp);
	bcnodname = ivector(1,nnodbc);
	bcnod = ivector(1,nnodbc);
	bctyp = ivector(1,nnodbc);
	bcpotqj = vector(1,nnodbc);
	for(inodbc=1; inodbc<=nnodbc; inodbc++)	fscanf(ifp,"%i %i %f%*[^\n]", &bcnodname[inodbc],&bctyp[inodbc],&bcpotqj[inodbc]);
	fclose(ifp);

	//Read parameters for drawing segments and nodes
	xsl0 = vector(1, 3);
	xsl1 = vector(1, 3);
	xsl2 = vector(1, 3);
	ifp = fopen("ContourParams.dat", "r");
	fscanf(ifp, "%f %f %f%*[^\n]", &xsl0[1], &xsl0[2], &xsl0[3]);
	fscanf(ifp, "%f %f %f%*[^\n]", &xsl1[1], &xsl1[2], &xsl1[3]);
	fscanf(ifp, "%f %f %f%*[^\n]", &xsl2[1], &xsl2[2], &xsl2[3]);
	fclose(ifp);
}