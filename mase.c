#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/*SISTEM MASA NA OPRUZI */

int main (int argc, char *argv[]) {
	int i,j;
	integer n=4, info;
	double M[]={2,5,3,6};
	double K[]={24,-9,-5,0,-9,22,-8,-5,-5,-8,25,-7,0,-5,-7,18};
	double A[16];
	for (i=0;i<n;i++)
		for (j=0;j<n;j++) 
			A[i+j*n]=K[i+j*n]/(sqrt(M[i])*sqrt(M[j]));

	printf ("A=\n");
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) printf ("%f  ", A[i+j*n]); 
		printf ("\n");
		}

	char jobz='V',uplo='L';
	doublereal *w; w=malloc(n*sizeof(doublereal));
	integer lwork=3*n-1;
	doublereal *work; work=malloc(lwork*sizeof(doublereal));
	dsyev_(&jobz, &uplo, &n, A, &n, w, work, &lwork, &info);

	printf ("\nL=\n");
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) {
			if (i==j) printf ("%f ", w[i]); 
			else printf ("0.000000 "); }
		printf ("\n");
		}

	printf ("\nU=\n");
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) printf ("%f  ", A[i+j*n]); 
		printf ("\n");
		}

	return 0;
}
