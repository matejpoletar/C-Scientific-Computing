/*Problem najmanjih kvadrata
	|| A * x - b || -> min

Rješavanje QR faktorizacijom
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main (int argc, char *argv[]) {
	int i,j;
	double a[20];
	double b[]={3.5,4.9,6.8,9.3,10.9,13.4,15.1,16.7,19,21.2};
	for(i=0;i<10;i++) 
		a[i]=1;
	for(i=10;i<20;i++) 
		a[i]=i-9;
	integer n=10, m=2;

	doublereal *tau; 
	tau=malloc(n*sizeof(doublereal));
	doublereal *work; 
	work=malloc(n*sizeof(doublereal));
	integer info;

	dgeqrf_(&n,&m,a,&n,tau,work,&n,&info);
	doublereal r[4]={0};
	r[0]=a[0];
	r[2]=a[10];
	r[3]=a[11];

	printf ("R1=\n");
	for (j=0;j<m;j++) {
		for(i=0;i<m;i++) printf ("%f ", r[j+m*i]); 
	printf ("\n");
	}

	dorgqr_(&n,&m,&m,a,&n,tau,work,&m,&info);
	printf ("Q=\n");
	for (j=0;j<n;j++) {
		for(i=0;i<m;i++) printf ("%f ", a[j+n*i]); 
		printf ("\n");
	}

	char transa='n',transb='t',side='l',uplo='u',diag='n';
	integer incx=1;
	doublereal alpha=1, beta=0;
	doublereal x[4];
	dgemv_(&transb,&n,&m,&alpha,a,&n,b,&incx,&beta,x,&incx);
	dtrsm_(&side,&uplo,&transa,&diag,&m,&m,&alpha,r,&m,x,&m);
	printf ("x=\n");
	for (j=0;j<m;j++) {
		printf ("%f ", x[j]); 
	printf ("\n");
	}

	printf ("Aproksimativni pravac je p(x)=%.4f+%.4fx\n",x[0],x[1]);
	return 0;
}
