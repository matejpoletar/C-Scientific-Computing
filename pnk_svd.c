/*
Problem najmanjih kvadrata
	|| A * x - b || -> min

Rješavanje SVD dekompozicijom 
x = V \Sigma^{-1} U^T b
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

	char jobu='S',jobvt='S';
	doublereal *s; 
	s=malloc(m*sizeof(doublereal));
	doublereal *u; 
	u=malloc(n*m*sizeof(doublereal));
	doublereal *vt; 
	vt=malloc(m*m*sizeof(doublereal));

	integer ldw=5*n,info;
	doublereal *wrk; 
	wrk=malloc(ldw*sizeof(doublereal));
	dgesvd_(&jobu, &jobvt, &n, &m, a, &n,s,u,&n,vt,&m,wrk,&ldw,&info);
	
	printf ("diag(Sigma) = \n");
	for (j=0;j<m;j++) 
		printf ("%f ", s[j]); 
	printf ("\n\n U =\n");
	
	for (i=0;i<n;i++) {
		for(j=0;j<m;j++) 
			printf ("%f ", u[j*n+i]); 
		printf ("\n");
		}
	printf ("\n V =\n");

	for (i=0;i<m;i++) {
		for(j=0;j<m;j++) 
			printf ("%f ", vt[i+j*m]); 
		printf ("\n");
		}
	
	for (i=10;i<20;i++) 
		u[i]=-u[i];

	for (i=2;i<4;i++) 
		vt[i]=-vt[i];

	for (i=0;i<m;i++) 
		for(j=0;j<m;j++) 
			vt[j+i*m] /= s[i];

	char transa='n',transb='t';
	integer incx=1;
	doublereal uvt[100]={0};
	doublereal alpha=1, beta=0;
	dgemm_(&transa,&transb,&m,&n,&m,&alpha,vt,&m,u,&n,&beta,uvt,&n);
	doublereal x[2]={0};
	dgemv_(&transa,&m,&n,&alpha,uvt,&n,b,&incx,&beta,x,&incx);

	printf ("\n x=\n");
	for(j=0;j<m;j++) 
		printf ("%f ", x[j]); 
	printf ("\n");
	printf ("\nAproksimirajući pravac p(x) = %f x + %f\n", x[0], x[1]);

	return 0;
}
