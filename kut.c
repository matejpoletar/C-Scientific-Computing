/*Računanje glavnih kutova i glavnih vektora između dva potprostora*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main (int argc, char *argv[]) {
	doublereal a[16]={1,0,0,0,1,1,1,1};
	doublereal b[16]={1,-1,1,-1,0,1,0,1};
	int i,j;
	integer m=4,n=2;
	doublereal *tau; 
	tau=malloc(m*sizeof(doublereal));
	doublereal *work; 
	work=malloc(m*sizeof(doublereal));
	integer info;

	dgeqrf_(&m,&n,a,&m,tau,work,&m,&info);
	dorgqr_(&m,&m,&m,a,&m,tau,work,&m,&info);

	/*for (i=0;i<m;i++) {
		for (j=0;j<m;j++) printf ("%f ", a[i+j*m]); printf ("\n"); } printf ("\n");*/

	dgeqrf_(&m,&n,b,&m,tau,work,&m,&info);
	dorgqr_(&m,&m,&m,b,&m,tau,work,&m,&info);

	for (i=0;i<m;i++) {
		for (j=0;j<m;j++) printf ("%f ", b[i+j*m]); printf ("\n"); }

	char transa='n',transb='n';
	doublereal c[16]={0};
	doublereal alpha=1, beta=0;
	dgemm_(&transa,&transb,&m,&m,&m,&alpha,a,&m,b,&m,&beta,c,&m);

	for (i=0;i<m;i++) {
		for (j=0;j<m;j++) printf ("%f ", c[i+j*m]); printf ("\n"); }

	char jobu='A',jobvt='A';

	doublereal *s; 
	s=malloc(m*sizeof(doublereal));
	doublereal *u; 
	u=malloc(m*m*sizeof(doublereal));
	doublereal *vt; 
	vt=malloc(m*m*sizeof(doublereal));

	integer ldw=5*m;
	doublereal *wrk; wrk=malloc(ldw*sizeof(doublereal));
	dgesvd_(&jobu, &jobvt, &m, &m, c, &m,s,u,&m,vt,&m,wrk,&ldw,&info);

	for (j=0;j<m;j++) 
		printf ("%f ", s[j]); 
	printf ("\n");

return 0;
}
