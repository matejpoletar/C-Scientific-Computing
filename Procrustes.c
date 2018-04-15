/*
Rješavanje ortogonalnog Procrustes problema primjenom SVD

||A Q - B ||_F -> min za Q \in Orth
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

int main (int argc, char *argv[]) {
doublereal a[]={83.5,102.2,99.6,122.8};
int i,j;
integer n=2;
char jobu='A',jobvt='A';

doublereal *s; 
s=malloc(n*sizeof(doublereal));
doublereal *u; 
u=malloc(n*n*sizeof(doublereal));
doublereal *vt; 
vt=malloc(n*n*sizeof(doublereal));

integer ldw=5*n,info;
doublereal *work; 
work=malloc(ldw*sizeof(doublereal));
dgesvd_(&jobu, &jobvt, &n, &n, a, &n,s,u,&n,vt,&n,work,&ldw,&info);

char trans='n';
doublereal alpha=1;
doublereal q[4]={0};
dgemm_(&trans,&trans,&n,&n,&n,&alpha,u,&n,vt,&n,&alpha,q,&n);

printf ("Q=\n");
for (i=0;i<n;i++) {
	for (j=0;j<n;j++) printf ("%f ", q[i+j*n]); 
printf ("\n"); 
}

doublereal a1[]={1.2,2.9,5.2,6.8,2.1,4.3,6.1,8.1};
doublereal b[]={1,3,5,7,2,4,6,8};
integer m=4;
doublereal beta=-1;
dgemm_(&trans,&trans,&m,&n,&n,&alpha,a1,&m,q,&n,&beta,b,&m);

for (i=0;i<n;i++) {
	for (j=0;j<n;j++) printf ("%f ", b[i+j*n]); 
printf ("\n"); 
}

char norm='f';
doublereal norma=dlange_(&norm,&m,&n,b,&m,work);
//printf ("norma = %f\n", norma);

return 0;
}
