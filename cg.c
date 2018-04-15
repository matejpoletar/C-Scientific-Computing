#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/* Primjer metode Konjugiranih gradijenata - CG  */

doublereal cg (integer n, doublereal *a, doublereal *b, doublereal *x0, doublereal tol){

	int i,j;
	doublereal *r; 
	r=malloc(n*sizeof (doublereal));

	for (i=0; i<n; i++) 
		r[i]=b[i];

	char trans ='n'; 
	doublereal alpha=-1, beta=1; 
	integer incx=1;

	dgemv_(&trans, &n, &n, &alpha, a, &n, x0, &incx, &beta, r, &incx);

	doublereal *dd; 
	dd=malloc(n*sizeof(doublereal));
	doublereal *d; 
	d=malloc(n*sizeof(doublereal));

	for (i=0;i<n;i++) 
		d[i]=r[i];

	doublereal rr, alpha0, beta1;
	char norm='f'; 
	integer m=1; 
	doublereal norm_r;

	doublereal *work; 
	work=malloc(n*n*sizeof(doublereal));

	doublereal norm_b=dlange_(&norm,&n,&m,b,&n,work);

	do{
		rr=ddot_(&n,r,&incx,r,&incx);
		alpha=1; 
		beta=0;
		dgemv_(&trans, &n, &n, &alpha, a, &n, d, &incx, &beta, dd, &incx);
		doublereal dAd = ddot_(&n,d,&incx,dd,&incx);
		alpha0=rr/dAd;

		daxpy_(&n, &alpha0, d, &incx, x0, &incx);
		alpha0=-alpha0;

		daxpy_(&n, &alpha0, dd, &incx, r, &incx);

		beta1=ddot_(&n,r,&incx,r,&incx)/rr;
		beta1=1/beta1;

		daxpy_(&n, &beta1, r, &incx, d, &incx);

		beta1=1/beta1;
		dscal_(&n, &beta1, d, &incx);
		norm_r=dlange_(&norm,&n,&m,r,&n,work);

		for (i=0; i<n; i++) 
			printf ("%f ", x0[i]); 
		printf ("\n");
		} while (norm_r/norm_b > tol);

}

int main(integer argc, char *argv[]) {
	integer n=4, i ,j;
	doublereal a[]={101,-4,8,12,-4,20,-7,3,8,-7,78,32,12,3,32,113};
	doublereal x0[]={0,0,0,0};
	doublereal x[]={1,1,1,1};
	doublereal b[4];
	for (i=0;i<n;i++) {b[i]=0;
		for (j=0;j<n;j++) b[i]+=a[i+j*n]*x[j];}
	doublereal tol=1e-8;
	cg (n,a,b,x0,tol);
	return 0;
}
