#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/*
PREDKONDICIONIRANA METODA KONJUGIRANIH GRADIJENATA
*/

void cg (integer n, doublereal *a, doublereal *b, doublereal *x0, doublereal tol){
	int i,j;
	doublereal *r; 
	r=malloc(n*sizeof (doublereal));

	for (i=0;i<n;i++) r[i]=b[i];

	char trans ='n'; doublereal alpha=-1, beta=1; integer incx=1;
	dgemv_(&trans, &n, &n, &alpha, a, &n, x0, &incx, &beta, r, &incx);

	doublereal *dd; 
	dd=malloc(n*sizeof(doublereal));
	doublereal *d; 
	d=malloc(n*sizeof(doublereal));

	for (i=0;i<n;i++) d[i]=r[i];

	doublereal rr, alpha0, beta1;
	char norm='f'; 
	integer m=1;
	doublereal *work; 
	work=malloc(n*n*sizeof(doublereal));

	doublereal norm_b=dlange_(&norm,&n,&m,b,&n,work);
	doublereal norm_r=dlange_(&norm,&n,&m,r,&n,work);

	int k=0;
	while (norm_r/norm_b > tol){
		rr=ddot_(&n,r,&incx,r,&incx);
		alpha=1; beta=0;

		dgemv_(&trans, &n, &n, &alpha, a, &n, d, &incx, &beta, dd, &incx);
		doublereal dAd =ddot_(&n,d,&incx,dd,&incx);

		alpha0=rr/dAd;
		daxpy_(&n, &alpha0, d, &incx, x0, &incx);
		alpha0=-alpha0;
		daxpy_(&n, &alpha0, dd, &incx, r, &incx);

		beta1=ddot_(&n,r,&incx,r,&incx)/rr;
		doublereal *rr; 
		rr=malloc(n*sizeof(doublereal));

		dcopy_(&n,r,&incx,rr,&incx);
		daxpy_(&n, &beta1, d, &incx, rr, &incx);
		dcopy_(&n,rr,&incx,d,&incx);
		norm_r=dlange_(&norm,&n,&m,r,&n,work);
		k++;

	}
	printf ("%d \n", k); 
	for (i=0;i<n;i++) printf ("%f, ", x0[i]);
printf ("\n");
}


void ic (integer n, doublereal *a){
int i,j,k;
for (i=0; i<n; i++){
	for (k=0; k<i; k++)
		a[i+i*n]-=a[k+i*n]*a[k+i*n];
	a[i+i*n]=sqrt(a[i+i*n]);
	for (j=i+1; j<n; j++){
		if (a[i+j*n]!=0){
		for (k=0; k<i; k++)
			a[i+j*n]-=a[k+i*n]*a[k+j*n];
			a[i+j*n]/=a[i+i*n];
			}
		}
	}

}


void diag_pcg (integer n, doublereal *a, doublereal *b, doublereal *x0, doublereal tol){
	integer incx=1, incy=1;
	int i, k=0;
	doublereal alpha=-1, beta=1, alpha_k, beta_k;
	char uplo='u',trans='t',diag='n';

	doublereal *r; 
	r=malloc(n*sizeof(doublereal)); /* r je rezidual*/

	dcopy_(&n,b,&incx,r,&incy);
	dgemv_(&trans, &n, &n, &alpha, a, &n, x0, &incx, &beta, r, &incx); /*r=b-Ax*/

	doublereal *p; 
	p=malloc(n*sizeof(doublereal));

	for (i=0;i<n;i++) p[i]=r[i]/sqrt(a[i+n*i]);

	doublereal *d; 
	d=malloc(n*sizeof(doublereal));
	dcopy_(&n,p,&incx,d,&incy);

	integer l=1;
	char norm='f';

	doublereal *work; 
	work=malloc(n*n*sizeof(doublereal));

	doublereal norm_b=dlange_(&norm,&n,&l,b,&n,work);
	doublereal norm_r=dlange_(&norm,&n,&l,r,&n,work);

	doublereal *dd; 
	dd=malloc(n*sizeof(doublereal));
	doublereal *pom; 
	pom=malloc(n*sizeof(doublereal));
	doublereal rp;

while (norm_r/norm_b > tol) {
	alpha=1; beta=0;
	dgemv_(&trans, &n, &n, &alpha, a, &n, d, &incx, &beta, dd, &incx); /* racuna A*d */
	doublereal dAd=ddot_(&n,d,&incx,dd,&incx); /*Racuna dT*(A*d)*/ 
	
	rp=ddot_(&n,r,&incx,p,&incx);
	alpha_k=rp/dAd;
	daxpy_(&n, &alpha_k, d, &incx, x0, &incx);
	alpha_k=-alpha_k;
	daxpy_(&n, &alpha_k, dd, &incx, r, &incx);
	dcopy_(&n,r,&incx,p,&incy);
	for (i=0;i<n;i++) 
		p[i]/=sqrt(a[i+n*i]);
	beta_k=ddot_(&n,r,&incx,p,&incx)/rp;
	
	dcopy_(&n,p,&incx,pom,&incy);
	daxpy_(&n,&beta_k,d,&incx,pom,&incx);
	dcopy_(&n,pom,&incx,d,&incy);
	norm_r=dlange_(&norm,&n,&l,r,&n,work);
	k++;
	}
	printf ("%d \n",k); 
	for (i=0;i<n;i++) printf ("%f  ", x0[i]);
printf ("\n");
}


void pcg (integer n, doublereal *m, doublereal *a, doublereal *b, doublereal *x0, doublereal tol){
integer incx=1, incy=1; int i;
doublereal alpha=-1, beta=1, alpha_k, beta_k;
char uplo='u',trans='t',diag='n';

doublereal *r; 
r=malloc(n*sizeof(doublereal)); /* r je rezidual*/

dcopy_(&n,b,&incx,r,&incy);
dgemv_(&trans, &n, &n, &alpha, a, &n, x0, &incx, &beta, r, &incx); /*r=b-Ax*/
doublereal *p; 
p=malloc(n*sizeof(doublereal));

dcopy_(&n,r,&incx,p,&incy);
dtrsv_(&uplo,&trans,&diag,&n,m,&n,p,&incx); /*rjesava Mp=r*/

trans='n';
dtrsv_(&uplo,&trans,&diag,&n,m,&n,p,&incx);
doublereal *d; 
d=malloc(n*sizeof(doublereal));

dcopy_(&n,p,&incx,d,&incy);
integer l=1;
int k=0;
char norm='f';

doublereal *work; 
work=malloc(n*n*sizeof(doublereal));

doublereal norm_b=dlange_(&norm,&n,&l,b,&n,work);
doublereal norm_r=dlange_(&norm,&n,&l,r,&n,work);

doublereal *dd; 
dd=malloc(n*sizeof(doublereal));

doublereal *pom; 
pom=malloc(n*sizeof(doublereal));

doublereal rp;

while (norm_r/norm_b > tol) {
	alpha=1; beta=0;
	dgemv_(&trans, &n, &n, &alpha, a, &n, d, &incx, &beta, dd, &incx); /* racuna A*d */
	doublereal dAd=ddot_(&n,d,&incx,dd,&incx); /*Racuna dT*(A*d)*/ 
	rp=ddot_(&n,r,&incx,p,&incx);
	alpha_k=rp/dAd;
	daxpy_(&n, &alpha_k, d, &incx, x0, &incx);
	alpha_k=-alpha_k;
	daxpy_(&n, &alpha_k, dd, &incx, r, &incx);
	dcopy_(&n,r,&incx,p,&incy);
	trans='t';
	dtrsv_(&uplo,&trans,&diag,&n,m,&n,p,&incx); 
	trans='n';
	dtrsv_(&uplo,&trans,&diag,&n,m,&n,p,&incx);
	beta_k=ddot_(&n,r,&incx,p,&incx)/rp;
	dcopy_(&n,p,&incx,pom,&incy);
	daxpy_(&n,&beta_k,d,&incx,pom,&incx);
	dcopy_(&n,pom,&incx,d,&incy);
	norm_r=dlange_(&norm,&n,&l,r,&n,work);
	k++;
	}
	printf ("%d \n",k); 
	for (i=0;i<n;i++) printf ("%f  ", x0[i]);
printf ("\n");
}

int main(integer argc, char *argv[]) {
/* Ucitavanje matrice  A*/
	FILE *f;
	doublereal *a;
	int i,j,n=100;
	a=malloc(n*n*sizeof(doublereal));
	f=fopen("stieltjes_matr.txt","r");
		for (i=0;i<n*n;i++){
			fscanf(f,"%lf",a+i);
		}
	fclose(f);

/*pocetna iteracija x0 i egzaktno rjesenje x*/
	doublereal x0[100]={0};
	doublereal x[100]; 
	for(i=0;i<n;i++) x[i]=1;

/*Racunanje b*/ 
	doublereal b[100];
	for (i=0;i<n;i++) {b[i]=0;
		for (j=0;j<n;j++) b[i]+=a[i+j*n]*x[j];}

	doublereal tol=1e-8;
	doublereal *m; 
	m=malloc(n*n*sizeof (doublereal));

	for (i=0;i<n*n;i++) m[i]=a[i];
	ic(n,m); /* m je gornjetrokutasta, u biti R-faktor Choleskog*/
	printf("Broj iteracija:\nNeprekondicionirani sustav:");
	cg(n,a,b,x0,tol);

	for (i=0;i<n;i++) x0[i]=0;
	printf ("\nDijagonalno prekondicionirani sustav:");
	diag_pcg(n,a,b,x0,tol);

	for (i=0;i<n;i++) x0[i]=0;
	printf ("\nPrekondicionirani sustav:");
	pcg(n,m,a,b,x0,tol);

	return 0;
}

