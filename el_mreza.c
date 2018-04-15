
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/* Primjena Gauss-Seidelove i SOR metode na rješavanje 
problema električne mreže pri računanju potencijala
u čvornim točkama */


doublereal sor_norma(char norm, doublereal *a, integer n, doublereal *M, doublereal *T, doublereal omega){
	char U='U', L='L'; integer n1=n, i,j;

	dlacpy_(&L,&n1,&n1,a,&n1,M,&n1);
	dlacpy_(&U,&n1,&n1,a,&n1,T,&n1); /* zasada T=N*/ 
	for (i=0;i<n1;i++) M[i*n1+i]/=omega;
	for (i=0;i<n1;i++) T[i*n1+i]*=((omega-1)/omega);
	for (j=0;j<n1;j++)
		for(i=j;i<n1;i++)T[i*n1+j]=-T[i*n1+j];

	/*for(j=0;j<n1;j++){
			for(i=0;i<n1;i++){
				printf("%.2f ",M[i*n1+j]);}
			printf("\n");}
	for(j=0;j<n1;j++){
			for(i=0;i<n1;i++){
				printf("%.2f ",T[i*n1+j]);}
			printf("\n");} */
	 
	char transa='n'; doublereal alpha=1.00;
	dtrsm_(&L,&L,&transa,&transa,&n1,&n1,&alpha,M,&n1,T,&n1);

	/*for(j=0;j<n1;j++){
			for(i=0;i<n1;i++){
				printf("%.2f ",T[i*n1+j]);}
			printf("\n");} */

	doublereal *work;
	work=malloc(n*n*sizeof(doublereal));
	doublereal norma=dlange_(&norm,&n1,&n1,T,&n1,work);
	return norma;}

	int sor_konvergencija (doublereal norma){
	if (norma<1) return 1;
	return 0;
}

doublereal sor_rjesavac(doublereal *M,doublereal *T, doublereal *x0, integer n, doublereal *c, doublereal eps,doublereal *a, doublereal *b){

	/*racunanje cSOR*/
	char L='L', transa='n'; doublereal alpha=1.00; integer b_columns=1, n1=n, i, j;
	dtrsm_(&L,&L,&transa,&transa,&n1,&b_columns,&alpha,M,&n1,c,&n1);
	 
	doublereal *x1; 
	x1=malloc(n*sizeof(doublereal));
	for (i=0;i<n;i++) x1[i]=0;
	for (i=0;i<n;i++){
			for(j=0;j<n;j++) x1[i]+=T[j*n+i]*x0[j];
			x1[i]+=c[i];}

	doublereal *nor; nor=malloc(n*sizeof(doublereal));
	doublereal nor_b=0;
	doublereal nor_a;
	for (i=0;i<n;i++) nor_b+=b[i]*b[i]; nor_b=sqrt(nor_b);

	int l=1;
	while((l==1 || (nor_a/nor_b)>=eps) ){l++;
	for(i=0;i<n;i++) {x0[i]=x1[i]; x1[i]=0;}
	for (i=0;i<n;i++){
			for(j=0;j<n;j++) x1[i]+=T[j*n+i]*x0[j];
			x1[i]+=c[i];}

	for (i=0;i<n;i++) nor[i]=0;
	for (i=0;i<n;i++){
			for(j=0;j<n;j++) nor[i]+=a[j*n+i]*x1[j];
			nor[i]=b[i]-nor[i];}
	nor_a=0;
	for (i=0;i<n;i++) nor_a+=nor[i]*nor[i]; nor_a=sqrt(nor_a);}
	printf ("Potrebno je %d iteracija za dostizanje tocnosti %f \nAproksimacija rjesenja: \nx=(",l, eps); 
	for(i=0;i<n-1;i++) printf ("%.0f, ", x1[i]); 
	printf ("%.0f)\n", x1[n-1]);
}

int main(integer argc, char *argv[]) {
	char norm='i';
	integer n=6, i ,j;
	doublereal eps=1e-8;
	doublereal *M, *T;
	M=malloc(n*n*sizeof(doublereal));
	T=malloc(n*n*sizeof(doublereal));

	doublereal x0[6]={0};
	doublereal a[]={11,-20,0,0,0,-2,-5,41,-3,0,-3,0,0,-15,7,-1,0,0,0,0,-4,2,-10,0,0,-6,0,-1,28,-15,-1,0,0,0,-15,47};
	doublereal b[]={500,0,0,0,0,0};
	doublereal c[]={500,0,0,0,0,0};

	doublereal omega=1; printf ("Gauss-Seidel: \n");
	doublereal normaT=sor_norma(norm,a,n,M,T,omega);
	sor_rjesavac (M,T,x0,n,c,eps,a,b);

	doublereal *M1, *T1;
	M1=malloc(n*n*sizeof(doublereal));
	T1=malloc(n*n*sizeof(doublereal));

	doublereal x00[6]={0};
	doublereal a1[]={11,-20,0,0,0,-2,-5,41,-3,0,-3,0,0,-15,7,-1,0,0,0,0,-4,2,-10,0,0,-6,0,-1,28,-15,-1,0,0,0,-15,47};
	doublereal b1[]={500,0,0,0,0,0};
	doublereal c1[]={500,0,0,0,0,0};

	omega=1.35; printf ("SOR: \n");
	normaT=sor_norma(norm,a,n,M1,T1,omega);
	sor_rjesavac (M1,T1,x00,n,c1,eps,a1,b1);

	return 0;
}


