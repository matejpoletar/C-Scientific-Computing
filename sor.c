#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/* SOR metoda za linearne sustave 
primjer: Stieltjesova matrica */

doublereal sor_norma(char norm, doublereal *a, integer n, doublereal *M, doublereal *T){
	char U='U', L='L'; integer n1=n, i,j;
	doublereal omega=1.05;
	dlacpy_(&L,&n1,&n1,a,&n1,M,&n1);
	dlacpy_(&U,&n1,&n1,a,&n1,T,&n1); /* zasada T=N*/ 
	for (i=0;i<n1;i++) 
		M[i*n1+i]/=omega;
	for (i=0;i<n1;i++) 
		T[i*n1+i]*=((omega-1)/omega);
	for (j=0;j<n1;j++)
		for(i=j;i<n1;i++)
			T[i*n1+j]=-T[i*n1+j];

	/*for(j=0;j<n1;j++){
			for(i=0;i<n1;i++){
				printf("%.2f ",M[i*n1+j]);}
			printf("\n");}
	for(j=0;j<n1;j++){
			for(i=0;i<n1;i++){
				printf("%.2f ",T[i*n1+j]);}
			printf("\n");} */
 
	char transa='n'; 
	doublereal alpha=1.00;
	dtrsm_(&L,&L,&transa,&transa,&n1,&n1,&alpha,M,&n1,T,&n1);

	/*for(j=0;j<n1;j++){
			for(i=0;i<n1;i++){
				printf("%.2f ",T[i*n1+j]);}
			printf("\n");} */

	doublereal *work;
	work=malloc(n*n*sizeof(doublereal));
	doublereal norma=dlange_(&norm,&n1,&n1,T,&n1,work);

	return norma;
}

int sor_konvergencija (doublereal norma){
	if (norma<1) 
		return 1;
	else
		return 0;
}

doublereal sor_rjesavac(doublereal *M,doublereal *T, doublereal *x0, integer n, doublereal *c, doublereal eps,doublereal *a, doublereal *b){

/*racunanje cSOR*/
char L='L', transa='n'; 
doublereal alpha=1.00; 
integer b_columns=1, n1=n, i, j;
dtrsm_(&L,&L,&transa,&transa,&n1,&b_columns,&alpha,M,&n1,c,&n1);
 
doublereal *x1; 
x1=malloc(n*sizeof(doublereal));

for (i=0;i<n;i++) x1[i]=0;

for (i=0;i<n;i++){
	for(j=0;j<n;j++) x1[i]+=T[j*n+i]*x0[j];
	x1[i]+=c[i];
}

doublereal *nor; 
nor=malloc(n*sizeof(doublereal));
doublereal nor_b=0;
doublereal nor_a;
for (i=0;i<n;i++) 
	nor_b+=b[i]*b[i]; 
nor_b=sqrt(nor_b);

int l=1;
while((l==1 || (nor_a/nor_b)>=eps) ){
	l++;
	for(i=0;i<n;i++) {x0[i]=x1[i]; x1[i]=0;}
for (i=0;i<n;i++){
		for(j=0;j<n;j++) x1[i]+=T[j*n+i]*x0[j];
		x1[i]+=c[i];
}

for (i=0;i<n;i++) nor[i]=0;
for (i=0;i<n;i++){
		for(j=0;j<n;j++) nor[i]+=a[j*n+i]*x1[j];
		nor[i]=b[i]-nor[i];
}
nor_a=0;

for (i=0;i<n;i++) 
	nor_a+=nor[i]*nor[i]; 
	nor_a=sqrt(nor_a);
}
printf ("Potrebno je %d iteracija za dostizanje tocnosti %f \n Aproksimacija rjesenja: \n x=\n ",l, eps); 
for(i=0;i<n;i++) printf ("%.6f\n", x1[i]);

}

int main(integer argc, char *argv[]) {
char norm;
printf ("Koja norma? (norma1=1, norma_beskonacno=I, Frobeniusova norma=F) \n");
scanf ("%c", &norm);
integer n=100, i ,j;
doublereal eps=1e-5;
doublereal *M, *T;
M=malloc(n*n*sizeof(doublereal));
T=malloc(n*n*sizeof(doublereal));

/* Učitavanje matrice A 100x100*/
FILE *f;
double *a;
a=malloc(n*n*sizeof(double));
f=fopen("stieltjes_matr.txt","r");
for (i=0;i<n*n;i++){
fscanf(f,"%lf",a+i);
}
fclose(f);

/*izračun b*/
doublereal x[100];
for (i=0;i<n;i++) x[i]=1;
doublereal *b;
b=malloc(n*sizeof(doublereal));
for (i=0;i<n;i++){b[i]=0;
for (j=0;j<n;j++) b[i]+=a[j*n+i]*x[i];}

doublereal *c;
c=malloc(n*sizeof(doublereal));
for(i=0;i<n;i++) c[i]=b[i];

/*generiranje nulte iteracije*/
integer idist=1;
integer iseed[]={1,2,3,5};
doublereal *x0;
x0=malloc(n*sizeof(doublereal));
dlarnv_(&idist, iseed, &n, x0);
	
doublereal normaT=sor_norma(norm,a,n,M,T);
printf("Norma operatora T je %.2f \n", normaT);
if (sor_konvergencija (normaT)) {printf ("Metoda konvergira. \n");
sor_rjesavac (M,T,x0,n,c,eps,a,b);} 
else printf ("Metoda ne konvergira \n"); 
}
