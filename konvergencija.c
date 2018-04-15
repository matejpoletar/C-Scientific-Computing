/*
==============================================================
Završni zadatak:
		USPOREDBA BRZINE KONVERGENCIJE 
		JACOBIJEVE I GAUSS-SEIDELOVE METODE
==============================================================

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

void Jacobi (integer n, doublereal *a, doublereal *b, doublereal *x) {

int i,j,k;
doublereal *x0; 

x0=malloc(n*sizeof(doublereal));
doublereal *e; 
e=malloc(n*sizeof(doublereal));
doublereal *work; 
work=malloc(n*sizeof(doublereal));

double p,q;
integer m=1;
char norm='i';
doublereal norma;

for (k=0;k<100;k++){
	for (i=0;i<n;i++) x0[i]=x[i];

for(i=0;i<n;i++){
	x[i]=b[i]; 
	for (j=0;j<i;j++)
		x[i]-=a[i+j*n]*x0[j];
	for (j=i+1;j<n;j++)
		x[i]-=a[i+j*n]*x0[j];
	x[i]/=a[i+n*i];
	}
}

for(i=0;i<n;i++) e[i]=1-x[i];

norma=dlange_(&norm,&n,&m,e,&n,work);
q=k+1;
p=pow(norma,1/q); 
printf ("Brzina konvergencije: %f\n",p);

}

void GS (integer n, doublereal *a, doublereal *b, doublereal *x) {

int i,j,k;
doublereal *e; 
e=malloc(n*sizeof(doublereal));
doublereal *work; 
work=malloc(n*sizeof(doublereal));

double p,q;
integer m=1;
char norm='i';
doublereal norma;

for (k=0;k<100;k++){
	for(i=0;i<n;i++){
		x[i]=b[i]; 
	for (j=0;j<i;j++)
		x[i]-=a[i+j*n]*x[j];
	for (j=i+1;j<n;j++)
		x[i]-=a[i+j*n]*x[j];
	x[i]/=a[i+n*i];
	}
}

for(i=0;i<n;i++) e[i]=1-x[i];

norma=dlange_(&norm,&n,&m,e,&n,work);
q=k+1;
p=pow(norma,1/q); 
printf ("Brzina konvergencije: %f\n",p);

}

int main(integer argc, char *argv[]) {
integer n; 
int i,j,k;

printf ("Unesite red matrice (Bit će generirana slučajna matrica reda n iz uniformne distribucije) \nn=");
scanf ("%ld", &n);

doublereal *s; 
s=malloc(n*n*sizeof(doublereal)); /*slucajna matrica*/
doublereal *a; 
a=malloc(n*n*sizeof(doublereal)); /*modificirana sl.matrica=matrica sustava*/ 

integer idist=2; 
integer iseed[]={1,3,5,7};
integer n2=n*n;

dlarnv_(&idist, iseed, &n2, s);
doublereal *rj; 
rj=malloc(n*sizeof(doublereal));

for (i=0;i<n;i++) rj[i]=1;

doublereal *x; 
x=malloc(n*sizeof(doublereal));
doublereal *b; 
b=malloc(n*sizeof(doublereal));

/*for(i=0;i<n;i++) printf ("%f ", b[i]);printf("\n\n");
for (i=0;i<n;i++) 
	for (j=0;j<n;j++) printf ("%f ", a[i*n+j]); printf("\n\n"); */

int metoda;
printf ("\nKoju metodu želite testirati? 1=Jacobi 2=Gauss-Seidel?\n");
scanf ("%d", &metoda);

int mode;
printf ("\nMode dominacije (1, 2 ili 3) ?\n");
scanf ("%d", &mode); 

double suma_l, suma_u, sigma, f, sigma1, sigma2, korak;

printf ("\nUnesite raspon sigma vrijednosti (sigma>1):\nsigma(min)=");
scanf("%lf",&sigma1);

printf ("sigma(max)=");
scanf("%lf",&sigma2);

printf ("\nUnesite korak kojim će se povećavati parametar sigma:");
scanf("%lf",&korak);

int br=(sigma2-sigma1)/korak+1;

switch (mode){
	case 1:
	/*MODE=1*/
	suma_l=0; suma_u=0;
	for (i=0;i<n;i++) 
		for (j=0;j<=i;j++)
	suma_l+=abs(s[i*n+j]);
	for (i=0;i<n;i++) 
		for (j=i+1;j<n;j++)
	suma_u+=abs(s[i+j*n]);
	for(k=0;k<br;k++){
		sigma=sigma1+k*korak;
		f=sigma*suma_u/suma_l;
	for (i=0;i<n;i++) 
		for (j=0;j<n;j++) {
	if (j<=i) a[i*n+j]=s[i*n+j]*f; else a[i*n+j]=s[i*n+j];}
	printf ("sigma=%.2f  ",sigma);

	for (i=0;i<n;i++){b[i]=0;
	for (j=0;j<n;j++) b[i]+=(a[j*n+i]*rj[j]);}
	for (i=0;i<n;i++) x[i]=0;

	if (metoda==1) {Jacobi (n,a,b,x); printf ("\n");}
	else {GS (n,a,b,x); printf ("\n");}
	}break;

/*for (i=0;i<n;i++) 
	for (j=0;j<n;j++) printf ("%f ", a[i*n+j]); printf("\n\n");*/

case 2:
/*MODE=2*/
	for(k=0;k<br;k++){
	sigma=sigma1+k*korak;
	for (i=0;i<n-1;i++) {
		suma_l=0; suma_u=0;
		for (j=0;j<=i;j++) suma_l+=abs(s[i+j*n]);
		for (j=i+1;j<n;j++) suma_u+=abs(s[i+j*n]);
			f=sigma*suma_u/suma_l;
		for (j=0;j<n;j++) {
			if (j<=i) 
				a[i+j*n]=s[i+j*n]*f; 
			else 
				a[i+j*n]=s[i+j*n];}}
	for (j=0;j<n;j++) a[n-1+j*n]=s[n-1+j*n];
	printf ("sigma=%.2f ",sigma);

	for (i=0;i<n;i++){b[i]=0;
	for (j=0;j<n;j++) b[i]+=(a[j*n+i]*rj[j]);}
	for (i=0;i<n;i++) x[i]=0;

	if (metoda==1) {Jacobi (n,a,b,x); printf ("\n");}
	else {GS (n,a,b,x); printf ("\n");}
} break;

case 3:
/*MODE=3*/
	for(k=0;k<br;k++){
	sigma=sigma1+k*korak;
	for (i=0;i<n;i++) 
		for (j=0;j<n;j++) {
	if (j<=i) { 
	if (s[i+j*n]>0) f=1; else f=-1;
		a[i+j*n]=f*abs(s[j+i*n])*sigma; }
	else { 	
	if (s[i+j*n]>0) f=1; else f=-1;
		a[i+j*n]=f*abs(s[i+j*n]);}}
	printf ("sigma=%.2f ",sigma);

	for (i=0;i<n;i++){b[i]=0;
	for (j=0;j<n;j++) b[i]+=(a[j*n+i]*rj[j]);}
	for (i=0;i<n;i++) x[i]=0;
	if (metoda==1) {Jacobi (n,a,b,x); printf ("\n");}
	else {GS (n,a,b,x); printf ("\n");} 
}break;
}

return 0;
}

