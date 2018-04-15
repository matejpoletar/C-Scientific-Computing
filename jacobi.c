#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/* 

========= SVOJSTVENI PROBLEM - Jacobijeva metoda ==========

*/

void jacobi_sd (int n, double A[4][4], double tol) {
int p,q,k,i,j,sign;
double tau,t,c,s,app,apq,aqq,pom,S;
S=0;
for (i=0;i<n;i++) 
	for (j=0;j<n;j++) if(i!=j) S+=pow(A[i][j],2);
S=sqrt(S);

while (S>tol){
	S=0;
	for (i=0;i<n;i++) 
		for (j=0;j<n;j++) if(i!=j) S+=pow(A[i][j],2);
		S=sqrt(S);	
		for (p=0;p<n-1;p++){
			for (q=p+1;q<n;q++){
				if (A[p][q]!=0){
					tau=(A[q][q]-A[p][p])/(2*A[p][q]);
					if(tau>0) 
						sign =1; 
					else 
						sign=-1;

					t=sign/(abs(tau)+sqrt(1+pow(tau,2)));
					c=1/sqrt(1+pow(t,2));
					s=t*c;

					app=A[p][p];apq=A[p][q];aqq=A[q][q];
					app=app-t*apq;
					aqq=aqq+t*apq;

					for (k=0;k<n;k++){
						pom=A[k][p];
						A[k][p]=c*pom-s*A[k][q];
						A[k][q]=s*pom+c*A[k][q];
						A[p][k]=A[k][p];A[q][k]=A[k][q];
					}
						A[p][q]=0;A[q][p]=0;
						A[p][p]=app; A[q][q]=aqq;
						}			
				}
			}
		}
}

int main (int argc, char *argv[]) {
	int i,j,k;
	integer n=4;
	double A[4][4],du[16], d[]={-10,-5,0.1,0.2};
	double tol=4*pow(10,-16);

	doublereal *u; 
	u=malloc(n*n*sizeof(doublereal));
	doublereal *tau; 
	tau=malloc(n*sizeof(doublereal));
	doublereal *wrk; 
	wrk=malloc(n*sizeof(doublereal));

	/*Generiranje slucajne matrice u*/
	integer idist=1; integer iseed[]={1,3,5,7};
	integer m=n*n;
	dlarnv_(&idist, iseed, &m, u);

	/*Racunanje QR faktorizacije*/
	integer info;
	dgeqrf_(&n,&n,u,&n,tau,wrk,&n,&info);
	dorgqr_(&n,&n,&n,u,&n,tau,wrk,&n,&info); /*u je ortogonalna */

	/*Generiranje matrice A=UDU'*/
	for (i=0;i<n;i++)
		for (j=0;j<n;j++) du[i+j*n]=u[i+j*n]*d[j];
	for (k=0;k<n;k++)
		for (i=0;i<n;i++) {
			A[i][k]=0;
		for (j=0;j<n;j++) 
			A[i][k]+=u[i+j*n]*du[k+j*n];}
	doublereal *a; 
	a=malloc(n*n*sizeof(doublereal));

	/*Kopira polje A[4][4] u a i izracunava F-normu*/
	for (i=0;i<n;i++)
		for (j=0;j<n;j++) a[i+j*n]=A[i][j];
	char norm='f';
	doublereal norma;
	norma=dlange_(&norm,&n,&n,a,&n,wrk);

	tol=tol*norma;
	jacobi_sd (n,A,tol);

	printf ("Svojstvene vrijednosti su: ");
	for (i=0;i<n;i++) printf ("%.2f  ", A[i][i]); 
	printf ("\n");

	return 0;
}
