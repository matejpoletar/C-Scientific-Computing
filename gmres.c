#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "f2c.h"
#include "fblaswr.h"
#include "clapack.h"

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;
static doublereal c_b20 = 0.;

doublereal *a; 
void matvec (doublereal *alpha, doublereal *x, doublereal *beta, doublereal *y ){
char trans='n';
integer n=100, incx=1;
dgemv_(&trans,&n,&n,alpha,a,&n,x,&incx,beta,y,&incx);}

void psolve (doublereal *x, doublereal *b){
integer n=100,incx=1;
dcopy_(&n,b,&incx,x,&incx);
}

int gmres_(n, b, x, restrt, work, ldw, h, ldh, iter, resid, matvec, psolve, 
           info)
   integer *n, *restrt, *ldw, *ldh, *iter, *info;
   doublereal *b, *x, *work, *h, *resid;
   int (*matvec) (), (*psolve) ();
{
    /* System generated locals */
    integer work_dim1, work_offset, h_dim1, h_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int drot_();
    static doublereal bnrm2;
    extern doublereal dnrm2_();
    static integer i, k, r, s, v, w, y;
    extern /* Subroutine */ int dscal_(), basis_(), dcopy_(), drotg_();
    static integer maxit;
    static doublereal rnorm, aa, bb;
    static integer cs, av, sn;
    extern /* Subroutine */ int update_();
    static doublereal tol;

    /* Parameter adjustments */
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    work_dim1 = *ldw;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --x;
    --b;

    /* Executable Statements */
    *info = 0;

/*     Test the input parameters. */

    if (*n < 0) {
	*info = -1;
    } else if (*ldw < max(1,*n)) {
	*info = -2;
    } else if (*iter <= 0) {
	*info = -3;
    } else if (*ldh < *restrt + 1) {
	*info = -4;
    }
    if (*info != 0) {
	return 0;
    }

    maxit = *iter;
    tol = *resid;

/*     Alias workspace columns. */

    r = 1;
    s = r + 1;
    w = s + 1;
    y = w;
    av = y;
    v = av + 1;

/*     Store the Givens parameters in matrix H. */

    cs = *restrt + 1;
    sn = cs + 1;

/*     Set initial residual (AV is temporary workspace here). */

    dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
    if (dnrm2_(n, &x[1], &c__1) != 0.) {

/*        AV is temporary workspace here. */

	dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
	(*matvec)(&c_b7, &x[1], &c_b8, &work[av * work_dim1 + 1]);
    }
    (*psolve)(&work[r * work_dim1 + 1], &work[av * work_dim1 + 1]);
    bnrm2 = dnrm2_(n, &b[1], &c__1);
    if (bnrm2 == 0.) {
	bnrm2 = 1.;
    }
    if (dnrm2_(n, &work[r * work_dim1 + 1], &c__1) / bnrm2 < tol) {
	goto L70;
    }

    *iter = 0;

L10:

    i = 0;

/*        Construct the first column of V. */

    dcopy_(n, &work[r * work_dim1 + 1], &c__1, &work[v * work_dim1 + 1], &
	    c__1);
    rnorm = dnrm2_(n, &work[v * work_dim1 + 1], &c__1);
    d__1 = 1. / rnorm;
    dscal_(n, &d__1, &work[v * work_dim1 + 1], &c__1);

/*        Initialize S to the elementary vector E1 scaled by RNORM. */

    work[s * work_dim1 + 1] = rnorm;
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	work[k + s * work_dim1] = 0.;
/* L20: */
    }

L30:

    ++i;
    ++(*iter);

    (*matvec)(&c_b8, &work[(v + i - 1) * work_dim1 + 1], &c_b20, &work[av * 
	    work_dim1 + 1]);
    (*psolve)(&work[w * work_dim1 + 1], &work[av * work_dim1 + 1]);

/*           Construct I-th column of H orthnormal to the previous */
/*           I-1 columns. */

    basis_(&i, n, &h[i * h_dim1 + 1], &work[v * work_dim1 + 1], ldw, &work[w *
	     work_dim1 + 1]);

/*           Apply Givens rotations to the I-th column of H. This */
/*           "updating" of the QR factorization effectively reduces */
/*           the Hessenberg matrix to upper triangular form during */
/*           the RESTRT iterations. */

    i__1 = i - 1;
    for (k = 1; k <= i__1; ++k) {
	drot_(&c__1, &h[k + i * h_dim1], ldh, &h[k + 1 + i * h_dim1], ldh, &h[
		k + cs * h_dim1], &h[k + sn * h_dim1]);
/* L40: */
    }

/*           Construct the I-th rotation matrix, and apply it to H so that
 */
/*           H(I+1,I) = 0. */

    aa = h[i + i * h_dim1];
    bb = h[i + 1 + i * h_dim1];
    drotg_(&aa, &bb, &h[i + cs * h_dim1], &h[i + sn * h_dim1]);
    drot_(&c__1, &h[i + i * h_dim1], ldh, &h[i + 1 + i * h_dim1], ldh, &h[i + 
	    cs * h_dim1], &h[i + sn * h_dim1]);

/*           Apply the I-th rotation matrix to [ S(I), S(I+1) ]'. This */
/*           gives an approximation of the residual norm. If less than */
/*           tolerance, update the approximation vector X and quit. */

    drot_(&c__1, &work[i + s * work_dim1], ldw, &work[i + 1 + s * work_dim1], 
	    ldw, &h[i + cs * h_dim1], &h[i + sn * h_dim1]);
    *resid = (d__1 = work[i + 1 + s * work_dim1], abs(d__1)) / bnrm2;
int iteracija=*iter;
printf ("%d.iter: %lf \n",iteracija, *resid);
    if (*resid <= tol) {
	update_(&i, n, &x[1], &h[h_offset], ldh, &work[y * work_dim1 + 1], &
		work[s * work_dim1 + 1], &work[v * work_dim1 + 1], ldw);
	goto L70;
    }
    if (*iter == maxit) {
	goto L50;
    }
    if (i < *restrt) {
	goto L30;
    }

L50:

/*        Compute current solution vector X. */

    update_(restrt, n, &x[1], &h[h_offset], ldh, &work[y * work_dim1 + 1], &
	    work[s * work_dim1 + 1], &work[v * work_dim1 + 1], ldw);

/*        Compute residual vector R, find norm, then check for tolerance. 
*/
/*        (AV is temporary workspace here.) */

    dcopy_(n, &b[1], &c__1, &work[av * work_dim1 + 1], &c__1);
    (*matvec)(&c_b7, &x[1], &c_b8, &work[av * work_dim1 + 1]);
    (*psolve)(&work[r * work_dim1 + 1], &work[av * work_dim1 + 1]);
    work[i + 1 + s * work_dim1] = dnrm2_(n, &work[r * work_dim1 + 1], &c__1);
    *resid = work[i + 1 + s * work_dim1] / bnrm2;
iteracija=*iter;
printf ("%d.iter: %f  \n",iteracija, *resid);
    if (*resid <= tol) {
	goto L70;
    }
    if (*iter == maxit) {
	goto L60;
    }

/*        Restart. */

    goto L10;

L60:

/*     Iteration fails. */

    *info = 1;
    return 0;

L70:

/*     Iteration successful; return. */

    return 0;

/*     End of GMRES */

} /* gmres_ */


/*     =============================================================== */
/* Subroutine */ int update_(i, n, x, h, ldh, y, s, v, ldv)
integer *i, *n;
doublereal *x, *h;
integer *ldh;
doublereal *y, *s, *v;
integer *ldv;
{
    /* System generated locals */
    integer h_dim1, h_offset, v_dim1, v_offset;

    /* Local variables */
    extern /* Subroutine */ int dgemv_(), dcopy_(), dtrsv_();



/*     This routine updates the GMRES iterated solution approximation. */


/*     .. Executable Statements .. */

/*     Solve H*Y = S for upper triangualar H. */

    /* Parameter adjustments */
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --s;
    --y;
    h_dim1 = *ldh;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    --x;

    /* Function Body */
    dcopy_(i, &s[1], &c__1, &y[1], &c__1);
    dtrsv_("UPPER", "NOTRANS", "NONUNIT", i, &h[h_offset], ldh, &y[1], &c__1 );

/*     Compute current solution vector X = X + V*Y. */

    dgemv_("NOTRANS", n, i, &c_b8, &v[v_offset], ldv, &y[1], &c__1, &c_b8, &x[
	    1], &c__1 );

    return 0;

} /* update_ */


/*     ========================================================= */
/* Subroutine */ int basis_(i, n, h, v, ldv, w)
integer *i, *n;
doublereal *h, *v;
integer *ldv;
doublereal *w;
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    doublereal d__1;

    /* Local variables */
    extern doublereal ddot_(), dnrm2_();
    static integer k;
    extern /* Subroutine */ int dscal_(), dcopy_(), daxpy_();



/*     Construct the I-th column of the upper Hessenberg matrix H */
/*     using the Gram-Schmidt process on V and W. */


    /* Parameter adjustments */
    --w;
    v_dim1 = *ldv;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --h;

    /* Function Body */
    i__1 = *i;
    for (k = 1; k <= i__1; ++k) {
	h[k] = ddot_(n, &w[1], &c__1, &v[k * v_dim1 + 1], &c__1);
	d__1 = -h[k];
	daxpy_(n, &d__1, &v[k * v_dim1 + 1], &c__1, &w[1], &c__1);
/* L10: */
    }
    h[*i + 1] = dnrm2_(n, &w[1], &c__1);
    dcopy_(n, &w[1], &c__1, &v[(*i + 1) * v_dim1 + 1], &c__1);
    d__1 = 1. / h[*i + 1];
    dscal_(n, &d__1, &v[(*i + 1) * v_dim1 + 1], &c__1);

    return 0;

}

int main(integer argc, char *argv[]) {
integer n=100, i ,j, k;
doublereal tol=pow(10,-5);
doublereal B[10000];
doublereal I[10000]={0};
a=malloc(n*n*sizeof(doublereal));
doublereal x0[100]={0};
doublereal x[100];
doublereal *q; q=malloc(n*n*sizeof(doublereal));
doublereal *tau; tau=malloc(n*sizeof(doublereal));
doublereal *wrk; wrk=malloc(n*sizeof(doublereal));

/*Generiranje slucajne matrice q*/
integer idist=1; integer iseed[]={1,3,5,7};
integer m=n*n;
dlarnv_(&idist, iseed, &m, q);

/*Racunanje QR faktorizacije*/
integer info;
dgeqrf_(&n,&n,q,&n,tau,wrk,&n,&info);

char side, trans;
side='R'; trans='T';
dorgqr_(&n,&n,&n,q,&n,tau,wrk,&n,&info); /*q je ortogonalna */

/*Racunanje g(i)*/
doublereal g[100]; doublereal b[100]={0};
for(i=0;i<n;i++) g[i]=sqrt(pow(100-i,2)-pow(99-i,2));

/*Racunanje vektora b*/ 
for (i=0;i<n;i++) 
	for (j=0;j<n;j++) b[i]+=q[i+n*j]*g[j];

/*Generiranje matrice B*/
for (i=0;i<n;i++) B[i]=b[i];
for (i=1;i<n;i++) 
	for (j=0;j<n;j++) B[i*n+j]=q[i*n-n+j];

/*Generiranje matrice A=BI*inv(B)*/
doublereal BI[10000];
for (i=0;i<n-1;i++) 
	for (j=0;j<n;j++) BI[i*n+j]=q[i*n+j];
for (i=0;i<n;i++) BI[i+9900]=b[i];

/*Racunanje inverza B dgetrf+dgetri*/
integer *ipiv; ipiv=malloc(n*n*sizeof(integer));
dgetrf_(&n,&n,B,&n,ipiv,&info);
dgetri_(&n,B,&n,ipiv,wrk,&n,&info); /*B je inverz originalne matrice B*/

for (k=0;k<n;k++) 
for (i=0;i<n;i++) {a[i+k*n]=0;
	for (j=0;j<n;j++) a[i+k*n]+=BI[i+j*n]*B[k*n+j];}

integer restrt=n;
doublereal *work; work=malloc(n*n*sizeof(doublereal));
integer ldh=restrt+2;
doublereal *h; h=malloc(ldh*ldh*sizeof(doublereal));
integer iter=n;

gmres_(&n,b,x0,&restrt,work,&n,h,&ldh,&iter,&tol,matvec,psolve,&info);

integer ilo=1,ihi=n;
dgehrd_(&n,&ilo,&ihi,a,&n,tau,work,&n,&info);
char job='E',compz='N';
doublereal *z;
doublereal *wr; wr=malloc(n*sizeof(doublereal));
doublereal *wi; wi=malloc(n*sizeof(doublereal));
dhseqr_(&job,&compz,&n,&ilo,&ihi,a,&n,wr,wi,z,&n,work,&n,&info); 

printf ("\n");
double fi;
for (i=0;i<n;i++) printf ("%ld. %f + %fi \n",i+1,wr[i],wi[i]);
for (i=0;i<n;i++) {
fi=atan(wi[i]/wr[i]);
printf ("%ld. %f + %fi \n",i+1,cos(100*fi),sin(100*fi));
}


}

