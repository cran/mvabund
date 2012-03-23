// GLM estimation 
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 16-333-2011

#include "resampTest.h"
//#include "time.h"
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>  
#include <gsl/gsl_sf_psi.h>

// Note try to use gsl functions as much as poosible to increase speed and stabilitiy

glm::glm(const reg_Method *mm)
   : mmRef(mm), Yref(NULL), Xref(NULL), Oref(NULL), 
     Beta(NULL), varBeta(NULL), Mu(NULL), Eta(NULL), Res(NULL), 
     Var(NULL), wHalf(NULL), sqrt1_Hii(NULL), phi(NULL),
     ll(NULL), dev(NULL), aic(NULL), iterconv(NULL) 
{ 
     mintol = mmRef->tol;
     lTol=-log(mintol);
     maxiter = 50;
//     printf("mintol=%.5f, lTol=%.3f\n", mintol, lTol);
     //    printf("glm constructor called.\n"); 
}

PoissonGlm::PoissonGlm(const reg_Method *mm):glm(mm)
{
//    printf("Poisson constructor called.\n");
}

LogiGlm::LogiGlm(const reg_Method *mm):PoissonGlm(mm) {
//    printf("Logi constructor called.\n");
}

NBinGlm::NBinGlm(const reg_Method *mm):PoissonGlm(mm){ 
//   printf("NB constructor called.\n");
}

glm::~glm() {
//   printf("glm destructor called.\n");
}


PoissonGlm::~PoissonGlm() {
//   printf("Poisson destructor called.\n");
}

LogiGlm::~LogiGlm() {
//   printf("Logi destructor called.\n");
}

NBinGlm::~NBinGlm() {
//   printf("NB destructor called.\n");
}

void glm::releaseGlm(void)
{ 
    if (Beta!=NULL)
       gsl_matrix_free(Beta);
    if (Mu!=NULL)
        gsl_matrix_free(Mu);
    if (Eta!=NULL)
        gsl_matrix_free(Eta);
    if (Res!=NULL)
        gsl_matrix_free(Res);
    if (Var!=NULL)
        gsl_matrix_free(Var);
    if (wHalf!=NULL)
        gsl_matrix_free(wHalf);
    if (sqrt1_Hii!=NULL)
        gsl_matrix_free(sqrt1_Hii);
    if (varBeta!=NULL)
        gsl_matrix_free(varBeta);
    if (phi!=NULL)
        delete[] phi;
    if (ll!=NULL)	
        delete[] ll;
    if (dev!=NULL)	
        delete[] dev;
    if (iterconv!=NULL)	
        delete[] iterconv;
    if (aic!=NULL)	
        delete[] aic;
//    printf("glm object released.\n");
}

void glm::initialGlm(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O) 
{ 
    releaseGlm();

    Yref = Y;
    Oref = O;
    Xref = X;

    nRows = Y->size1;
    nVars = Y->size2;
    nParams = X->size2;

    unsigned int i, j;
    phi = new double [nVars];
    ll = new double [nVars];
    dev = new double [nVars];
    aic = new double [nVars];
    iterconv = new unsigned int [nVars]; 

    //Xref = gsl_matrix_alloc(nRows, nParams);
    Beta = gsl_matrix_alloc(nParams, nVars);
    Mu = gsl_matrix_alloc(nRows, nVars);
    Eta = gsl_matrix_alloc(nRows, nVars);
    Res = gsl_matrix_alloc(nRows, nVars);
    Var = gsl_matrix_alloc(nRows, nVars);
    wHalf = gsl_matrix_alloc(nRows, nVars);
    sqrt1_Hii = gsl_matrix_alloc(nRows, nVars);
    varBeta = gsl_matrix_alloc(nParams, nVars);

    gsl_matrix_set_zero (Beta);
    gsl_matrix_set_zero (varBeta);
//  Note: setting the initial value is important
//  e.g., using mean(Y) for binomial regression doesn't work
//    gsl_matrix *t1;
//    t1 = gsl_matrix_alloc(nRows, 1);
//    gsl_matrix_set_all (t1, 1.0); // intercept
//    GetMean(t1, Y, Mu);
//    gsl_matrix_free(t1);
//
//  Use binomial$initialize: MuStart = (Y+0.5)/2 
//  It seems to work for poisson and negative.binomial as well
    gsl_matrix_memcpy(Mu, Y);
    gsl_matrix_add_constant(Mu, 0.5);
    gsl_matrix_scale(Mu, 0.5);
    for (j=0; j<nVars; j++) {
        phi[j] = 0;
        ll[j] = 0;
        dev[j] = 0;
        aic[j] = 0;
        iterconv[j] = 0;
        for (i=0; i<nRows; i++) 
            gsl_matrix_set(Eta, i, j, link(gsl_matrix_get(Mu, i, j)));
    }    
    rdf = nRows - nParams;
}


int glm::copyGlm(glm *src)
{    
    initialGlm(src->Yref, src->Xref, src->Oref);

    // copy properties
    Xref = gsl_matrix_alloc(src->nRows, src->nParams);
    gsl_matrix_memcpy(Xref, src->Xref);
    gsl_matrix_memcpy(Beta, src->Beta);
    gsl_matrix_memcpy(Mu, src->Mu);
    gsl_matrix_memcpy(Eta, src->Eta);
    gsl_matrix_memcpy(Res, src->Res);
    gsl_matrix_memcpy(Var, src->Var);
    gsl_matrix_memcpy(wHalf, src->wHalf);
    gsl_matrix_memcpy(sqrt1_Hii, src->sqrt1_Hii);
    gsl_matrix_memcpy(varBeta, src->varBeta);
    
    for (unsigned int i=0; i<nVars; i++) {
        phi[i] = src->phi[i];	
        ll[i] = src->ll[i];
        dev[i] = src->dev[i];
        iterconv[i] = src->iterconv[i];
	aic[i] = src->aic[i];
    }
    
    return SUCCESS;    
}



int PoissonGlm::EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, double *a)
{
    initialGlm(Y, X, O);

    unsigned int i, j;   
    double yij, mij, vij, wij, rij, tol, hii;
    gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *TMP = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);   
    gsl_vector *t = gsl_vector_alloc(nParams);
    gsl_vector_view wj, vj, dj, hj;

    for (j=0; j<nVars; j++) {
       if ( a!=NULL ) phi[j]=a[j]; 
       // estimate mu and beta   
       iterconv[j] = betaEst(j, maxiter, &tol, phi[j]);        
       // other properties based on mu
       wj = gsl_matrix_column (wHalf, j);
       for (i=0; i<nRows; i++) {
            mij = gsl_matrix_get(Mu, i, j);
            // get variance
            vij = varfunc( mij, phi[j] );
            gsl_matrix_set(Var, i, j, vij); 
            // get weight
            wij = weifunc(mij, phi[j]);           
            gsl_vector_set(&wj.vector, i, sqrt(wij));             
            // get (Pearson) residuals
            yij = gsl_matrix_get(Y, i, j);
            rij = (yij-mij)/sqrt(vij);
            gsl_matrix_set(Res, i, j, rij);        
            // get elementry log-likelihood	   
            ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
            dev[j] = dev[j] + devfunc( yij, mij, phi[j] );            
       }      
       aic[j]=-ll[j]+2*(nParams);

       // Get X * W^1/2
       wj = gsl_matrix_column (wHalf, j);
       for (i=0; i<nParams; i++) 
           gsl_matrix_set_col (WX, i, &wj.vector);
       gsl_matrix_mul_elements (WX, X);
       // X^T * W * X
       gsl_matrix_set_zero (XwX);
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
       gsl_linalg_cholesky_decomp (XwX); // provided XwX is non-singular        
       // Calc varBeta
       gsl_linalg_cholesky_invert (XwX);
       vj = gsl_matrix_column (varBeta, j);
       dj = gsl_matrix_diagonal (XwX);
       gsl_vector_memcpy (&vj.vector, &dj.vector);

       // Calc sqrt(1-hii)
       // hii is diagonal element of W^1/2*X*(X'WX)^-1*X^T*W^1/2
       gsl_blas_dsymm (CblasRight, CblasLower, 1.0, XwX, WX, 0.0, TMP);
       gsl_matrix_mul_elements (TMP, WX);
       gsl_vector_set_all (t, 1.0);
       hj = gsl_matrix_column (sqrt1_Hii, j);
       gsl_blas_dgemv (CblasNoTrans, 1.0, TMP, t, 0.0, &hj.vector); 
       for (i=0; i<nRows; i++) {
           hii = gsl_vector_get(&hj.vector, i);
           gsl_vector_set(&hj.vector, i, sqrt(1-MIN((1-mintol),hii))); 
       }
   } 
   // standardize perason residuals by rp/sqrt(1-hii) 
   gsl_matrix_div_elements (Res, sqrt1_Hii);
   subtractMean(Res);  // have mean subtracted

   gsl_matrix_free(XwX);
   gsl_matrix_free(WX);
   gsl_matrix_free(TMP);
   gsl_vector_free(t);

   return SUCCESS;    
}


int PoissonGlm::betaEst( unsigned int id, unsigned int iter, double *tol, double a )
{
   unsigned int i, isConv=FALSE, step=0;
   double oij, wij, zij;   
   double eij, mij, yij;   
   double dev_old, diff;
   gsl_vector *z = gsl_vector_alloc(nRows);
   gsl_vector *y_m = gsl_vector_alloc(nRows);
   gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams); 
   gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);
   gsl_vector *Xwz = gsl_vector_alloc(nParams);

   gsl_vector_view yj, mj, ej, bj, oj, wj, Xwi;
   yj=gsl_matrix_column(Yref, id);
   mj=gsl_matrix_column(Mu, id);
   ej=gsl_matrix_column(Eta, id);
   bj=gsl_matrix_column (Beta, id);   
   wj=gsl_matrix_column (wHalf, id);   
   if ( Oref!=NULL ) oj = gsl_matrix_column(Oref, id);
   // IRLS
//   unsigned int debug = FALSE;
//   double det;
   while ( isConv != TRUE ) {
       step++;
       dev_old = dev[id];
       gsl_matrix_memcpy (WX, Xref);
       for (i=0; i<nRows; i++) {
           eij = gsl_vector_get(&ej.vector, i);
           mij = gsl_vector_get(&mj.vector, i);
           yij = gsl_vector_get(&yj.vector, i);           
           if (Oref == NULL) oij = 0;
           else oij = gsl_vector_get(&oj.vector, i); 
           // update weight 
	   wij = sqrt(weifunc(mij, a));
	   // update z
	   zij = eij + (yij-mij)/rcpLinkDash(mij) - oij;
           gsl_vector_set( z, i, wij*zij ); // z = wHalf * z
           // get Xw
	   Xwi = gsl_matrix_row (WX, i);
           gsl_vector_scale (&Xwi.vector, wij); // Xw = WHalf * X 
/*           if (wij < mintol) {
              printf("wHalf(%d, %d) < mintol\n", i, id);
              debug == TRUE;
           }
*/
        }
        // XwX can be singular so use QR LSQ instead
//       invLSQ(WX, z, &bj.vector);   

       gsl_matrix_set_zero (XwX);
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);

 /*      if (debug == TRUE) {
          displaymatrix(WX, "WX");
          printf("det(XwX) = %.8f\n", calcDet(XwX));
       }
 */        
       // solve X^T * W * X * bj = X^T * W * z
       gsl_linalg_cholesky_decomp (XwX); // provided XwX is non-singular
       // X^T * W * z = (Xw)^T * z
       gsl_blas_dgemv (CblasTrans, 1.0, WX, z, 0.0, Xwz);
       gsl_linalg_cholesky_solve (XwX, Xwz, &bj.vector);
 
       // update eta = X*beta + offset
       gsl_blas_dgemv (CblasNoTrans, 1.0, Xref, &bj.vector, 0.0, &ej.vector);
       if ( Oref!=NULL ) gsl_vector_add (&ej.vector, &oj.vector);	   


       // update mu and deviance
       dev[id] = 0;
       for (i=0; i<nRows; i++) {
            eij=gsl_vector_get(&ej.vector, i);
            eij=(eij<-lTol)?-lTol:((eij>lTol)?lTol:eij);
            gsl_vector_set(&ej.vector, i, eij);
            gsl_vector_set(&mj.vector, i, invLink(eij));
            yij = gsl_vector_get(&yj.vector, i);
            dev[id] = dev[id] + devfunc(yij, mij, a);
        }        
	// Test convergence as the glm function in R
        // *tol = ABS(dev[id]-dev_old)/(ABS(dev[id])+0.1); 
        diff = dev[id]-dev_old;        
        *tol = GSL_MAX(diff, -diff)/(GSL_MAX(dev[id], -dev[id])+0.1);
        if ( (*tol < mintol) | (step == iter )) break;
   } 

   gsl_vector_free(z);
   gsl_vector_free(y_m);
   gsl_matrix_free(WX); 
   gsl_matrix_free(XwX); 
   gsl_vector_free(Xwz); 

   return step;
}


int NBinGlm::nbinfit(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O)
{   
    initialGlm(Y, X, O);

    unsigned int i, j, isConv;
    double yij, mij, vij, hii;
    double a, tol, fA, fAdash;
    double initphi=1e-4;
    gsl_vector_view b0j, m0j, e0j, v0j;
    gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *TMP = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);   
    gsl_vector *t = gsl_vector_alloc(nParams);
    gsl_vector_view wj, vj, dj, hj;

    // Get initial estimates from Poisson    
    PoissonGlm fit0( mmRef );
    fit0.initialGlm( Y, X, O );    

    for (j=0; j<nVars; j++) {  
        // Get initial beta estimates from Poisson
        fit0.betaEst(j, maxiter, &tol, 0);
	b0j = gsl_matrix_column(fit0.Beta, j);
	gsl_matrix_set_col(Beta, j, &b0j.vector);
        m0j = gsl_matrix_column(fit0.Mu, j);
	gsl_matrix_set_col(Mu, j, &m0j.vector);
        e0j = gsl_matrix_column(fit0.Eta, j);
	gsl_matrix_set_col(Eta, j, &e0j.vector);
        v0j = gsl_matrix_column(fit0.varBeta, j);
	gsl_matrix_set_col(varBeta, j, &v0j.vector);
        dev[j] = fit0.dev[j];

        // Get initial phi estimates 
        if (mmRef->estiMethod==CHI2) {
            a = fit0.getDisper(j); 
            phi[j] = (a<1)? 0:1/a;
	    isConv = (a<1)? TRUE:FALSE;
	    while (isConv != TRUE) {
                iterconv[j]++;
                // 1-step update beta
	        betaEst(j, 1, &tol, phi[j]); 
                // 1-step update phi
		phi[j] = phi[j]*getDisper(j);
	        if ((tol<mintol) | (iterconv[j]==maxiter) | (phi[j]<0)) 
	           break;
        }   }
        else if (mmRef->estiMethod==NEWTON) {
            getfAfAdash(initphi, j, &fA, &fAdash);
            a = (fA>0)?-fA/fAdash:0;
            phi[j] = (a<0)? mintol:a; 
	    isConv = (fA>0)? FALSE:TRUE;
	    while ( isConv != TRUE ) {
                iterconv[j]++;	    
                // 1-step update beta
	        betaEst(j, 1, &tol, phi[j]);
                // 1-step update phi
                getfAfAdash(phi[j], j, &fA, &fAdash);
                phi[j] = phi[j]-fA/fAdash;
                // check convergence
	        if ((tol<mintol*10) | (iterconv[j]==maxiter) | (phi[j] < 0) | (phi[j]!=phi[j])) break;
       }   }
       // restore poisson if phi[j]<0 or nan
       if ( (phi[j] < 0) | (phi[j]!=phi[j]) ) { 
//            printf("%.2f, fA=%.2f, fAdash=%.2f\n", phi[j], fA, fAdash);
            phi[j]=0;
            gsl_matrix_set_col (Beta, j, &b0j.vector);
            gsl_matrix_set_col (Mu, j, &m0j.vector);
            gsl_matrix_set_col (Eta, j, &e0j.vector);
            gsl_matrix_set_col (varBeta, j, &v0j.vector);
            dev[j]=fit0.dev[j];
       }
       // other properties based on mu and phi
       for (i=0; i<nRows; i++) {
           yij = gsl_matrix_get(Y, i, j);
           mij = gsl_matrix_get(Mu, i, j);
           // get variance
           vij = varfunc( mij, phi[j] );
           gsl_matrix_set(Var, i, j, vij); 
           // get weight
           gsl_matrix_set(wHalf, i, j, sqrt(weifunc(mij, phi[j]))); 
           // get (Pearson) residuals
           gsl_matrix_set(Res, i, j, (yij-mij)/sqrt(vij));        
           // get elementry log-likelihood
           ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
       }
       aic[j]=-ll[j]+2*(nParams+1);

       // Get X * W^1/2
       wj = gsl_matrix_column (wHalf, j);
       for (i=0; i<nParams; i++) 
           gsl_matrix_set_col (WX, i, &wj.vector);
       gsl_matrix_mul_elements (WX, X);
       // X^T * W * X
       gsl_matrix_set_zero (XwX);
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
       gsl_linalg_cholesky_decomp (XwX); // provided XwX is non-singular        
       // Calc varBeta
       gsl_linalg_cholesky_invert (XwX);
       vj = gsl_matrix_column (varBeta, j);
       dj = gsl_matrix_diagonal (XwX);
       gsl_vector_memcpy (&vj.vector, &dj.vector);
//       displayvector(&vj.vector, "vj");

       // Calc sqrt(1-hii)
       // hii is diagonal element of W^1/2*X*(X'WX)^-1*X^T*W^1/2
       gsl_blas_dsymm (CblasRight, CblasLower, 1.0, XwX, WX, 0.0, TMP);
       gsl_matrix_mul_elements (TMP, WX);
       gsl_vector_set_all (t, 1.0);
       hj = gsl_matrix_column (sqrt1_Hii, j);
       gsl_blas_dgemv (CblasNoTrans, 1.0, TMP, t, 0.0, &hj.vector); 
       for (i=0; i<nRows; i++) {
           hii = gsl_vector_get(&hj.vector, i);
           gsl_vector_set(&hj.vector, i, sqrt(1-MIN((1-mintol),hii))); 
       }    
   } // end nVar for j loop
   gsl_matrix_div_elements (Res, sqrt1_Hii);
   subtractMean(Res);

   fit0.releaseGlm();  
   gsl_matrix_free(XwX);
   gsl_matrix_free(WX);
   gsl_matrix_free(TMP);
   gsl_vector_free(t);

   return SUCCESS;    
}


double PoissonGlm::getDisper( unsigned int id ) const
{
    unsigned int i, df, nNonZero=0;
    double ss2, yij, mij, chi2=0;

    gsl_vector_view yj = gsl_matrix_column (Yref, id);
    gsl_vector_view mj = gsl_matrix_column (Mu, id);
    for (i=0; i<nRows; i++) {
        yij = gsl_vector_get (&yj.vector, i);
        mij = gsl_vector_get (&mj.vector, i);
	ss2 = (yij-mij)*(yij-mij); // ss = (y-mu)^2
	if ( mij < mintol ) mij = 1;
	else  nNonZero++;	   
//	printf("ss2=%.2f, var=%.2f ", ss2, varfunc(mij, phi[id]));
        chi2 = chi2 + ss2/varfunc(mij, phi[id]); // dist dependant
    }
    if (nNonZero > nParams) 
        df = nNonZero - nParams; 
    else df = 1;
//    df = nRows - nParams;    
    return chi2/df;
}


double NBinGlm::llfunc ( double yi, double mui, double a  ) const
{
    double l=0, k=0, p=0;

    if (a < 0 ) 
       GSL_ERROR("Error in llfunc, phi should be non-negative", GSL_ERANGE); 
    else if (a<mintol) { // equivalent to Poisson log likelihood
       if (yi < 0)
          GSL_ERROR("Error in llfunc, y should be non-negative", GSL_ERANGE);
       else if (yi<mintol) l = -mui;
       else l = (yi*log(mui)-mui-gsl_sf_lngamma(yi+1));
    }
    else {
       if (yi < 0)
          GSL_ERROR("Error in llfunc, y should be non-negative", GSL_ERANGE);
       else if (yi<mintol)  l = -(log(1+mui*a))/a;    
       else {
          k = 1/a;
          l = gsl_sf_lngamma(yi+k)-gsl_sf_lngamma(k)-gsl_sf_lngamma(yi+1);
          p = 1/(1+GSL_MAX(mui, mintol)*a);
          l = l + log(p)*k + yi*log( ((1-p)<mintol) ? 1:(1-p) );  
       }
    }

    return 2*l;
}


int NBinGlm::getfAfAdash(double a, unsigned int id, double *fAPtr, double *fAdashPtr )
{
    unsigned int i;
    double yij, mij, dl, ddl, k=1/a;
    *fAPtr = 0;
    *fAdashPtr = 0;
    for ( i=0; i<nRows; i++ ) {
        yij = gsl_matrix_get(Yref, i, id);
	mij = gsl_matrix_get(Mu, i, id);
//	printf("%.1f ", mij);
        dl=gsl_sf_psi(yij+k)-gsl_sf_psi(k)-log(mij+k)+log(k)-(yij-mij)/(mij+k); // dl/da
	*fAPtr = *fAPtr + dl;  // sum
	ddl=gsl_sf_psi_1(yij+k)-gsl_sf_psi_1(k)+a*mij/(mij+k)+(yij-mij)/(mij+k)/(mij+k); // dl^2/d^2a
	*fAdashPtr = *fAdashPtr + ddl + 2*a*dl;	//sum
    }

    *fAPtr = - k*k*(*fAPtr);
    *fAdashPtr = exp(4*log(k))*(*fAdashPtr);

    return SUCCESS;
    
}

/*
void glm::display(void)
{   
    unsigned int j;
    if ( mmRef->model == LM )
       printf("Linear regression:\n");
    else if ( mmRef->model == POISSON )
       printf("Poisson regression:\n");
    else if ( mmRef->model == LOGIT )
       printf("Logistic regression:\n");
    else if ( mmRef->model == NB ) {
       printf("Negative Binomial regression ");	
       switch (mmRef->estiMethod) {
          case NEWTON:
              printf("(Newton-ML):\n");
	      break;
	  case CHI2:
	      printf("(Chi2):\n");
	      break;
          case FISHER:
              printf("(Fisher Scoring):\n");
	      break;
          default: 
              printf("no such method available.\n");
              exit(1);
       }
    }   
    else {
        printf("regression not available.\n");
	exit(1);
    }

//    printf("Two-log-like=\n " );
//    for ( j=0; j<nVars; j++ ) printf("%.2f ", ll[j]);	
//    printf("\n");
//    printf("AIC=\n " );
//    for ( j=0; j<nVars; j++ ) printf("%.2f ", aic[j]);	
//    printf("\n");
//    printf("# of convergence\n");    
//    for ( j=0; j<nVars; j++ )
//        printf("%d ", iterconv[j]); 
//    printf("\n");	       
//    printf("Residual deviance=\n " );
//    for ( j=0; j<nVars; j++ ) printf("%.2f ", dev[j]);
//    printf("\n");	       
    if ( mmRef->model == NB ) {
        printf("\nphi=\n ");
        for (j=0; j<nVars; j++ ) printf("%.2f ", phi[j]);
    }
    printf("\n");	       
//    if (Oref != NULL)
//       displaymatrix(Oref, "O");
//    displaymatrix(Xref, "X");
//    displaymatrix(Eta, "Eta");
//    displaymatrix(Beta, "Beta");
//    displaymatrix(varBeta, "varBeta");
//    displaymatrix(Mu, "Mu");
//    displaymatrix(Var, "Var");    
//    displaymatrix(Res, "Res");    
//    displaymatrix(sqrt1_Hii, "sqrt1_Hii");    
//    displaymatrix(wHalf, "wHalf");
    
}
*/
