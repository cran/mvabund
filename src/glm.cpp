// GLM estimation 
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 16-333-2011

#include "resampTest.h"
#include <string.h>
#include <gsl/gsl_sf_gamma.h>  
#include <gsl/gsl_sf_psi.h>

// Note try to use gsl functions as much as poosible to increase speed and stabilitiy

glm::glm(const reg_Method *mm)
   : mmRef(mm), Yref(NULL), Xref(NULL), Oref(NULL), 
     Beta(NULL), varBeta(NULL), Mu(NULL), Eta(NULL), Res(NULL), 
     Var(NULL), wHalf(NULL), sqrt1_Hii(NULL), PitRes(NULL), 
     phi(NULL), ll(NULL), dev(NULL), aic(NULL), iterconv(NULL) 
{ 
     n=mmRef->n;
     maxiter = 50;
     eps = 1e-5;
     if (mmRef->model==BIN) maxtol=n-eps;
     else maxtol = 1/eps;
}

PoissonGlm::PoissonGlm(const reg_Method *mm):glm(mm){
}

BinGlm::BinGlm(const reg_Method *mm):PoissonGlm(mm) {
}

NBinGlm::NBinGlm(const reg_Method *mm):PoissonGlm(mm){ 
}

glm::~glm() {
}


PoissonGlm::~PoissonGlm() {
}

BinGlm::~BinGlm() {
}

//LogiGlm::~LogiGlm() {
//}

NBinGlm::~NBinGlm() {
}

void glm::releaseGlm(void)
{ 
    if (Xref!=NULL) 
       gsl_matrix_free(Xref);
    if (Yref!=NULL) 
       gsl_matrix_free(Yref);
    if (Oref!=NULL) 
       gsl_matrix_free(Oref);
    if (Beta!=NULL)
       gsl_matrix_free(Beta);
    if (varBeta!=NULL)
        gsl_matrix_free(varBeta);
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
    if (PitRes!=NULL)
        gsl_matrix_free(PitRes);
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
}

void glm::initialGlm(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B) 
{ 
    releaseGlm();

    nRows = Y->size1;
    nVars = Y->size2;
    nParams = X->size2;

    unsigned int i, j;
    phi = new double [nVars];
    ll = new double [nVars];
    dev = new double [nVars];
    aic = new double [nVars];
    iterconv = new unsigned int [nVars]; 

    Xref = gsl_matrix_alloc(nRows, nParams);
    gsl_matrix_memcpy (Xref, X);
    Yref = gsl_matrix_alloc(nRows, nVars);
    gsl_matrix_memcpy (Yref, Y);
    if (O==NULL) Oref=NULL;
    else {
        Oref = gsl_matrix_alloc(nRows, nVars);
        gsl_matrix_memcpy(Oref, O);
    }
    Beta = gsl_matrix_alloc(nParams, nVars);
    varBeta = gsl_matrix_alloc(nParams, nVars);
    Mu = gsl_matrix_alloc(nRows, nVars);
    Eta = gsl_matrix_alloc(nRows, nVars);
    Res = gsl_matrix_alloc(nRows, nVars);
    Var = gsl_matrix_alloc(nRows, nVars);
    wHalf = gsl_matrix_alloc(nRows, nVars);
    sqrt1_Hii = gsl_matrix_alloc(nRows, nVars);
    PitRes = gsl_matrix_alloc(nRows, nVars);

    gsl_matrix_set_zero (varBeta);

    for (j=0; j<nVars; j++) {
        phi[j] = 0;
        ll[j] = 0;
        dev[j] = 0;
        aic[j] = 0;
        iterconv[j] = 0;
    }
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
    if (B!=NULL) {
       gsl_matrix_memcpy(Beta, B);
       gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0,X,Beta,0.0,Eta);
       for (i=0; i<nRows; i++)
       for (j=0; j<nVars; j++)
           gsl_matrix_set(Mu, i, j, invLink(gsl_matrix_get(Eta, i, j)));
    }  
    else if (O!=NULL) {
       gsl_matrix_memcpy(Eta, O);
       for (i=0; i<nRows; i++)
       for (j=0; j<nVars; j++)
           gsl_matrix_set(Mu, i, j, invLink(gsl_matrix_get(Eta, i, j)));
       gsl_matrix_set_zero(Beta);   
    }
    else {
       gsl_matrix_memcpy(Mu, Y);
       gsl_matrix_add_constant(Mu, 0.5);
       gsl_matrix_scale(Mu, 0.5);
       for (i=0; i<nRows; i++)
       for (j=0; j<nVars; j++)
          gsl_matrix_set(Eta, i, j, link(gsl_matrix_get(Mu, i, j)));
       gsl_matrix_set_zero (Beta); // intercept
       gsl_vector_view b0=gsl_matrix_column(Beta, 0);
       gsl_vector_set_all(&b0.vector, 1.0);
    }

//   displaymatrix(Eta, "Eta");
//   displaymatrix(Mu, "Mu");

   rdf = nRows - nParams;
}


int glm::copyGlm(glm *src)
{    
    initialGlm(src->Yref, src->Xref, src->Oref, NULL);

    // copy properties
    gsl_matrix_memcpy(Beta, src->Beta);
    gsl_matrix_memcpy(varBeta, src->varBeta);
    gsl_matrix_memcpy(Mu, src->Mu);
    gsl_matrix_memcpy(Eta, src->Eta);
    gsl_matrix_memcpy(Res, src->Res);
    gsl_matrix_memcpy(Var, src->Var);
    gsl_matrix_memcpy(wHalf, src->wHalf);
    gsl_matrix_memcpy(sqrt1_Hii, src->sqrt1_Hii);
    gsl_matrix_memcpy(PitRes, src->PitRes);
    
    for (unsigned int i=0; i<nVars; i++) {
        phi[i] = src->phi[i];	
        ll[i] = src->ll[i];
        dev[i] = src->dev[i];
        iterconv[i] = src->iterconv[i];
	aic[i] = src->aic[i];
    }
    
    return SUCCESS;    
}

int PoissonGlm::EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B, double *a)
{
    initialGlm(Y, X, O, B);

    gsl_set_error_handler_off();
    gsl_rng *rnd=gsl_rng_alloc(gsl_rng_mt19937);
    unsigned int i, j;   
    int status;
    double yij, mij, vij, wij, tol, hii, uij, wei;
    gsl_vector_view Xwi, Xi, vj, hj, dj;

    gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *TMP = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);   

    for (j=0; j<nVars; j++) {
       if ( a!=NULL ) phi[j]=a[j]; 
       // estimate mu and beta   
       iterconv[j] = betaEst(j, maxiter, &tol, phi[j]); 
       if (iterconv[j]>maxiter) 
           printf("EstIRLS reached max iterations\n");
       gsl_matrix_memcpy (WX, X);
       for (i=0; i<nRows; i++) {
            mij = gsl_matrix_get(Mu, i, j);
            // get variance
            vij = varfunc( mij, phi[j] );
            gsl_matrix_set(Var, i, j, vij); 
            // get weight
            wij = sqrt(weifunc(mij, phi[j]));  
            gsl_matrix_set(wHalf, i, j, wij);             
            // get (Pearson) residuals
            yij = gsl_matrix_get(Y, i, j);
            gsl_matrix_set(Res, i, j, (yij-mij)/sqrt(vij));        
            // get PIT residuals for discrete data
            wei = gsl_rng_uniform_pos (rnd); // wei ~ U(0, 1)
            uij = wei*cdf(yij, mij, phi[j]);
            if (yij>0) uij=uij+(1-wei)*cdf((yij-1),mij,phi[j]);
            gsl_matrix_set(PitRes, i, j, uij);
            // get elementry log-likelihood	   
            ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
            // W^1/2 X
            Xwi = gsl_matrix_row (WX, i);
            gsl_vector_scale(&Xwi.vector, wij);            
       }      
       aic[j]=-ll[j]+2*(nParams);

       // X^T * W * X
       gsl_matrix_set_zero(XwX);
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
       if (calcDet(XwX)<eps) {  // XwX+eps*I
          dj=gsl_matrix_diagonal(XwX);
          gsl_vector_add_constant(&dj.vector, eps);
       } 
       status=gsl_linalg_cholesky_decomp (XwX);
       if (status) {
          printf("Singular info matrix in EstIRLS\n");
       }
       gsl_linalg_cholesky_invert (XwX);

       // Calc varBeta
       dj = gsl_matrix_diagonal (XwX);
       vj = gsl_matrix_column (varBeta, j);       
       gsl_vector_memcpy (&vj.vector, &dj.vector);

       // hii is diagonal element of H=X*(X'WX)^-1*X'*W
       hj = gsl_matrix_column (sqrt1_Hii, j);
       gsl_blas_dsymm(CblasRight,CblasLower,1.0,XwX,Xref,0.0,TMP); // X*(X'WX)^-1
       for (i=0; i<nRows; i++) {
           Xwi=gsl_matrix_row(TMP, i);
           Xi=gsl_matrix_row(Xref, i);
           wij=gsl_matrix_get(wHalf, i, j);
           gsl_blas_ddot(&Xwi.vector, &Xi.vector, &hii);
           gsl_vector_set(&hj.vector, i, GSL_MAX(eps, sqrt(1-wij*wij*hii)));
       } 
   } 
   // standardize perason residuals by rp/sqrt(1-hii) 
   gsl_matrix_div_elements (Res, sqrt1_Hii);
   subtractMean(Res);  // have mean subtracted

   gsl_matrix_free(XwX);
   gsl_matrix_free(WX);
   gsl_matrix_free(TMP);
   gsl_rng_free(rnd);

   return SUCCESS;    
}


int PoissonGlm::betaEst( unsigned int id, unsigned int iter, double *tol, double phi )
{
   gsl_set_error_handler_off();
   int status, isValid;
   unsigned int i, j, ngoodobs;
   unsigned int step, step1; 
   double wij, zij, eij, mij, yij;   
   double dev_old, dev_grad=1.0;
   unsigned int *good = new unsigned int [nRows];
   gsl_vector_view Xwi;
   gsl_matrix *WX, *XwX;
   gsl_vector *z, *Xwz;
   gsl_vector *coef_old = gsl_vector_alloc(nParams);
   gsl_vector_view bj=gsl_matrix_column (Beta, id);

   // Main Loop of IRLS begins
   XwX = gsl_matrix_alloc(nParams, nParams);
   Xwz = gsl_vector_alloc(nParams);
   step=0;
   *tol = 1.0;
   gsl_vector_memcpy (coef_old, &bj.vector);
   while ( step<iter ) {
       // Decide good observations
       ngoodobs=0;
       for (i=0; i<nRows; i++) {
           good[i]=TRUE;
           eij = gsl_matrix_get(Eta, i, id);
           mij = gsl_matrix_get(Mu, i, id);  
           if ((mij<eps)|(mij>maxtol)) 
	      good[i]=FALSE;
  //         if ((mmRef->model==BIN)&(rcpLinkDash(mij)<eps*maxtol))
  //            good[i]=FALSE;
  //         if ((mmRef->model!=BIN)&(rcpLinkDash(mij)<eps))
           if (rcpLinkDash(mij)<eps) good[i]=FALSE;
           if (good[i]==TRUE) ngoodobs++;
       }  
    
       if (ngoodobs==0) {
 //          printf("mij(rcpLinkDash)\n");
//           for (i=0; i<nRows; i++) {
//              mij = gsl_matrix_get(Mu, i, id);
 //             printf("%.8f(%.8f)\t", mij, rcpLinkDash(mij));
 //          }
 //          printf("\nError: ngoodobs=0\n");      
          for (i=0; i<nRows; i++) {
              good[i]=TRUE;
	      if (rcpLinkDash(mij)<eps*0.1) good[i]=FALSE;
	      if (good[i]==TRUE) ngoodobs++;
	  }
       }

//       gsl_vector_set_zero (&wj.vector);
       z = gsl_vector_alloc(ngoodobs);
       WX = gsl_matrix_alloc(ngoodobs, nParams); 
       j=0;  // index for z and WX
       for (i=0; i<nRows; i++) { // (y-m)/g'
           if ( good[i]==TRUE ) {
              // z=xb+(y-m)/g'-o
              yij = gsl_matrix_get(Yref, i, id);
              mij = gsl_matrix_get(Mu, i, id);
              eij = gsl_matrix_get(Eta, i, id);
              zij = eij + (yij-mij)/rcpLinkDash(mij);
              if (Oref!=NULL) 
                 zij = zij - gsl_matrix_get(Oref, i, id);
              // wt=sqrt(weifunc);
              wij = sqrt(weifunc(mij, phi));
              // W^1/2*z[good]
              gsl_vector_set(z, j, wij*zij); 
              // W^1/2*X[good]
              Xwi = gsl_matrix_row (Xref, i);
              gsl_matrix_set_row (WX, j, &Xwi.vector);
              Xwi = gsl_matrix_row (WX, j);
              gsl_vector_scale(&Xwi.vector, wij); 
              j++;
          }
       }
       // in glm2, solve WXb=Wz, David suggested not good 
       // So back to solving X'WXb=X'Wz
       gsl_matrix_set_zero (XwX);
       gsl_blas_dsyrk (CblasLower,CblasTrans,1.0,WX,0.0,XwX); 
       if (calcDet(XwX)<eps){ // test if singular or nan
          gsl_vector_view dj=gsl_matrix_diagonal(XwX);
          gsl_vector_add_constant(&dj.vector, eps);
       }   
       status=gsl_linalg_cholesky_decomp(XwX);
       if (status) {
          printf("Singular XwX in betaEst\n");
       }
       gsl_blas_dgemv(CblasTrans,1.0,WX,z,0.0,Xwz);
       gsl_linalg_cholesky_solve (XwX, Xwz, &bj.vector);

       // Given bj, update eta, mu
       dev_old = dev[id];
       isValid=predict(bj, coef_old, id, phi);
       dev_grad=(dev[id]-dev_old)/(ABS(dev[id])+0.1);
       *(tol)=ABS(dev_grad);  

       step1 = 0;
       // If divergent or increasing deviance, half step
       // (step>1) -> (step>0) gives weired results for NBin fit       
       // below works for boundary values, esp BIN fit but not NBin fit
       while (((dev[id]>10)|(dev_grad>-0.25))&(step>1)){
            gsl_vector_add (&bj.vector, coef_old);
            gsl_vector_scale (&bj.vector, 0.5);
            dev_old=dev[id];
            isValid=predict(bj, coef_old, id, phi);
            dev_grad=(dev[id]-dev_old)/(ABS(dev[id])+0.1); 
            *tol=ABS(dev_grad);
            if (*tol<eps) break;
            step1++;
            if (step1>maxiter) {
               printf("\t Internal loop reached max iter, d_dev=%.4f\n", dev_grad);
               break;
            }
       }

       gsl_vector_free(z);
       gsl_matrix_free(WX); 
       if ((isValid==TRUE))
          gsl_vector_memcpy (coef_old, &bj.vector);
       step++;
       if (*tol<eps) break;

   } 

//   if ((step>1)&(isValid==FALSE)) { // restore previous valid one
//      update(coef_old, id);
//      dev[id]=0;
//      for ( i=0; i<nRows; i++ ) {
//          yij = gsl_matrix_get(Yref, i, id);
//          mij = gsl_matrix_get(Mu, i, id);         
//          dev[id] = dev[id] + devfunc(yij, mij, phi);
//      }
//   }

   delete[] good;
   gsl_matrix_free(XwX); 
   gsl_vector_free(Xwz); 
   gsl_vector_free(coef_old); 

   return step;
}

int PoissonGlm::update(gsl_vector *bj, unsigned int id)
{
    int isValid=TRUE;
    unsigned int i;
    double eij, mij;
    gsl_vector_view xi;
 
    for (i=0; i<nRows; i++) {
       xi = gsl_matrix_row (Xref, i);
       gsl_blas_ddot (&xi.vector, bj, &eij);
       if (Oref!=NULL) 
          eij = eij+gsl_matrix_get(Oref, i, id);
       mij = invLink(eij);
       if ((mij<eps)|(mij>maxtol)){
          mij=(mij<eps)?eps:((mij>maxtol)?maxtol:mij);
          eij=link(mij);
          isValid = FALSE;
       }
       gsl_matrix_set(Eta, i, id, eij);
       gsl_matrix_set(Mu, i, id, mij);
   } 

   return isValid;
}

int PoissonGlm::predict(gsl_vector_view bj, gsl_vector *coef_old, unsigned int id, double phi) 
{
    unsigned int i;
    double yij, mij;
    
    int isValid=update(&bj.vector, id);

/*  unsigned int step=0;
    double bij, diff, delta;
    while ( isValid==FALSE ) {
         // Step too large, half step 
         gsl_vector_add (&bj.vector, coef_old);
         gsl_vector_scale (&bj.vector, 0.5);
         isValid = update(&bj.vector, id);
         step++;
         // If difference small, stop
         delta=0;
         for (i=0; i<nParams; i++) {
             bij = gsl_vector_get(&bj.vector, i);
             diff=bij-gsl_vector_get(coef_old, i);
             delta=delta+ABS(diff/(ABS(bij)+0.1));
         }
         if (delta<eps) { // accuracy 
//            printf("\t\t Internal loop II(%d), delta=%.4f\n", step, delta);
             break;
         }
         if (step>maxiter) {
             printf("\t\t Internal loop II(%d) reached maximum iter %d, delta=%.4f\n", isValid, step, delta);
             break;
         }
     }
*/
     // Given valid mj, estimate deviance
     dev[id]=0;
     for ( i=0; i<nRows; i++ ) {
         yij = gsl_matrix_get(Yref, i, id);
         mij = gsl_matrix_get(Mu, i, id);         
         dev[id] = dev[id] + devfunc(yij, mij, phi);
     }

     return isValid;
}


int NBinGlm::nbinfit(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, gsl_matrix *B)
{   
    gsl_set_error_handler_off();

    initialGlm(Y, X, O, B);

    gsl_rng *rnd=gsl_rng_alloc(gsl_rng_mt19937);
    unsigned int i, j, isConv;
    double yij, mij, vij, hii, uij, wij, wei;
    double a, phi_old, fA, fAdash, tol;
//    double tol_grad, tol_old;
    int status;
    gsl_vector_view b0j, m0j, e0j, v0j;
    gsl_matrix *WX = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *TMP = gsl_matrix_alloc(nRows, nParams);   
    gsl_matrix *XwX = gsl_matrix_alloc(nParams, nParams);   
    gsl_vector_view Xwi, Xi, vj, dj, hj;

    // Get initial estimates from Poisson    
    PoissonGlm fit0( mmRef );
    fit0.initialGlm( Y, X, O, B);    

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
        iterconv[j]=0;
        if (mmRef->estiMethod==CHI2) {
            a = fit0.getDisper(j); 
            phi[j] = (a<1)? 0:1/a;
	    isConv = (a<1)? TRUE:FALSE;
	    while (isConv != TRUE) {
                iterconv[j]++;
	        betaEst(j, 1, &tol, phi[j]);  // 1-step beta
		phi[j] = phi[j]*getDisper(j); // 1-step phi
	        if (tol<eps) break; // Normal break;
                if ( (phi[j]<0) | (phi[j]!=phi[j])) break; // invalide phi 
                if (iterconv[j]>maxiter) {
                   printf("NBinFit(Chi2) reached maximum iter\n"); 
                   break;
                }
            }
        }   
        else if (mmRef->estiMethod==NEWTON) {         
            getfAfAdash(eps, j, &fA, &fAdash);
            a = (fA>0)?-fA/fAdash:0;
            phi[j] = (a<0)? eps:a; 
	    isConv = (fA>0)? FALSE:TRUE;
            phi_old = phi[j];
	    while ( isConv != TRUE ) {
                iterconv[j]++;	    
	        betaEst(j, 1, &tol, phi[j]); // 1-step beta (better than m-step beta)
                getfAfAdash(phi[j], j, &fA, &fAdash); // 1-step phi
                if ( (ABS(fAdash)<eps)|(ABS(fA/fAdash)>0.5) ) 
                   phi[j] = (phi[j]+phi_old)/2; // half step
                else {
                   phi_old = phi[j];
                   phi[j] = phi[j]-fA/fAdash;
                }        
	        if ( tol<eps ) break; // Normal break;
                if ( (phi[j]<0) | (phi[j]!=phi[j]) ) break; // invalid phi
                if (iterconv[j]>maxiter) {                
//                   printf("NBinFit(ML) reached %d, phi[j]=%.4f, fA=%.4f, (.)=%.4f, tol=%.6f, tol_old=%.6f, tol_grad=%.6f\n", iterconv[j], phi[j], fA, (fA/fAdash), tol, tol_old, tol_grad); 
                   break; 
                } 
           }
       }  

       // restore poisson if phi[j]<0 or nan
       if ( (phi[j] < 0) | (phi[j]!=phi[j]) ) { 
            phi[j]=0;
            gsl_matrix_set_col (Beta, j, &b0j.vector);
            gsl_matrix_set_col (Mu, j, &m0j.vector);
            gsl_matrix_set_col (Eta, j, &e0j.vector);
            gsl_matrix_set_col (varBeta, j, &v0j.vector);
            dev[j]=fit0.dev[j];
       }
       // other properties based on mu and phi
       gsl_matrix_memcpy(WX, Xref);
       for (i=0; i<nRows; i++) {
           yij = gsl_matrix_get(Y, i, j);
           mij = gsl_matrix_get(Mu, i, j);
           vij = varfunc( mij, phi[j] );
           gsl_matrix_set(Var, i, j, vij); 
           wij = sqrt(weifunc(mij, phi[j]));
           gsl_matrix_set(wHalf, i, j, wij); 
           gsl_matrix_set(Res, i, j, (yij-mij)/sqrt(vij));        
           ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
           // get PIT residuals for discrete data
           wei = gsl_rng_uniform_pos (rnd); // wei ~ U(0, 1)
           uij = wei*cdf(yij, mij, phi[j]);
           if (yij>0) uij=uij+(1-wei)*cdf((yij-1),mij,phi[j]);           
           gsl_matrix_set(PitRes, i, j, uij);
           // W^1/2 X
           Xwi = gsl_matrix_row (WX, i);
           gsl_vector_scale(&Xwi.vector, wij);
       }
       aic[j]=-ll[j]+2*(nParams+1);

       // X^T * W * X
       gsl_matrix_set_zero (XwX);
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, WX, 0.0, XwX);
       if (calcDet(XwX)<eps) { //XwX + eps*I
          dj=gsl_matrix_diagonal(XwX);
          gsl_vector_add_constant(&dj.vector, eps);
       }       
       status=gsl_linalg_cholesky_decomp (XwX);
       if (status) {
          displaymatrix(XwX, "XwX");
          printf("Singular info mat in NBinFit - calcDet(XwX)=%.4f\n", calcDet(XwX)); 
       }
       gsl_linalg_cholesky_invert (XwX); // (X'WX)^-1

       // Calc varBeta
       vj = gsl_matrix_column (varBeta, j);
       dj = gsl_matrix_diagonal (XwX);
       gsl_vector_memcpy (&vj.vector, &dj.vector);

       // hii is diagonal element of H=X*(X'WX)^-1*X'*W
       hj = gsl_matrix_column (sqrt1_Hii, j);
       gsl_blas_dsymm(CblasRight,CblasLower,1.0,XwX,Xref,0.0,TMP); // X*(X'WX)^-1
       for (i=0; i<nRows; i++) {
           Xwi=gsl_matrix_row(TMP, i);
           Xi=gsl_matrix_row(Xref, i);
           wij=gsl_matrix_get(wHalf, i, j);
           gsl_blas_ddot(&Xwi.vector, &Xi.vector, &hii);
           gsl_vector_set(&hj.vector, i, GSL_MAX(eps,sqrt(1-wij*wij*hii)));
       }
   } // end nVar for j loop
   gsl_matrix_div_elements (Res, sqrt1_Hii);
   subtractMean(Res);

   fit0.releaseGlm();  
   gsl_matrix_free(XwX);
   gsl_matrix_free(WX);
   gsl_matrix_free(TMP);
   gsl_rng_free(rnd);

   return SUCCESS;    
}


double PoissonGlm::getDisper( unsigned int id )const
{
    unsigned int i, df, nNonZero=0;
    double ss2, yij, mij, chi2=0;

    gsl_vector_view yj = gsl_matrix_column (Yref, id);
    gsl_vector_view mj = gsl_matrix_column (Mu, id);
    for (i=0; i<nRows; i++) {
        yij = gsl_vector_get (&yj.vector, i);
        mij = gsl_vector_get (&mj.vector, i);
	ss2 = (yij-mij)*(yij-mij); // ss = (y-mu)^2
	if ( mij < eps ) mij = 1;
	else  nNonZero++;	   
        if ( varfunc(mij, phi[id])>eps )
            chi2 = chi2 + ss2/varfunc(mij, phi[id]); // dist dependant
    }
    if (nNonZero > nParams) 
        df = nNonZero - nParams; 
    else df = 1;
//    df = nRows - nParams;    
    return chi2/df;
}


double NBinGlm::llfunc ( double yi, double mui, double a  )const 
{
    double l=0, k=0, p=0;
//    double eps=mmRef->tol;

    if (a < 0 ) 
       GSL_ERROR("Error in llfunc, phi should be non-negative", GSL_ERANGE); 
    else if (a<eps) { // equivalent to Poisson log likelihood
       if (yi < 0)
          GSL_ERROR("Error in llfunc, y should be non-negative", GSL_ERANGE);
       else if (yi<eps) l = -mui;
       else l = (yi*log(mui)-mui-gsl_sf_lngamma(yi+1));
    }
    else {
       if (yi < 0)
          GSL_ERROR("Error in llfunc, y should be non-negative", GSL_ERANGE);
       else if (yi<eps)  l = -(log(1+mui*a))/a;    
       else {
          k = 1/a;
          l = gsl_sf_lngamma(yi+k)-gsl_sf_lngamma(k)-gsl_sf_lngamma(yi+1);
          p = 1/(1+GSL_MAX(mui, eps)*a);
          l = l + log(p)*k + yi*log( ((1-p)<eps) ? 1:(1-p) );  
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
        dl=gsl_sf_psi(yij+k)-gsl_sf_psi(k)-log(mij+k)+log(k)-(yij-mij)/(mij+k); // dl/da
	*fAPtr = *fAPtr + dl;  // sum
	ddl=gsl_sf_psi_1(yij+k)-gsl_sf_psi_1(k)+a*mij/(mij+k)+(yij-mij)/(mij+k)/(mij+k); // dl^2/d^2a
	*fAdashPtr = *fAdashPtr + ddl + 2*a*dl;	//sum
    }

    *fAPtr = - k*k*(*fAPtr);
    *fAdashPtr = exp(4*log(k))*(*fAdashPtr);

    return SUCCESS;
    
}

void glm::display(void)
{   
    unsigned int j;
/*
    if ( mmRef->model == LM )
       printf("Linear regression:\n");
    else if ( mmRef->model == POISSON )
       printf("Poisson regression:\n");
    else if ( mmRef->model == BIN ) {
       if (n==1 ) printf("Logistic regression:\n");
       else printf("Binomial regression:\n");
    }
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
              printf("phi estimation method not available");
       }
    }   
    else { 
        printf("GLM regression method not available");
    }
*/
    printf("Two-log-like=\n " );
    for ( j=0; j<nVars; j++ ) printf("%.2f ", ll[j]);	
    printf("\n");
//    printf("AIC=\n " );
//    for ( j=0; j<nVars; j++ ) printf("%.2f ", aic[j]);	
//    printf("\n");
//    printf("# of convergence\n");    
//    for ( j=0; j<nVars; j++ )
//        printf("%d ", iterconv[j]); 
//    printf("\n");	       
//   printf("Residual deviance=\n " );
//    for ( j=0; j<nVars; j++ ) printf("%.2f ", dev[j]);
//    printf("\n");	       
//    if ( mmRef->model == NB ) {
//        printf("\nphi=\n ");
//        for (j=0; j<nVars; j++ ) printf("%.2f ", phi[j]);
//    }
//    printf("\n");	       
//    if (Oref != NULL)
//       displaymatrix(Oref, "O");
//    displaymatrix(Xref, "X");
//    displaymatrix(Eta, "Eta");
//    displaymatrix(Beta, "Beta");
//    displaymatrix(varBeta, "varBeta");
//    displaymatrix(Mu, "Mu");
//    displaymatrix(Var, "Var");    
//    displaymatrix(PitRes, "PitRes");    
//    displaymatrix(sqrt1_Hii, "sqrt1_Hii");    
//    displaymatrix(wHalf, "wHalf");
    
}
