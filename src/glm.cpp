// GLM estimation 
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// // 16-June-2011

#include "resampTest.h"
//#include "time.h"
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_sf_gamma.h>  
#include <gsl/gsl_sf_psi.h>

// Note try to use gsl functions as much as poosible to increase speed and stabilitiy

glm::glm(const reg_Method *mm)
   : mmRef(mm), Yref(NULL), Xref(NULL), Oref(NULL), 
     Beta(NULL), Mu(NULL), Eta(NULL), Res(NULL), 
     Var(NULL), wHalf(NULL), sqrt1_Hii(NULL), phi(NULL),
     ll(NULL), dev(NULL), aic(NULL), iterconv(NULL) 
{ 
     //    printf("glm constructor called.\n"); 
}

PoissonGlm::PoissonGlm(const reg_Method *mm):glm(mm)
{
    maxiter = 50;
    mintol = 1e-4;    
    lTol=-log(mintol);
//    printf("Poisson constructor called.\n");
}

LogiGlm::LogiGlm(const reg_Method *mm):PoissonGlm(mm) {
    maxiter = 1000;    
//    maxiter = 100;    
    mintol = 1e-10;
//    mintol = 1e-4;
    lTol=-log(mintol);
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

    //Xref = gsl_matrix_alloc(nRows, nParams);
    Beta = gsl_matrix_alloc(nParams, nVars);
    Mu = gsl_matrix_alloc(nRows, nVars);
    Eta = gsl_matrix_alloc(nRows, nVars);
    Res = gsl_matrix_alloc(nRows, nVars);
    Var = gsl_matrix_alloc(nRows, nVars);
    wHalf = gsl_matrix_alloc(nRows, nVars);
    sqrt1_Hii = gsl_matrix_alloc(nRows, nVars);

    //gsl_matrix_memcpy(Xref, X);
    gsl_matrix_set_zero (Beta);
    gsl_matrix *t1;
    t1 = gsl_matrix_alloc(nRows, 1);
    gsl_matrix_set_all (t1, 1.0); // intercept
    GetMean(t1, Y, Mu);
    gsl_matrix_free(t1);

    size_t i, j;
    double eij;
   // double lTol=-log(mintol);
    for (i=0; i<nRows; i++) 
    for (j=0; j<nVars; j++) {
       eij = link(gsl_matrix_get(Mu, i, j));
     //  eij=(eij<-lTol)?-lTol:((eij>lTol)?lTol:eij);
       gsl_matrix_set(Eta, i, j, eij);
    }

    phi = new double [nVars];
    ll = new double [nVars];
    dev = new double [nVars];
    aic = new double [nVars];
    iterconv = new int [nVars];    
   
    rdf = nRows - nParams;

}

void glm::InitResampGlm(gsl_matrix *T, gsl_matrix *X, gsl_matrix *O, glm *src, gsl_matrix *Hii)
{
    releaseGlm();

    Yref = NULL;
    Xref = X;
    Oref = O;

    nRows = X->size1;
    nParams = X->size2;
    nVars = T->size2;

    //Xref = gsl_matrix_alloc(nRows, nParams);
    Beta = gsl_matrix_alloc(nParams, nVars);
    Mu = gsl_matrix_alloc(nRows, nVars);
    Eta = gsl_matrix_alloc(nRows, nVars);
    Res = gsl_matrix_alloc(nRows, nVars);
    Var = gsl_matrix_alloc(nRows, nVars);
    wHalf = gsl_matrix_alloc(nRows, nVars);
    sqrt1_Hii = gsl_matrix_alloc(nRows, nVars);

    phi = new double [nVars];
    ll = new double [nVars];
    dev = new double [nVars];
    aic = new double[nVars];
    iterconv = new int [nVars]; 
   
    rdf = nRows - nParams;

    //gsl_matrix_memcpy(Xref, X);
    gsl_matrix_set_zero (Beta);
    gsl_matrix_memcpy(Mu, src->Mu);
    gsl_matrix_memcpy(Eta, src->Eta);
    gsl_matrix_memcpy(Var, src->Var);
    gsl_matrix_memcpy(wHalf, src->wHalf);
    gsl_matrix_memcpy(Res, T);
    if ( Hii!=NULL )
       gsl_matrix_memcpy(sqrt1_Hii, Hii);
    
}

int glm::copyGlm(glm *src)
{    
//    releaseGlm();

    initialGlm(src->Yref, src->Xref, src->Oref);

    // copy properties
    Xref = gsl_matrix_alloc(src->nRows, src->nParams);
//    gsl_matrix_set_zero(Xref);
//    displaymatrix(Xref, "Xref");
//    displaymatrix(src->Xref, "src->Xref");
    gsl_matrix_memcpy(Xref, src->Xref);
    gsl_matrix_memcpy(Beta, src->Beta);
    gsl_matrix_memcpy(Mu, src->Mu);
    gsl_matrix_memcpy(Eta, src->Eta);
    gsl_matrix_memcpy(Res, src->Res);
    gsl_matrix_memcpy(Var, src->Var);
    gsl_matrix_memcpy(wHalf, src->wHalf);
    gsl_matrix_memcpy(sqrt1_Hii, src->sqrt1_Hii);
    
    for (size_t i=0; i<nVars; i++) {
        phi[i] = src->phi[i];	
        ll[i] = src->ll[i];
        dev[i] = src->dev[i];
        iterconv[i] = src->iterconv[i];
	aic[i] = src->aic[i];
    }
    
    return SUCCESS;    
}


void glm::display(void)
{   
    size_t j;
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
//    for ( j=0; j<nVars; j++ ) printf("%.3f ", ll[j]);	
//    printf("\n");
//    printf("AIC=\n " );
//    for ( j=0; j<nVars; j++ ) printf("%.3f ", aic[j]);	
//    printf("\n");
//    printf("Residual deviance=\n " );
//    for ( j=0; j<nVars; j++ )
//        printf("%d ", iterconv[j]); 
//        printf("%.2f ", dev[j]);
//    printf("\n");	       
    if ( mmRef->model == NB ) {
        printf("\nphi=\n ");
        for (j=0; j<nVars; j++ ) printf("%.4f ", phi[j]);
    }
    printf("\n");	       
//    if (Oref != NULL)
//       displaymatrix(Oref, "O");
    displaymatrix(Yref, "Y");
    displaymatrix(Eta, "Eta");
    displaymatrix(Beta, "Beta");
    displaymatrix(Mu, "Mu");
//    displaymatrix(Var, "Var");    
//    displaymatrix(Res, "Res");    
//    displaymatrix(wHalf, "wHalf");
    
}


int PoissonGlm::EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, double *a)
{
    initialGlm(Y, X, O);
//    displaymatrix(Mu, "Mu0");
//    displaymatrix(Beta, "Beta0");

    size_t i, j;   
    double yij, mij, vij, wij, rij, tol;

//    clock_t clk_start = clock();
    for (j=0; j<nVars; j++) {
       ll[j] = 0;
       dev[j] = 0;
       iterconv[j] = 1;
       if ( a == NULL ) phi[j]=0; 
       else phi[j] = a[j];

       // estimate mu and beta   
       iterconv[j] = betaEst(j, maxiter, &tol, phi[j]);        

       // other properties based on mu
       for (i=0; i<nRows; i++) {
            mij = gsl_matrix_get(Mu, i, j);
            // get variance
            vij = varfunc( mij, phi[j] );
            gsl_matrix_set(Var, i, j, vij); 
            // get weight
            wij = weifunc(mij, phi[j]);           
            gsl_matrix_set(wHalf, i, j, sqrt(wij)); 
            // get (Pearson) residuals
            yij = gsl_matrix_get(Y, i, j);
            rij = (yij-mij)/sqrt(vij);
            gsl_matrix_set(Res, i, j, rij);        
            // get elementry log-likelihood	   
            ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
            dev[j] = dev[j] + devfunc( yij, mij, phi[j] );
       }      
       aic[j]=-ll[j]+2*(nParams);
   } 
//   float time_diff = (float)(clock() - clk_start)/(float)(CLOCKS_PER_SEC);
//   printf("\n [var%d] betaEst time: %.6f sec", j, time_diff);
   // standardize perason residuals by rp/sqrt(1-hii) 
   // hii is diagonal element in Hat = X*(X'WX)^-1*X'W
   getHat(Xref, wHalf, sqrt1_Hii);
   gsl_matrix_div_elements (Res, sqrt1_Hii);
   subtractMean(Res);  // have mean subtracted
//   displaymatrix(Res, "Res");

   return SUCCESS;    
}


int PoissonGlm::betaEst( size_t id, size_t iter, double *tol, double a )
{
   size_t i, isConv=FALSE, step=1;
   double eij, mij=0, yij, oij, wij, zij;

   gsl_vector *z = gsl_vector_alloc(nRows);
   gsl_vector *bj_old = gsl_vector_alloc(nParams); 
   gsl_matrix *Xw = gsl_matrix_alloc(nRows, nParams);   

   gsl_vector_view bj = gsl_matrix_column (Beta, id);   
   gsl_vector_view ej=gsl_matrix_column(Eta, id);
   // IRLS
   while ( isConv != TRUE ) {
       gsl_vector_memcpy( bj_old, &bj.vector );      
       gsl_matrix_memcpy( Xw, Xref );
       for (i=0; i<nRows; i++) {
           eij = gsl_vector_get(&ej.vector, i);
           eij=(eij<-lTol)?-lTol:((eij>lTol)?lTol:eij);
           mij  = invLink(eij);
           // update weight 
	   wij = sqrt(weifunc(mij, a));
	   // update z
           yij = gsl_matrix_get(Yref, i, id);
	   if (Oref==NULL) oij = 0;
	   else oij = gsl_matrix_get(Oref, i, id); 
	   zij = eij + (yij-mij)/rcpLinkDash(mij) - oij;
           gsl_vector_set( z, i, wij*zij ); // z = wHalf * z
           // get Xw
	   gsl_vector_view Xwi = gsl_matrix_row (Xw, i);
           gsl_vector_scale (&Xwi.vector, wij); // Xw = WHalf * X
        }
//	displaymatrix(Xw, "Xw");
	invLSQ(Xw, z, &bj.vector);  // least square inverse
	// update eta = X*beta + offset
        gsl_blas_dgemv (CblasNoTrans, 1.0, Xref, &bj.vector, 0.0, &ej.vector);
	if (Oref!=NULL) {
	   gsl_vector_view oj=gsl_matrix_column(Oref, id);
           gsl_vector_add (&ej.vector, &oj.vector);	   
        }	
//	printf("step=%d ", step);
//	displayvector(ej, "ej");

	// Test convergence
        *tol = 0;
	gsl_vector_sub (bj_old, &bj.vector);
        for (i=0; i<nParams; i++) 
	    *tol = *tol + ABS(gsl_vector_get(bj_old, i));
        if ( (*tol < mintol) | (step == iter )) break;
	step++;
   } 
   // update mu
   for (i=0; i<nRows; i++) {
        eij=gsl_matrix_get (Eta, i, id);
        eij=(eij<-lTol)?-lTol:((eij>lTol)?lTol:eij);
        gsl_matrix_set(Eta, i, id, eij);
        gsl_matrix_set(Mu, i, id, invLink(eij));
   }
//   gsl_vector_view mj=gsl_matrix_column(Mu, id);
//   displayvector(&mj.vector, "mj");

   gsl_vector_free(z);
   gsl_vector_free(bj_old);
   gsl_matrix_free(Xw); 

   return step;
}

double PoissonGlm::getDisper( size_t id ) const
{
    size_t i,nNonZero=0;
    double ss2, yij, mij, chi2=0;

    gsl_vector_view yj = gsl_matrix_column (Yref, id);
    gsl_vector_view mj = gsl_matrix_column (Mu, id);
    for (i=0; i<nRows; i++) {
        yij = gsl_vector_get (&yj.vector, i);
        mij = gsl_vector_get (&mj.vector, i);
	ss2 = (yij-mij)*(yij-mij); // ss = (y-mu)^2
	if ( mij < mintol ) mij = 1;
	else  nNonZero++;	   
//	printf("ss2=%.3f ", ss2);
        chi2 = chi2 + ss2/varfunc(mij, phi[id]); // dist dependant
    }
 //   printf("chi2=%.3f\n", chi2);
    return chi2/(nNonZero-nParams);
}


double NBinGlm::llfunc ( double yi, double mui, double a  ) const
{
    double l=0, k=1/a, p;
    if ((a==0) & (yi==0)) l = -mui;
    else if ((a==0) & (yi>0))  l = (yi*log(mui)-mui-gsl_sf_lngamma(yi+1));  
    else if ((a>0) & (yi==0))  l = k*(log(k)-log(mui+k));    
    else if ((a>0) & (yi>0)) {
       p = 1/(1+MAX(mui, 0)*a);
//       double pComp = ((1-p)<mintol) ? 1:(1-p);
       l = gsl_sf_lngamma(yi+k)-gsl_sf_lngamma(k)-gsl_sf_lngamma(yi+1);
       l = l + log(p)*k + yi*log( ((1-p)==0) ? 1:(1-p) );  
    }
    else 
       GSL_ERROR("Error in llfunc, y or phi should be non-negative", GSL_ERANGE); 
/*    
    if ( a == 0 )
        l = (yi*log(MAX(mui, 1))-mui-gsl_sf_lngamma(MAX(yi, 0)+1));
    else {
        p = 1/(1+MAX(mui, 0)*a);	
	l = gsl_sf_lngamma(MAX(yi+k, 1))-gsl_sf_lngamma(k)-gsl_sf_lngamma(MAX(yi+1, 1));
	l = l + log(p)*k + yi*log( ((1-p)==0) ? 1:(1-p) );
    }
*/    
    return 2*l;
}

int NBinGlm::getfAfAdash(double a, size_t id, double *fAPtr, double *fAdashPtr )
{
    size_t i;
    double yij, mij, dl, ddl, k=1/a;
    *fAPtr = 0;
    *fAdashPtr = 0;
    for ( i=0; i<nRows; i++ ) {
        yij = MAX(gsl_matrix_get(Yref, i, id), 0);
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

int NBinGlm::nbinfit(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O)
{   
    initialGlm(Y, X, O);

    size_t i, j, isConv;
    double yij, mij, vij;
    double disper, a, tol, fA, fAdash, crtold;
    gsl_vector_view b0j, m0j, e0j;

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
//        displayvector(&m0j.vector, "m0j");

        // Get initial phi estimates 
        phi[j]=0;
	iterconv[j]=0;
        if (mmRef->estiMethod==CHI2) {
            disper = fit0.getDisper(j); 
            phi[j] = (disper<1)? 0:1/disper;
	    isConv = (disper<1)? TRUE:FALSE;
	    while (isConv != TRUE) {
                crtold = phi[j];
//                crtold = disper;
	        betaEst(j, 1, &tol, phi[j]); // tol returns |beta-beta.old|
		disper = getDisper(j);
		phi[j] = phi[j]*disper;
                tol = ABS(phi[j]-crtold) + tol;
//                tol = ABS(disper - crtold);                 
	        if ((tol<mintol) | (iterconv[j]==maxiter) | (phi[j]<0)) 
	           break;
                iterconv[j]++;	    
	    }
        }
        else if (mmRef->estiMethod==NEWTON) {
            getfAfAdash(mintol, j, &fA, &fAdash);
            a = (fA>0)?-fA/fAdash:0;
            phi[j] = (a<0)? mintol:a; 
	    isConv = (fA>0)? FALSE:TRUE;
//            printf("fA=%.2f, fAdash=%.2f\n", fA, fAdash);
/*
	    dev[j]=0;
            ll[j]=0;	    
	    for (i=0; i<nRows; i++) {
                yij = gsl_matrix_get(Y, i, j);
                mij = gsl_matrix_get(Mu, i, j);
                dev[j] = dev[j] + devfunc( yij, mij, phi[j] );
//                ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
            }	
*/	    
	    while ( isConv != TRUE ) {
		crtold = phi[j];
//		crtold = dev[j];
//		crtold = ll[j];
	        betaEst(j, 1, &tol, phi[j]);
		gsl_vector_view bj=gsl_matrix_column(Beta, j);
                getfAfAdash(MAX(phi[j], mintol), j, &fA, &fAdash);
                phi[j] = phi[j]-fA/fAdash;
                tol = tol + ABS(phi[j]-crtold);
/*                dev[j] = 0;
//                ll[j]=0;                
	        for (i=0; i<nRows; i++) {
                    yij = gsl_matrix_get(Y, i, j);
                    mij = gsl_matrix_get(Mu, i, j);
                    dev[j] = dev[j] + devfunc( yij, mij, phi[j] );
//                    ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
                }		
		tol = ABS(dev[j]-crtold);
//		tol = ABS(ll[j]-crtold);
*/		                
	        if ((tol<mintol) | (iterconv[j]==maxiter) | (phi[j] < 0)) break;
                iterconv[j]++;	    
	    }
//            printf("phi[%d]=%.2f, ", j, phi[j]);
       }

       if ( phi[j] < 0) { // restore poisson
            phi[j]=0;
            gsl_matrix_set_col (Beta, j, &b0j.vector);
            gsl_matrix_set_col (Mu, j, &m0j.vector);
            gsl_matrix_set_col (Eta, j, &e0j.vector);
/*            dev[j]=0;
            for (i=0; i<nRows; i++) {
               yij = gsl_matrix_get(Y, i, j);
               mij = gsl_matrix_get(Mu, i, j);
               dev[j] = dev[j] + devfunc( yij, mij, phi[j] );
       //        ll[j] = ll[j] + llfunc( yij, mij, phi[j] );
           }			   
*/
       }
       // other properties based on mu and phi
       ll[j] = 0; 
       dev[j] = 0;
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
           dev[j] = dev[j] + devfunc( yij, mij, phi[j] );
       }
       aic[j]=-ll[j]+2*(nParams+1);

   } // end nVar for j loop
  

   // standardize pearson residual rp = rp/sqrt(1-hii)
   // hii is diagonal element in Hat = X*(X'WX)^-1*X'W
   getHat(Xref, wHalf, sqrt1_Hii);
   gsl_matrix_div_elements (Res, sqrt1_Hii);
   subtractMean(Res);
//   displaymatrix(Res, "standardized Pearson resdiduals");

   fit0.releaseGlm();  

   return SUCCESS;    
}

