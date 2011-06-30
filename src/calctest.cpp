// utility functions
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// // 16-June-2011

#include "resampTest.h"

int testStatCalc(mv_mat *H0, mv_mat *H1, mv_Method *mmRef, const int ifcalcH1det, double *stat, gsl_vector *statj)
{
	size_t i, j;
        int s;
        size_t nVars=H0->SS->size1;
        size_t nRows=H0->mat->size1;
	double ss0j, ss1j, sum;
        double logDetSS0, logDetSS1;

	// calc test stats
	if ( mmRef->test == LOGWILK ) {
	   // statj = nRows*(log(diag(ss0)) - log(diag(ss1)));
	    for ( j=0; j<statj->size; j++ ) {
                 ss0j = gsl_matrix_get(H0->SS, j, j);
                 ss1j = gsl_matrix_get(H1->SS, j, j);
	         gsl_vector_set(statj, j, nRows*(log(ss0j)-log(ss1j)));
	    } 
            if ( mmRef->corr == IDENTITY ){ 
            	 sum = 0.0;
                 for ( j=0; j<nVars; j++ )
	             sum += gsl_vector_get(statj, j);
            }
	    else {
	      	 logDetSS0 = log( calcDet(H0->SS) );		 
                 //printf("%4.2f ", *logDetSS0);
		 // used in summary to speed up 
	         if (ifcalcH1det == TRUE){
	            logDetSS1 = log( calcDet(H1->SS) );
		    H1->teststat = logDetSS1;
		 }
		 else
		    logDetSS1 = H1->teststat;
                 sum = nRows * (logDetSS0 - logDetSS1);
            }
	}
	else if ( mmRef->test == HOTELING ) { // Hotelling 
           // statj = diag(ss0)./diag(ss1)-1;  
	   for ( j=0; j<nVars; j++ ) {
               ss0j = gsl_matrix_get(H0->SS, j, j);
	       ss1j = gsl_matrix_get(H1->SS, j, j);
	       gsl_vector_set(statj, j, (double)ss0j/ss1j-1);
	   }
	   if ( mmRef->corr == IDENTITY ){ 
	       sum = 0.0;
               for ( j=0; j<nVars; j++ )
	           sum += gsl_vector_get(statj, j);
           }
	   else {
 	       gsl_matrix *Sa = gsl_matrix_alloc(nVars, nVars);
               gsl_matrix *LU=gsl_matrix_alloc(nVars, nVars);
               gsl_matrix_memcpy(LU, H1->SS); 
               gsl_matrix_memcpy(Sa, H0->SS);
               for (i=0;i<nVars; i++)
               for (j=i+1; j<nVars; j++){ 
                   gsl_matrix_set(LU, i, j, gsl_matrix_get(H1->SS, j, i));
                   gsl_matrix_set(Sa, i, j, gsl_matrix_get(H0->SS, j, i)); 
               }
               gsl_permutation *p = gsl_permutation_alloc(nVars);
               gsl_linalg_LU_decomp (LU, p, &s);
               for (j=0; j<nVars; j++) { 
                  gsl_vector_view xj=gsl_matrix_column (Sa, j);
                  gsl_linalg_LU_svx (LU, p, &xj.vector);
               } 
	       sum = 0.0;
	       for ( j=0; j<nVars; j++ ) 
	           sum += gsl_matrix_get(Sa, j, j)-1;
               //printf("multstat=%.4f\n", sum);

	       gsl_matrix_free(Sa);
               gsl_permutation_free(p);
               gsl_matrix_free(LU);
           }
	}
	else {
	   GSL_ERROR("Invalid test type", GSL_EINVAL);
	   return GSL_EINVAL;
	}

        // copy results
	*stat = sum;

	return 0;
}

int calcSS(gsl_matrix *Y, mv_mat *Hat, mv_Method *mmRef, const int ifcalcHat, const int ifcalcCoef, const int ifcalcSS)
{
    size_t j; 
    size_t nP=Hat->X->size2;
    size_t nRows=Hat->mat->size1;
    size_t nVars=Hat->SS->size1;

     // It is possible later to feed data with more varialbes (columns) 
     // than observations (rows). So better use SVD than QR
     if (ifcalcHat == TRUE) {
	gsl_matrix *U=gsl_matrix_alloc(nRows, nP);
        gsl_vector *t=gsl_vector_alloc(MIN(nRows, nP));

	gsl_vector *trj=gsl_vector_alloc(nRows);
	gsl_vector *tcj=gsl_vector_alloc(nP);

	gsl_matrix *Q=gsl_matrix_alloc(nRows, nRows);
        gsl_matrix *R=gsl_matrix_alloc(nRows, nP);

	gsl_matrix_memcpy(U, Hat->X);
	gsl_linalg_QR_decomp (U, t);

	gsl_vector_view yj, cj, rj;
        for (j=0; j<nVars; j++){
	    yj = gsl_matrix_column(Y, j);
	    if (ifcalcCoef == TRUE)
	       cj = gsl_matrix_column(Hat->Coef, j);
	    else
	       cj = gsl_vector_subvector (tcj, 0, nP);
	    rj = gsl_matrix_column(Hat->Res, j);
	    gsl_linalg_QR_lssolve (U, t, &yj.vector, &cj.vector, &rj.vector); 
        }

	gsl_linalg_QR_unpack (U, t, Q, R);	
	// Q'\R ie., solve Rx=Q'
        gsl_linalg_QR_decomp (R, t);
	gsl_vector_view mj, qj;
	gsl_matrix_view Xsub;
	for (j=0; j<nRows; j++){
	    qj=gsl_matrix_row(Q, j); // Q'[j, :]
	    gsl_linalg_QR_lssolve (R, t, &qj.vector, tcj, trj); 
	    // X*(Q'\R)
	    mj=gsl_matrix_subcolumn(Hat->mat, j, j, nRows-j);	    
	    Xsub=gsl_matrix_submatrix(Hat->X, j, 0, nRows-j, nP);
	    gsl_blas_dgemv (CblasNoTrans,1.0,&Xsub.matrix,tcj,0.0,&mj.vector);
	}
//        displaymatrix (Hat->Res, "Hat.Res from calcHat");
	gsl_matrix_free(U);
	gsl_matrix_free(Q);
	gsl_matrix_free(R);
	gsl_vector_free(t);
	gsl_vector_free(trj);
	gsl_vector_free(tcj);
    }

    // Compute SS
    gsl_matrix *I=NULL;
    if ( ifcalcSS == TRUE) {
       // generate the same residual as above
       if (ifcalcHat == FALSE) {
          I = gsl_matrix_alloc(nRows, nRows);
          gsl_matrix_set_identity(I);
          gsl_matrix_sub(I, Hat->mat);
          gsl_blas_dsymm (CblasLeft, CblasLower, 1.0, I, Y, 0.0, Hat->Res);
 //         displaymatrix (Hat->Res, "Hat.Res from (I-Hat)*Y");
       }
       // rcalc:  rcalc(Hat->Res, shrink_param, corr, Hat->SS);
        gsl_vector *e = gsl_vector_alloc(nRows);
        gsl_vector_set_all (e, 1.0); 

        gsl_matrix_set_zero(Hat->SS);
        double mean, dd;	
        // subtract mean = Rs^T e(12) from res
	for ( j=0; j<nVars; j++ ) {
            gsl_vector_view rj=gsl_matrix_column(Hat->Res, j);
            gsl_blas_ddot (&rj.vector, e, &mean); 
            gsl_vector_add_constant (&rj.vector, -mean/nRows); 
            if ( mmRef->corr == IDENTITY ) {
               gsl_blas_ddot (&rj.vector, &rj.vector, &dd); 
               if ( dd < TOL ) dd = 0.5*M_1_PI;
               gsl_matrix_set(Hat->SS, j, j, dd);
            }
	} 

        // compute the covariance matrix SS= (res'*res)/Rows-1
        if ( mmRef->corr != IDENTITY ){
           gsl_blas_dsyrk(CblasLower, CblasTrans,1.0, Hat->Res, 0.0, Hat->SS);
           gsl_matrix_scale ( Hat->SS, 1.0/(double)(nRows-1) );
           if ( mmRef->corr == SHRINK ) {
              gsl_vector_view dj=gsl_matrix_diagonal (Hat->SS);
              gsl_vector_scale (&dj.vector, 1.0/(mmRef->shrink_param));
           }
        } 
 
       gsl_vector_free(e);
       //displaymatrix(Hat->SS, "Hat->SS");
   }
   if ( I!=NULL )
      gsl_matrix_free(I);

   return 0;

}

double calcDet(gsl_matrix *SS)
{
     // SS is a nVars x nVars real symmetric matric
     int nVars = SS->size1;
     double result;
/*     // det(A) = prod(eig(A))
     gsl_eigen_symm_workspace *ws = gsl_eigen_symm_alloc(nVars);
     gsl_vector *eval=gsl_vector_alloc(nVars);    
     int j;
     double logDetSS=0.0;
     gsl_eigen_symm(SS, eval, ws);
//     displayvector(eval, "SS eigen values");
     for ( j=0; j<nVars; j++ )
         logDetSS += log(gsl_vector_get(eval, j));

     gsl_eigen_symm_free(ws);
     gsl_vector_free(eval);
     return logDetSS;
  */

   // OR: det(A) = diag(LU)
     // fill SS
     gsl_matrix *LU=gsl_matrix_alloc(nVars, nVars);
     gsl_matrix_memcpy(LU, SS); 
     int i, j, s;
     for (i=0;i<nVars; i++)
     for (j=i+1; j<nVars; j++)
         gsl_matrix_set(LU, i, j, gsl_matrix_get(LU, j, i));

     gsl_permutation *p = gsl_permutation_alloc(nVars);
     gsl_linalg_LU_decomp (LU, p, &s);
     result = gsl_linalg_LU_det (LU, s);

     gsl_permutation_free(p);
     gsl_matrix_free(LU);

     return result;

}

int is_sym_matrix(const gsl_matrix *mat)
{
     size_t i, j;
     if ( mat->size1 == mat->size2 ) {
        // test the upper triangle
        for ( i=0; i<mat->size1; i++)
            for ( j=i+1; j<mat->size2; j++ )
                if ( gsl_matrix_get(mat, i, j) != 0 )
	           return 0;
        return 1;
     }
     else
        return 0;
}

int subX(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    size_t j, k=0;
    size_t nP=ref->size;

    for (j=0; j<nP; j++) {
        if ( (size_t)gsl_vector_get(ref, j) > 0 ){
           gsl_vector_view col = gsl_matrix_column (X, j);
           gsl_matrix_set_col(Xi, k, &col.vector);
            k++;
	}    
     }
     return 0;
}

int subX1(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    size_t j, k=0;
    size_t nP=ref->size;

    for (j=0; j<nP; j++) {
        if ( (size_t)gsl_vector_get(ref, j) == 0 ){
            gsl_vector_view col = gsl_matrix_column (X, j);
            gsl_matrix_set_col(Xi, k, &col.vector);
            k++;
        }
    }

    return 0;
}


int subX2(gsl_matrix *X, int id, gsl_matrix *Xi)
{   
    size_t j, k=0;
    size_t nParam=X->size2;

    for (j=0; j<nParam; j++) {
        if ( j!=id ){
           gsl_vector_view col = gsl_matrix_column (X, j);
           gsl_matrix_set_col(Xi, k, &col.vector);
           k++;
	}   
    }
    return SUCCESS;
}

int subXrow(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    size_t j, k=0;
    size_t n=X->size1;

    for (j=0; j<n; j++) {       
        if ( (size_t)gsl_vector_get(ref, j)==0 ){
           gsl_vector_view row = gsl_matrix_row (X, j);
           gsl_matrix_set_row(Xi, k, &row.vector);
           k++;
        }
     }
     return 0;
}     

int subXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi)
{   
    size_t j, k=0;
    size_t n=X->size1;

    for (j=0; j<n; j++) {       
        if ( (size_t)gsl_vector_get(ref, j)>0 ){
           gsl_vector_view row = gsl_matrix_row (X, j);
           gsl_matrix_set_row(Xi, k, &row.vector);
           k++;
        }
     }
     return 0;
}     

int subXrow1(gsl_matrix *X, gsl_vector *ref0, gsl_vector *ref1, gsl_matrix *Xi)
{ 
    size_t j, k=0, id0, id1;
    size_t n=X->size1;

    for (j=0; j<n; j++) {       
        id0 = (size_t) gsl_vector_get(ref0, j);
        id1 = (size_t) gsl_vector_get(ref1, j);
        if ( (id0!=id1) & (id0==0) ){
           gsl_vector_view row = gsl_matrix_row (X, j);
           gsl_matrix_set_row(Xi, k, &row.vector);
           k++;
        }
     }
     return 0;
}     

int calcAdjustP(const int punit, const int nVars, double *bj, double *sj, double *pj, gsl_permutation *sortid)
{
    int k;
    
 //   printf("multiple test procedure [%d]\n", punit);
    if (punit == UNADJUST){
       for (k=0; k<nVars; k++)
       if (*(bj+k) >= *(sj+k))
          *(pj+k)=*(pj+k)+1;
     }
    else if (punit == SINGLESTEP) {
       gsl_vector_view bStatj=gsl_vector_view_array (bj, nVars);
       double maxstat=gsl_vector_max(&bStatj.vector);
       for (k=0; k<nVars; k++)
       if ( maxstat >= *(sj+k) )
           *(pj+k)=*(pj+k)+1;
     }
    else if (punit == FREESTEP){
       size_t sid, sid0;
       // successive maxima
       //gsl_permutation_fprintf(stdout, sortid, " %u");
 //      printf("\n");
       for (k=0; k<nVars; k++){
           sid = gsl_permutation_get(sortid, nVars-1-k);
//	   printf("%d ", (size_t)sid);
           if ( k>0 ) {
//	      printf("%d ", (size_t)sid0);
	      *(bj+sid)=MAX(*(bj+sid), *(bj+sid0));
           }   
//	   printf("%d ", (size_t)sid);
	   if ( *(bj+sid) >= *(sj+sid) ) 
	      *(pj+sid)=*(pj+sid)+1;
//	   if ( k> 0)    
//              printf("%.2f ", *(bj+sid0));   
//	   printf("%.2f ", *(bj+sid));   
//           printf("%.2f ", *(sj+sid));	  
//           printf("%d\n", (size_t)*(pj+sid));
	   sid0 = sid;   
    }  }
/*    else if (punit == STEPUP) {
       // successive minima
       for (k=0; k<nVars; k++){
           sid = gsl_permutation_get(sortid, k);
	   if ( k>0 ) {
	      *(bj+sid)=MIN(*(bj+sid), *(bj+sid0));
	   }
//           printf("%.3f ", (double)*(bj+sid));
	   //printf("%d ", (size_t)sid);
	   if (*(bj+sid) >= *(sj+sid))
	      *(pj+sid)=*(pj+sid)+1;
	   sid0 = sid;   
    }  }
    else if (punit == NOMONO) {
    }
//    printf("\n");
*/
    return 0;

}


int reinforceP(double *puj, size_t nVars, gsl_permutation *sortid)
{
// re-enforce monotonicity for univaraite tests
    size_t j, sid, sid0;
    for (j=1; j<nVars; j++) {
        sid=gsl_permutation_get(sortid, j);
	sid0=gsl_permutation_get(sortid, j-1);
	*(puj+sid) = MAX (*(puj+sid), *(puj+sid0));
    }	
    return 0;
}


int rcalc( gsl_matrix *Res, double shrink_param, int corr, gsl_matrix *SS)
{
    size_t j;
    size_t nRows = Res->size1;
    size_t nVars = Res->size2;
    
    gsl_vector *e = gsl_vector_alloc(nRows);
    gsl_vector_set_all (e, 1.0); 

    gsl_matrix_set_zero(SS);
    double mean, dd;	
    // subtract mean = Rs^T e(12) from res
    for ( j=0; j<nVars; j++ ) {
        gsl_vector_view rj=gsl_matrix_column(Res, j);
        gsl_blas_ddot (&rj.vector, e, &mean); 
        gsl_vector_add_constant (&rj.vector, -mean/nRows); 
        if ( corr == IDENTITY ) {
            gsl_blas_ddot (&rj.vector, &rj.vector, &dd); 
            if ( dd < 1e-10 ) dd = 0.5*M_1_PI;
            gsl_matrix_set(SS, j, j, dd);
        }
    } 
   // compute the covariance matrix SS= (res'*res)/Rows-1
    if ( corr != IDENTITY ){
        gsl_blas_dsyrk(CblasLower, CblasTrans,1.0, Res, 0.0, SS);
        gsl_matrix_scale ( SS, 1.0/(double)(nRows-1) );
        if ( corr == SHRINK ) {
           gsl_vector_view dj=gsl_matrix_diagonal (SS);
	   for (j=0; j<nVars; j++) {
	       dd = gsl_vector_get (&dj.vector, j);	       
	       if ( dd < 1e-10 ) // account for zero variances
	          gsl_vector_set (&dj.vector, j, 1.0/shrink_param); 
	       else  
	          gsl_vector_set (&dj.vector, j, dd/shrink_param); 
	   }    
        }
    } 
 
    gsl_vector_free(e);
    return 0;
}

int getHat(gsl_matrix *X, gsl_matrix *W, gsl_matrix *Hat) 
{
     double hii;
 
     size_t i, j;
     size_t nRows = X->size1;
     size_t nParam = X->size2;
     size_t nVars = Hat->size2;

     gsl_matrix *Xw = gsl_matrix_alloc(nRows, nParam);
     gsl_matrix *dW2 = gsl_matrix_alloc(nRows, nRows);
     gsl_matrix_set_zero(dW2);
     gsl_vector *t = gsl_vector_alloc(MIN(nRows, nParam));
     gsl_vector *twj = gsl_vector_alloc(nParam);
     gsl_vector *rj = gsl_vector_alloc(nRows);
     gsl_vector_view diagW, wj, qwj, xi;

     for ( j=0; j<nVars; j++) {
          gsl_matrix_memcpy (Xw, X);
          diagW = gsl_matrix_diagonal(dW2);
          wj = gsl_matrix_column(W, j);
          gsl_vector_memcpy(&diagW.vector, &wj.vector);
	  gsl_blas_dtrmm (CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, dW2, Xw); 
          gsl_linalg_QR_decomp (Xw, t); 
          // twj=(X'WX)^-1*XW = R\Q'wHalf where Xw = Q*R
          for (i=0; i<nRows; i++){
              qwj=gsl_matrix_column(dW2, i);
              gsl_linalg_QR_lssolve (Xw, t, &qwj.vector, twj, rj);

              xi=gsl_matrix_row(X, i);
	      gsl_blas_ddot(&xi.vector, twj, &hii);  //hii = xi*twj
              gsl_matrix_set(Hat, i, j, sqrt(1-hii));
          }         
     }
    
     gsl_matrix_free(Xw);
     gsl_matrix_free(dW2);
     gsl_vector_free(t);
     gsl_vector_free(twj);
     gsl_vector_free(rj);

     return 0;
}

int invLSQ(gsl_matrix *A, gsl_vector *b, gsl_vector *x)
{ 
   size_t dim1 = A->size1;
   size_t dim2 = A->size2;
   gsl_vector *t = gsl_vector_alloc(MIN(dim1, dim2));
   gsl_vector *r = gsl_vector_alloc(dim1);

   gsl_linalg_QR_decomp (A, t);
   gsl_linalg_QR_lssolve (A, t, b, x, r);

   gsl_vector_free(t);
   gsl_vector_free(r);
   
   return 0;
}

// the manyglm version of rcalc
int GetR(gsl_matrix *Res, size_t corr, gsl_vector *glmshrink, size_t k, gsl_matrix *R)
{
    if (corr == IDENTITY) {
       gsl_matrix_set_identity (R);
    }
    else {
       size_t j;
       size_t nRows = Res->size1;
       size_t nVars = Res->size2;
       double std, lambda;
       gsl_matrix *Sd =gsl_matrix_alloc(nVars, nVars);
       gsl_vector *s = gsl_vector_alloc(nVars);

       // get covariance matrix SS
       // note that residuals have already had means subtracted
       gsl_matrix_set_zero (R);
       gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0, Res, Res,0.0, R);
       gsl_matrix_scale (R, 1.0/(double)(nRows-1)); // SS

       gsl_matrix_set_all (Sd, 1.0);
       for ( j=0; j<nVars; j++ ) {               
           std = sqrt(gsl_matrix_get(R, j, j));          
           // account for zero variance
           gsl_vector_set(s, j, (std>TOL)?std:1/(2*M_1_PI)); 
       }
              
       gsl_blas_dsyr (CblasLower, 1.0, s, Sd);       
       gsl_matrix_set_zero(Sd);       
       gsl_blas_dger (1.0, s, s, Sd);
                     
       gsl_matrix_div_elements (R, Sd); // RR = SS./Sd              
       gsl_vector_view d = gsl_matrix_diagonal (R);       
       gsl_vector_set_all (&d.vector, 1.0); // ensure diag(R) = 1       
       //  displaymatrix(R, "corr-coefficents");
                                                        
       if ( corr == SHRINK ) { 
          lambda = gsl_vector_get(glmshrink, k);
          gsl_matrix_scale (R, lambda);
          gsl_vector_add_constant(&d.vector,(double)(1.0-lambda));          
       }
                                                        
       gsl_matrix_free(Sd);       
       gsl_vector_free(s);       
       
    }
    
    return SUCCESS;
    
}

int subtractMean(gsl_matrix *dat)
{

    size_t nRows = dat->size1;
    size_t nVars = dat->size2;

//    gsl_matrix *t1, *Mean;
//    t1 = gsl_matrix_alloc(nRows, 1);
//    gsl_matrix_set_all (t1, 1.0);
//    Mean = gsl_matrix_alloc(nRows, nVars);
//    GetMean(t1, dat, Mean);
//    gsl_matrix_sub(dat, Mean);

    gsl_matrix *tmp = gsl_matrix_alloc(nRows, nVars);
    gsl_vector *t1 = gsl_vector_alloc(nRows);
    gsl_vector_set_all (t1, 1.0);
    
    gsl_vector *mean = gsl_vector_alloc(nVars);    
    gsl_blas_dgemv (CblasTrans, 1.0, dat, t1, 0.0, mean);
    gsl_vector_scale ( mean, (double)(1/nRows) );

    gsl_matrix_set_zero (tmp);
    gsl_blas_dger (1.0, t1, mean, tmp); 

    // subtract mean from dat
    gsl_matrix_sub(dat, tmp);

    gsl_vector_free(mean);
    gsl_vector_free(t1); 
    gsl_matrix_free(tmp);

    return SUCCESS;
}

