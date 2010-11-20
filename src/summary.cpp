// calculates a summary of a multivariate linear model
// Author: Yi Wang (yi dot wang at computer dot org)
// Last modified on 20-April-2010


#include "resampTest.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_statistics.h> // only gsl_stats_vairaince used in rcalc 
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_double.h>

// Use gsl functions as much as poosible to increase speed and stabilitiy

Summary::Summary(mv_Method *mm, gsl_matrix *Y, gsl_matrix *X):mmRef(mm), Yref(Y), Xref(X)
{
    size_t i, j;
    nRows=Yref->size1, nVars=Yref->size2, nParam=Xref->size2;
 
    // initialize public variables: stats    
    multstat=(double *)malloc((nParam+1)*sizeof(double));
    Pmultstat=(double *)malloc((nParam+1)*sizeof(double));
    unitstat=gsl_matrix_alloc(nParam+1, nVars);
    Punitstat=gsl_matrix_alloc(nParam+1, nVars);

    for (j=0; j<nParam+1; j++)
         *(Pmultstat+j)=0.0; 	
    gsl_matrix_set_zero (Punitstat);
    
    // initialize private variables 
    bMultStat = (double *)malloc((nParam+1)*sizeof(double));
    bUnitStat = gsl_matrix_alloc(nParam+1, nVars);

    Hats = (mv_mat *)malloc((nParam+2)*sizeof(mv_mat));
    sortid=(gsl_permutation **)malloc((nParam+1)*sizeof(gsl_permutation *));
    gsl_vector *ref=gsl_vector_alloc(nParam);
    gsl_vector_set_all (ref, 1.0);

    for (i=0; i<nParam+2; i++ ) {
        Hats[i].mat=gsl_matrix_alloc(nRows, nRows);
        Hats[i].SS=gsl_matrix_alloc(nVars, nVars);
        Hats[i].Res=gsl_matrix_alloc(nRows, nVars);
        Hats[i].Y = gsl_matrix_alloc(nRows, nVars);

        // calc test statistics
	if ( i == 0 ) { // the complete model
           Hats[i].X = gsl_matrix_alloc(nRows, nParam);
           gsl_matrix_memcpy (Hats[i].X, Xref);
	   Hats[i].Coef = gsl_matrix_alloc(nParam, nVars);
           calcSS(Yref, &(Hats[i]),mmRef,TRUE,TRUE,TRUE);
        }
	else {
	   if ( i == 1 ) { // response mean
              Hats[i].X = gsl_matrix_alloc(nRows, 1);
              gsl_matrix_set_all (Hats[i].X, 1.0);
	      Hats[i].Coef = gsl_matrix_alloc(1, nVars);
           }
	   else { // subtract a term to test significance
	      if ( nParam == 1 ) // intercept model
	      {
                  Hats[i].X = gsl_matrix_alloc(nRows, 1);
                  Hats[i].Coef=gsl_matrix_alloc(1, nVars);
		  gsl_matrix_memcpy(Hats[i].X, Xref);
	      }
	      else {
                  Hats[i].X = gsl_matrix_alloc(nRows, nParam-1);
                  Hats[i].Coef=gsl_matrix_alloc(nParam-1, nVars);
                  gsl_vector_set (ref, i-2, 0);
                  subX(Xref, ref, Hats[i].X);
                  gsl_vector_set (ref, i-2, 1);
	      }
	   }	 
          calcSS(Yref, &(Hats[i]),mmRef,TRUE,TRUE,TRUE);
          gsl_vector_view sj=gsl_matrix_row(unitstat, i-1);
          testStatCalc(&(Hats[i]),&(Hats[0]),mmRef,TRUE,(multstat+i-1),&sj.vector);
          // sortid
          sortid[i-1] = gsl_permutation_alloc(nVars);
          gsl_sort_vector_index (sortid[i-1], &sj.vector); 
          // rearrange sortid in descending order
          gsl_permutation_reverse (sortid[i-1]);
//	     gsl_permutation_fprintf(stdout, sortid[i-1], " %u");
//	     printf("\n");
       }  
   }
   
    // fit = Y- resi 
    for (i=0; i<nParam+2; i++) {
        gsl_matrix_memcpy (Hats[i].Y, Yref);
        gsl_matrix_sub (Hats[i].Y, Hats[i].Res);	
    }   

    calcR2();

    // initialize resampling indices 
    bootID = NULL;

    gsl_vector_free(ref);
//    printf("Summary test initialized.\n");
}

Summary::~Summary(){
}

void Summary::releaseSummary()
{
    free(multstat);
    free(Pmultstat);
    gsl_matrix_free(unitstat);
    gsl_matrix_free(Punitstat);
// The above commented out s.t. the results are passed out to Rcpp wrapper

    size_t i;
    for ( i=0; i<nParam+2; i++ ){
        gsl_matrix_free(Hats[i].mat);
        gsl_matrix_free(Hats[i].SS);
	gsl_matrix_free(Hats[i].Res);
        gsl_matrix_free(Hats[i].Coef);
	gsl_matrix_free(Hats[i].X);
	gsl_matrix_free(Hats[i].Y);
    }
    free(bMultStat);
    gsl_matrix_free(bUnitStat);

    if ( bootID != NULL )
       gsl_matrix_free(bootID);

    for (i=0; i<nParam+1; i++)
        gsl_permutation_free(sortid[i]);
    free(sortid);

//    printf("Summary test released.\n");

}

void Summary::display(void)
{
    size_t i, j, k, nk, lk;
    printf("Summary of fitting (resampling under H1):\n");
    // significance test
    printf("Significance Test:\n");
    if (mmRef->test == LOGWILK)
        printf("Explan\t LR value\tPr(>LR)\n");
    else
        printf("Explan\t F value\tPr(>F)\n");

    for ( i=0; i<nParam; i++ )
       printf("  %d\t %.3f\t\t %.3f\n", i, multstat[i+1], Pmultstat[i+1]);

    if ( mmRef->punit!=NONE) {
       // Significance univariate tests
       printf("Univariate Tests for significance:\n");
       nk = (int)floor((double)nParam/WRAP);    
       for (k=0; k<nk; k++) {
           printf("\t\t");    
           for (j=k*WRAP; j<(k+1)*WRAP; j++) {
               printf("[Explan %d]\t", j);
           }	
           printf("\n");
           for ( i=0; i<nVars; i++ ){ 
               printf("[Respons %d]\t", i);
               for (j=k*WRAP; j<(k+1)*WRAP; j++) {
                   printf("%.3f(%.3f)\t", gsl_matrix_get(unitstat, j+1, i), gsl_matrix_get(Punitstat, j+1, i));
               }
               printf("\n");
            }
            printf("\n");
        }
        lk = nParam-nk*WRAP;
        if (lk>0) {
           printf("\t\t");    
           for (j=0; j<lk; j++) {
               printf("[Explan %d]\t", nk*WRAP+j); 
           }
           printf("\n");
           for ( i=0; i<nVars; i++ ){ 
               printf("[Respons %d]\t", i);
               for (j=0; j<lk; j++) {
                   printf("%.3f(%.3f)\t", gsl_matrix_get(unitstat, nk*WRAP+j+1, i), gsl_matrix_get(Punitstat, nk*WRAP+j+1, i));
               }
               printf("\n");
           }
        }
    }
    // Overall statistics
    printf("\n###################\n");
    printf("\nOverall statistics:\n");
    if ( mmRef->rsquare == HOOPER )
       printf("Hooper's R-squared: ");
    else
       printf("Vector's R-squared: ");
    printf("%.3f\n", R2);
    if (mmRef->test == LOGWILK)
       printf("Likelihood Ratio ");
    else
       printf("Lawley-Hotelling Trace ");
    printf("(p-value): %.3f (%.3f)\n", multstat[0], Pmultstat[0]);    
    if (mmRef->punit!=NONE){
        printf("Univariate Tests for overall:\n");
        // note the change of the display direction
        nk = (int)floor((double)nVars/WRAP);    
        for (k=0; k<nk; k++) {
            printf("\t\t");    
            for (j=k*WRAP; j<(k+1)*WRAP; j++) {
                printf("[Explain %d]\t", j);   
            } 	
            printf("\n");
            if (mmRef->test==LOGWILK)
	       printf("LR (Pr):\t");
            else
	       printf("F(Pr):\t");
            for (j=k*WRAP; j<(k+1)*WRAP; j++) {
                printf("%.3f(%.3f)\t", gsl_matrix_get(unitstat, 0, j), gsl_matrix_get(Punitstat, 0, j));
            }
            printf("\n");
         }
         lk = nVars-nk*WRAP;
         if (lk>0) {
            printf("\t\t");    
            for (j=0; j<lk; j++)
                printf("[Explain %d]\t", nk*WRAP+j); 
	    printf("\n");
            if (mmRef->test==LOGWILK)
	        printf("LR (Pr):\t");
	    else
	        printf("F(Pr):\t");
            for (j=0; j<lk; j++)
                printf("%.3f(%.3f)\t", gsl_matrix_get(unitstat,0,nk*WRAP+j), gsl_matrix_get(Punitstat,0,nk*WRAP+j));
            printf("\n");
         }    
     }	 
}

int Summary::resampTest(void)
{
  //  printf("Start resampling test ...\n");
    size_t i, j, id;
    size_t maxiter= mmRef->nboot; 
//    printf("maxiter=%d\n", maxiter);
    double hii, score;

    gsl_matrix *bX, *bY;
    bY = gsl_matrix_alloc(nRows, nVars);
    bX = gsl_matrix_alloc(nRows, nParam);

    gsl_rng *rnd;
    const gsl_rng_type *T; 
    gsl_rng_env_setup();
    T=gsl_rng_default;
    rnd=gsl_rng_alloc(T);

    // initialize permid
    size_t *permid=NULL;
    if ( bootID == NULL ) {
       if ( mmRef->resamp == PERMUTE ) {
           permid = (size_t *)malloc(nRows*sizeof(size_t));
	   for (i=0; i<nRows; i++)
	       permid[i]=i;
       }
    }
    nSamp = 0;
    // resampling options 
    if (mmRef->resamp == CASEBOOT) {
       for (i=0; i<maxiter; i++) {
           for ( j=0; j<nRows; j++ ){
	       // resampling index
 	       if (bootID == NULL) 
	          id = gsl_rng_uniform_int(rnd, nRows);
               else 
	          id = (size_t) gsl_matrix_get(bootID, i, j);
               // resample Y and X
               gsl_vector_view Yj=gsl_matrix_row(Yref, id);
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view Xj=gsl_matrix_row(Xref, id);
               gsl_matrix_set_row (bX, j, &Xj.vector); 
	    }
           smrycase(bY, bX);
        } 
     } 
    else if (mmRef->resamp == RESIBOOT) {
        for (i=0; i<maxiter; i++){ 
           for (j=0; j<nRows; j++){
           // resampling index
 	       if (bootID == NULL) 
	          id = gsl_rng_uniform_int(rnd, nRows);
               else 
	          id = (size_t) gsl_matrix_get(bootID, i, j);
               // bootr by resampling resi=(Y-fit)
               gsl_vector_view Yj=gsl_matrix_row(Yref, id);
               gsl_vector_view Fj=gsl_matrix_row(Hats[0].Y, id); //fit_alt
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view bootr=gsl_matrix_row(bY, j);
               gsl_vector_sub (&bootr.vector, &Fj.vector);  
               if (mmRef->student==TRUE) {
                  hii = gsl_matrix_get(Hats[0].mat, id, id);
                  gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
               } 
	   }  
	  // bY = bootr for resampling under H1 
//	  displaymatrix(Hats[0].Y, "fit");
          smryresi(bY);
       } 
    }
   else if (mmRef->resamp == SCOREBOOT) {
       for (i=0; i<maxiter; i++) {
          for ( j=0; j<nRows; j++ ) {
             // random score
	     if ( bootID == NULL )
	         score = gsl_ran_ugaussian (rnd); 
	     else
	         score = (double)gsl_matrix_get(bootID, i, j);
             // bootr = (Y - fit)*score 
             gsl_vector_view Yj=gsl_matrix_row(Yref, j);
             gsl_vector_view Fj=gsl_matrix_row(Hats[0].Y, j);
             gsl_matrix_set_row (bY, j, &Yj.vector);
             gsl_vector_view bootr=gsl_matrix_row(bY, j);
             gsl_vector_sub (&bootr.vector, &Fj.vector); 
             if (mmRef->student==TRUE) {
                hii = gsl_matrix_get(Hats[0].mat, j, j);
                gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
             }
             gsl_vector_scale (&bootr.vector, score);
 	   } 
	  // bY = bootr for resampling under H1 
          smryresi(bY);
       } 
    } 
   else if ( mmRef->resamp == PERMUTE ) { 
       gsl_matrix_add_constant (Punitstat, 1.0); 
       for (i=0; i<maxiter-1; i++) { //999
          if (bootID == NULL ) 
             gsl_ran_shuffle(rnd, permid, nRows, sizeof(size_t));
          // get bootr by permuting resi:Y-fit
          for (j=0; j<nRows; j++){
 	      if (bootID == NULL) 
	         id = permid[j];
              else 
	         id = (size_t) gsl_matrix_get(bootID, i, j);
              // bootr by resampling resi=(Y-fit)
              gsl_vector_view Yj=gsl_matrix_row(Yref, id);
              gsl_vector_view Fj=gsl_matrix_row(Hats[0].Y, id);
              gsl_matrix_set_row (bY, j, &Yj.vector);
              gsl_vector_view bootr=gsl_matrix_row(bY, j);
              gsl_vector_sub (&bootr.vector, &Fj.vector); 
              if (mmRef->student==TRUE) {
                 hii = gsl_matrix_get(Hats[0].mat, id, id);
                 gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
              }
          }
	  // bY = bootr for resampling under H1 
          smryresi(bY);
       }  
    }
   else 
       GSL_ERROR("Invalid resampling option", GSL_EINVAL);

   // p-values 
   size_t sid, sid0;
   double *pj;
   for (i=0; i<nParam+1; i++) {
//      printf("Pmultstat[%d]=%.3f\n", i, Pmultstat[i]);
//      printf("nSamp=%d\n", nSamp);
      Pmultstat[i] = (double) Pmultstat[i]/nSamp;
      pj = gsl_matrix_ptr (Punitstat, i, 0);
      if ( mmRef->punit == FREESTEP ) {
         for (j=1; j<nVars; j++) {
	     sid=gsl_permutation_get(sortid[i], j);
	     sid0=gsl_permutation_get(sortid[i], j-1);
	     *(pj+sid) = MAX (*(pj+sid), *(pj+sid0));
	  }
      }
      if ( mmRef->punit == STEPUP ) {
         for (j=2; j<nVars; j++) {
	     sid=gsl_permutation_get(sortid[i], nVars-j);
	     sid0=gsl_permutation_get(sortid[i], nVars-j+1);
	     *(pj+sid) = MIN (*(pj+sid), *(pj+sid0));
	  }
      }
      for (j=0; j<nVars; j++)
         *(pj+j) = (double)*(pj+j)/nSamp;
   }

   // free memory
   gsl_matrix_free(bX);
   gsl_matrix_free(bY);
   gsl_rng_free(rnd);
   if (permid!=NULL)
      free(permid);

   return 0;

}

int Summary::smrycase(gsl_matrix *bY, gsl_matrix *bX)
{
   size_t j;
   // if Y col is all zeros
   for ( j=0; j<nVars; j++ ){
       gsl_vector_view colj = gsl_matrix_column(bY, j);
       if ( gsl_vector_isnull(&colj.vector) == TRUE )
          return GSL_ERANGE;
   }

   size_t i;
   double *sj, *pj;

   //Y = X*coef
   gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1.0,bX,Hats[0].Coef,0.0,Hats[0].Y);
   // bZ = bY - bX*coef;
   gsl_matrix_sub (Hats[0].Y, bY);
   gsl_matrix_scale (Hats[0].Y, -1.0);

   // calc overall stats
   gsl_matrix_memcpy (Hats[0].X, bX);
   calcSS(Hats[0].Y, &(Hats[0]), mmRef, TRUE, FALSE, TRUE);
   calcSS(Hats[0].Y, &(Hats[1]), mmRef, FALSE, FALSE, TRUE);
   gsl_vector_view buj=gsl_matrix_row(bUnitStat, 0);
   testStatCalc(&(Hats[1]),&(Hats[0]), mmRef,TRUE,bMultStat,&buj.vector);

   // calc overal test unistat 
   if ( bMultStat[0] >= multstat[0] )
       Pmultstat[0]++;
   sj = gsl_matrix_ptr (unitstat, 0, 0);
   pj = gsl_matrix_ptr (Punitstat, 0, 0);
   calcAdjustP(mmRef->punit, nVars, &buj.vector, sj, pj, sortid[0]);

   // calc significance stats
   gsl_vector *ref=gsl_vector_alloc(nParam);
   gsl_vector_set_all (ref, 1.0);
   for (i=1; i<nParam+1; i++) {    
      gsl_vector_set (ref, i-1, 0);
      subX(bX, ref, Hats[i+1].X);
      gsl_vector_set (ref, i-1, 1);
     
      calcSS(Hats[0].Y, &(Hats[i+1]), mmRef, TRUE, FALSE, TRUE);
      buj = gsl_matrix_row (bUnitStat, i);
      testStatCalc(&(Hats[i+1]),&(Hats[0]),mmRef,FALSE,(bMultStat+i),&buj.vector);

      // count data related to P-values
      if ( bMultStat[i] >= multstat[i])
         Pmultstat[i]++;
      sj = gsl_matrix_ptr (unitstat, i, 0);
      pj = gsl_matrix_ptr (Punitstat, i, 0);
      calcAdjustP(mmRef->punit, nVars, &buj.vector, sj, pj, sortid[i]);
   }  
   nSamp++;
   gsl_vector_free(ref);
   return 0;

}

int Summary::smryresi(gsl_matrix *bY)
{
    size_t i;
    // count the right-hand tails
    calcSS(bY, &(Hats[0]), mmRef, FALSE, FALSE, TRUE);
    calcSS(bY, &(Hats[1]), mmRef, FALSE, FALSE, TRUE);
    gsl_vector_view buj=gsl_matrix_row(bUnitStat, 0);
    testStatCalc(&(Hats[1]),&(Hats[0]),mmRef,TRUE,(bMultStat+0),&buj.vector);
    // calc overal test unistat 
 //   printf("%.3f %.3f\n", bMultStat[0], multstat[0]);
    if ( bMultStat[0] >= multstat[0] ){
       Pmultstat[0]++;
     }

    double *sj = gsl_matrix_ptr (unitstat, 0, 0);
    double *pj = gsl_matrix_ptr (Punitstat, 0, 0);
      
    calcAdjustP(mmRef->punit, nVars, &buj.vector, sj, pj, sortid[0]);
    // calc significance stats
    for (i=1; i<nParam+1; i++) {    
       calcSS(bY, &(Hats[i+1]), mmRef, FALSE, FALSE, TRUE);
       gsl_vector_view buj=gsl_matrix_row (bUnitStat, i);
  //  displaymatrix (Hats[i+1].SS, "hypo.SS");
       testStatCalc(&(Hats[i+1]),&(Hats[0]),mmRef,FALSE,(bMultStat+i),&buj.vector);
       // count data related to P-values
       if ( bMultStat[i] >= multstat[i] )
          Pmultstat[i]++;
       sj = gsl_matrix_ptr (unitstat, i, 0);
       pj = gsl_matrix_ptr (Punitstat, i, 0);
       calcAdjustP(mmRef->punit, nVars, &buj.vector, sj, pj, sortid[i]);
    }
//   displaymatrix(Punitstat, "Punitstat");
  
   nSamp++;
   return 0;
}

int Summary::calcR2(void)
{
    gsl_matrix *mss=gsl_matrix_alloc (nVars, nVars);
    gsl_matrix *tss=gsl_matrix_alloc (nVars, nVars);
    gsl_matrix *YminusM=gsl_matrix_alloc(nRows, nVars);
    gsl_matrix *FminusM=gsl_matrix_alloc(nRows, nVars);


    // The interpreted SS (mss)
    gsl_vector *e=gsl_vector_alloc (nRows);
    gsl_vector_set_all (e, 1.0);

    // get fit
    gsl_matrix_memcpy (FminusM, Hats[0].Y);
    gsl_matrix_memcpy (YminusM, Yref);
    size_t i, j;
    double mean;
    for (j=0; j<nVars; j++) {
       // get mean(Y)
       gsl_vector_view mj=gsl_matrix_column (Yref, j);
       gsl_blas_ddot (&mj.vector, e, &mean);
       // fitMinusMeanY
       gsl_vector_view y0j=gsl_matrix_column (FminusM, j);
       gsl_vector_add_constant (&y0j.vector, -mean/nRows);
       //YminusMeanY
       gsl_vector_view y1j=gsl_matrix_column (YminusM, j);
       gsl_vector_add_constant (&y1j.vector, -mean/nRows); 
    }
    //mss = rcalc(fitMinusMean, COR_TYPE, lambda);
    rcalc(FminusM, mmRef->shrink_param, mmRef->corr, mss);
    //tss = rcalc(YMinusMean, COR_TYPE, lambda);
    rcalc(YminusM, mmRef->shrink_param, mmRef->corr, tss);

    gsl_vector_free(e);
    e=gsl_vector_alloc(nVars);
    gsl_vector_set_all (e, 1.0);
    // calc R squared
    if ( mmRef->rsquare == HOOPER ) {
       // fill mss;
       for (i=0; i<nVars; i++) {
          for (j=i+1; j<nVars; j++){
              gsl_matrix_set(mss, i, j, gsl_matrix_get(mss, j, i));
	      gsl_matrix_set(tss, i, j, gsl_matrix_get(tss, j, i));
           }
       }
       // tmp=inv(tss)*mss;
       int s;
       R2=0;
       gsl_permutation *p = gsl_permutation_alloc(nVars);
       gsl_linalg_LU_decomp (tss, p, &s);
       for (j=0; j<nVars; j++) {
           gsl_vector_view mssj=gsl_matrix_column (mss, j);
	   gsl_linalg_LU_svx (tss, p, &mssj.vector);
	   R2=R2+gsl_vector_get(&mssj.vector, j);          
       }
       R2 = R2/nVars;
       gsl_permutation_free(p);       
     }
    else if ( mmRef->rsquare == VECTOR )
       // R2 = det(mss)/det(tss); 
       R2 = (double)calcDet(mss)/calcDet(tss);      
    else 
       GSL_ERROR("Invalid R2 option", GSL_EINVAL);

    gsl_vector_free(e);
    gsl_matrix_free(mss);
    gsl_matrix_free(tss);
    gsl_matrix_free(YminusM);
    gsl_matrix_free(FminusM);

    return 0;
}
