// Resampling-based hypothesis test for comparing multivariate linear models
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 16-Nov-2009


#include "resampTest.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h> 
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_double.h>

// Use gsl functions as much as poosible to increase speed and stabilitiy

AnovaTest::AnovaTest(mv_Method *mm, gsl_matrix *Y, gsl_matrix *X, gsl_matrix *isXvarIn):mmRef(mm), Yref(Y), Xref(X), inRef(isXvarIn)
{
    int hid, aid;
    size_t i, j, count;
    nModels=inRef->size1, nParam=Xref->size2;
    nRows=Yref->size1, nVars=Yref->size2; 

//  printf("initialize public variables: stats\n");
    if ( nModels == 1 ){
        *multstat = 0.0;
	*Pmultstat = 0.0;
	*dfDiff = 0;
	statj = gsl_matrix_alloc(1, nVars);
	Pstatj = gsl_matrix_alloc(1, nVars);
    }	
    else { 
        multstat=(double *)malloc((nModels-1)*sizeof(double));
	Pmultstat = (double *)malloc((nModels-1)*sizeof(double));
        for (j=0; j<nModels-1; j++)
	    *(Pmultstat+j)=0.0; 
	dfDiff = (int *)malloc((nModels-1)*sizeof(int));

	statj = gsl_matrix_alloc(nModels-1, nVars);
	Pstatj = gsl_matrix_alloc(nModels-1, nVars);
	gsl_matrix_set_zero(Pstatj);
    }
    
//    printf("initialize private variables\n");
    bStatj = gsl_vector_alloc(nVars);
    Hats = (mv_mat *)malloc(nModels*sizeof(mv_mat)); 
    sortid = (gsl_permutation **)malloc((nModels-1)*sizeof(gsl_permutation *));
    for (i=0; i<nModels; i++ ) {
        // Hats[i]
        Hats[i].mat=gsl_matrix_alloc(nRows, nRows);
        Hats[i].SS=gsl_matrix_alloc(nVars, nVars);
        Hats[i].Res=gsl_matrix_alloc(nRows, nVars);
        Hats[i].Y = gsl_matrix_alloc(nRows, nVars);
	count = 0;
	for (j=0; j<nParam; j++){
	    count+=(int)gsl_matrix_get(inRef, i, j);
	}
//	printf("count=%d \n", count);
	Hats[i].X = gsl_matrix_alloc(nRows, count);
	Hats[i].Coef=gsl_matrix_alloc(count, nVars);
        gsl_vector_view refi=gsl_matrix_row(inRef, i);
	subX(Xref, &refi.vector, Hats[i].X);
        calcSS(Yref, &(Hats[i]), mmRef, TRUE, TRUE, TRUE);
 
//      printf(" - calc teststats\n");
        if (i>0){ 
           hid = i; aid = i-1;
           if ( mmRef->resamp != CASEBOOT ) {
              // fit = Y- resi 
              gsl_matrix_memcpy (Hats[i].Y, Yref);
              gsl_matrix_sub (Hats[i].Y, Hats[i].Res);
           }
           gsl_vector_view statij = gsl_matrix_row(statj, i-1);
           //displaymatrix(Hats[aid].SS, "SS1");

           testStatCalc(&(Hats[hid]), &(Hats[aid]), mmRef, TRUE, (multstat+i-1), &statij.vector); 
           dfDiff[aid] = Hats[aid].Coef->size1 - Hats[hid].Coef->size1; 
           //displayvector(&statij.vector, "statj");

           // sortid
           sortid[i-1] = gsl_permutation_alloc(nVars);
           gsl_sort_vector_index (sortid[i-1], &statij.vector); 
           // rearrange sortid in descending order
           gsl_permutation_reverse (sortid[i-1]);
       }
    }  

    // initialize resampling indices 
//    getBootID(); done in R
    bootID = NULL;
}

AnovaTest::~AnovaTest(){
}

void AnovaTest::releaseTest()
{
    free(multstat);
    free(Pmultstat);
    free(dfDiff);

    gsl_matrix_free(statj);
    gsl_matrix_free(Pstatj);
// The above commented out s.t. the results are passed out to R

    size_t i;
    for ( i=0; i<nModels; i++ ){
        gsl_matrix_free(Hats[i].mat);
        gsl_matrix_free(Hats[i].SS);
	gsl_matrix_free(Hats[i].Res);
        gsl_matrix_free(Hats[i].Coef);
	gsl_matrix_free(Hats[i].X);
	gsl_matrix_free(Hats[i].Y);
    }
    gsl_vector_free(bStatj);

    if ( bootID != NULL )
       gsl_matrix_free(bootID);

    for (i=0; i<nModels-1; i++)
       gsl_permutation_free(sortid[i]);
    free(sortid);

//    printf("Anova test released.\n");

}

void AnovaTest::display(void)
{
    size_t i, j, k;
    printf("Anova Table (resampling under H0):\n");
    printf("Hypo\t Alter\t df\t TestStat\t P-value\n");
    for ( i=0; i<nModels-1; i++ ){
       printf("Model%u\t Model%u\t %u\t %.3f\t\t %.3f\n", (unsigned int)i+1, (unsigned int)i, dfDiff[i], multstat[i], Pmultstat[i]);
    }

    if (mmRef->punit!=NONE){
       printf("Univariate Tests:\n");
       for (k=0; k<3; k++){
           printf("Response:\t");    
           for (j=k*4; j<(k+1)*4; j++)
               printf("variable%u\t", (unsigned int)j);    	
           printf("\n");
           for ( i=0; i<nModels-1; i++ ){       
               printf("Model %u:\t", (unsigned int)i+1);
               for (j=k*4; j<(k+1)*4; j++)
                   printf("%.3f(%.3f)\t", gsl_matrix_get(statj, i, j), gsl_matrix_get(Pstatj, i, j));
               printf("\n");
            }
	   printf("\n");
        }
     }	
}

int AnovaTest::resampTest(void)
{
//    printf("Start resampling test ...\n");
    size_t i, j, p, id;
    size_t maxiter=mmRef->nboot; 
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
       if ( mmRef->resamp == PERMUTE ){
          permid = (size_t *)malloc(nRows*sizeof(size_t));
          for (i=0; i<nRows; i++)
              permid[i] = i;
    } }
    // resampling options 
    if (mmRef->resamp == CASEBOOT) {
       nSamp = 0;
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
           anovacase(bY, bX);
           nSamp++;
        }
    } 
    else if (mmRef->resamp == RESIBOOT) {
        nSamp = 0;
        for (i=0; i<maxiter; i++) {
          for (p=1; p<nModels; p++) { 
            for (j=0; j<nRows; j++){
               // resampling index
 	       if (bootID == NULL) 
	          id = gsl_rng_uniform_int(rnd, nRows);
               else 
	          id = (size_t) gsl_matrix_get(bootID, i, j);
               // bootr by resampling resi=(Y-fit)
               gsl_vector_view Yj=gsl_matrix_row(Yref, id);
               gsl_vector_view Fj=gsl_matrix_row(Hats[p].Y, id);
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view bootr=gsl_matrix_row(bY, j);
               gsl_vector_sub (&bootr.vector, &Fj.vector);  
               if (mmRef->student==TRUE) {
                  hii = gsl_matrix_get(Hats[p].mat, id, id);
                  gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
               } 
               // bY = Y + bootr
               Yj=gsl_matrix_row(Hats[p].Y, j);
               gsl_vector_add (&bootr.vector, &Yj.vector);
	    } 
            anovaresi(bY, p);
         }
        nSamp++;
    } }
   else if (mmRef->resamp == SCOREBOOT) {
       nSamp = 0;
       for (i=0; i<maxiter; i++) {
         for (p=1; p<nModels; p++) {
           for ( j=0; j<nRows; j++ ) {
               // random score
	       if ( bootID == NULL )
	          score = gsl_ran_ugaussian (rnd); 
	       else
	          score = (double)gsl_matrix_get(bootID, i, j);
               // bootr = (Y - fit)*score 
               gsl_vector_view Yj=gsl_matrix_row(Yref, j);
               gsl_vector_view Fj=gsl_matrix_row(Hats[p].Y, j);
               gsl_matrix_set_row (bY, j, &Yj.vector);
               gsl_vector_view bootr=gsl_matrix_row(bY, j);
               gsl_vector_sub (&bootr.vector, &Fj.vector); 
               if (mmRef->student==TRUE) {
                  hii = gsl_matrix_get(Hats[p].mat, j, j);
                  gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
               }
                // bY = Y + bootr
               gsl_vector_scale (&bootr.vector, score);
               gsl_vector_add (&bootr.vector, &Fj.vector);
 	   } 
          anovaresi(bY, p);
        } 
        nSamp++;
   } }
   else if ( mmRef->resamp == PERMUTE ) { 
       gsl_matrix_add_constant (Pstatj, 1.0); 
       for (p=0; p<nModels-1; p++)
           Pmultstat[p]=1.0;       // include itself
        nSamp = 1;
        for (i=0; i<maxiter-1; i++) { //999
            for (p=1; p<nModels; p++){ 
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
                    gsl_vector_view Fj=gsl_matrix_row(Hats[p].Y, id);
                    gsl_matrix_set_row (bY, j, &Yj.vector);
                    gsl_vector_view bootr=gsl_matrix_row(bY, j);
                    gsl_vector_sub (&bootr.vector, &Fj.vector); 
                    if (mmRef->student==TRUE) {
                        hii = gsl_matrix_get(Hats[p].mat, id, id);
                        gsl_vector_scale (&bootr.vector, 1/sqrt(1-hii));
                    }
                    // bY = Y + bootr
                    Yj=gsl_matrix_row(Hats[p].Y, j);
                    gsl_vector_add (&bootr.vector, &Yj.vector);
                 }
                 anovaresi(bY, p);
           }
           nSamp++;
       }      
   }
   else 
       GSL_ERROR("Invalid resampling option", GSL_EINVAL);

   // p-values 
   size_t sid, sid0;
   double *pj;  
   for (i=0; i<nModels-1; i++) { 
        Pmultstat[i]=(double) Pmultstat[i]/nSamp;
        pj = gsl_matrix_ptr (Pstatj, i, 0);
        if ( mmRef->punit == FREESTEP ){ 
           for (j=1; j<nVars; j++){
               sid = gsl_permutation_get(sortid[i], j);
	       sid0 = gsl_permutation_get(sortid[i], j-1);
	       *(pj+sid)=MAX(*(pj+sid), *(pj+sid0)); 
	   }  
        }
        if ( mmRef->punit == STEPUP ){ 
           for (j=2; j<nVars; j++){
               sid = gsl_permutation_get(sortid[i], nVars-j);
	       sid0 = gsl_permutation_get(sortid[i], nVars-j+1);
	       *(pj+sid) = MIN(*(pj+sid), *(pj+sid0)); 
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

int AnovaTest::anovacase(gsl_matrix *bY, gsl_matrix *bX)
{
   size_t j;
   // if Y col is all zeros
   for ( j=0; j<nVars; j++ ){
       gsl_vector_view colj = gsl_matrix_column(bY, j);
       if ( gsl_vector_isnull(&colj.vector) == TRUE )
          return GSL_ERANGE;
   }

   size_t i, k, hid, aid;
   double *sj, *pj, *bj;
   // Hats.X 
   for (i=0; i<nModels; i++){
       k=0;
       for (j=0; j<nParam; j++)
       if (gsl_matrix_get(inRef, i, j)==TRUE){
	  gsl_vector_view Xj=gsl_matrix_column(bX, j);
          gsl_matrix_set_col(Hats[i].X, k, &Xj.vector);
	  k++;
       } 
       //Y = X*coef
       gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1.0,Hats[i].X,Hats[i].Coef,0.0,Hats[i].Y); 
       //Y = bY - Y;
       gsl_matrix_sub (Hats[i].Y, bY);
       gsl_matrix_scale (Hats[i].Y, -1.0);

       // calc teststats
       if (i>0) {
          hid = i; aid = i-1;  
          calcSS(Hats[aid].Y, &(Hats[hid]), mmRef, TRUE, FALSE, TRUE);
          calcSS(Hats[aid].Y, &(Hats[aid]), mmRef, TRUE, FALSE, TRUE);
          testStatCalc(&(Hats[hid]), &(Hats[aid]), mmRef, TRUE, &(bMultStat), bStatj);

          // count data related to P-values
          if (bMultStat >= multstat[i-1])
             Pmultstat[i-1]++;
          // get result ptr corresponds to model i
          sj = gsl_matrix_ptr (statj, i-1, 0);
	  pj = gsl_matrix_ptr (Pstatj, i-1, 0);
          bj = gsl_vector_ptr (bStatj, 0);          
          calcAdjustP(mmRef->punit, nVars, bj, sj, pj, sortid[i-1]);
      }
   }   

  return 0;
}

int AnovaTest::anovaresi(gsl_matrix *bY, const size_t i)
{
    size_t hid=i, aid = i-1;

    // count the right-hand tails
    calcSS(bY, &(Hats[aid]), mmRef, FALSE, FALSE, TRUE);
    calcSS(bY, &(Hats[hid]), mmRef, FALSE, FALSE, TRUE);
    testStatCalc(&(Hats[hid]), &(Hats[aid]), mmRef, TRUE, &(bMultStat), bStatj);
 
    // count data related to P-values
    if (bMultStat >= multstat[i-1])
        Pmultstat[i-1]++;
    // get result ptr corresponds to model i
    double *sj = gsl_matrix_ptr (statj, i-1, 0);
    double *pj = gsl_matrix_ptr (Pstatj, i-1, 0);
    double *bj = gsl_vector_ptr (bStatj, 0);
    calcAdjustP(mmRef->punit, nVars, bj, sj, pj, sortid[i-1]);    
       
   return 0;
}
