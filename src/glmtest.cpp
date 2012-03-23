// hypothesis testing, including summary and anova 
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// // 16-3-2012

#include "resampTest.h"
//#include "time.h"

GlmTest::GlmTest(const mv_Method *tm, glm *fit):tm(tm), fit(fit)
{  
    nRows = fit->nRows;
    nVars = fit->nVars;
    nParam = fit->nParams;

    smryStat = NULL;
    Psmry = NULL;

    anovaStat = NULL;
    Panova = NULL;

    Xin = NULL;
    GrpXs = NULL;
    GrpOs = NULL;

    XBeta = NULL;
    Sigma = NULL;
    bootID = NULL;

    // Prepared for geeCalc
    L = gsl_matrix_alloc(nParam, nParam);
    gsl_matrix_set_identity (L);
    Rlambda = gsl_matrix_alloc(nVars, nVars);
    bRlambda = gsl_matrix_alloc(nVars, nVars);
    Wj=gsl_matrix_alloc(nRows, nRows);

    rnd=gsl_rng_alloc(gsl_rng_mt19937);    
    if ( tm->resamp==PERMUTE ) {
        permid = (unsigned int *)malloc(nRows*sizeof(unsigned int));
	for (unsigned int i=0; i<nRows; i++) permid[i]=i;
    }	
    else permid=NULL;

    if (tm->resamp==MONTECARLO){
        XBeta = gsl_matrix_alloc(nRows, nVars);
        Sigma = gsl_matrix_alloc(nVars, nVars);
    }
    aic = new double [nVars];
}

GlmTest::~GlmTest(void){
}

void GlmTest::releaseTest(void)
{
    if (smryStat != NULL) gsl_matrix_free(smryStat);
    if (Psmry != NULL) gsl_matrix_free(Psmry);

    if (anovaStat != NULL) gsl_matrix_free(anovaStat);
    if (Panova != NULL) gsl_matrix_free(Panova);

    gsl_matrix_free(L);
    gsl_matrix_free(Rlambda);   
    gsl_matrix_free(bRlambda);   
    gsl_matrix_free(Wj);

    gsl_rng_free(rnd);
    if ( XBeta != NULL ) gsl_matrix_free(XBeta);
    if ( Sigma != NULL ) gsl_matrix_free(Sigma);

    if ( bootID != NULL ) gsl_matrix_free(bootID);
    if ( permid != NULL ) free(permid);

    delete[] aic;   

}

int GlmTest::summary(void)
{
    unsigned int k;
    unsigned int mtype = fit->mmRef->model-1;
    PoissonGlm pNull(fit->mmRef), pAlt(fit->mmRef);
    LogiGlm lNull(fit->mmRef), lAlt(fit->mmRef);
    NBinGlm nbNull(fit->mmRef), nbAlt(fit->mmRef);
    glm *PtrNull[3] = { &pNull, &nbNull, &lNull };
    glm *bAlt[3] = { &pAlt, &nbAlt, &lAlt };
    gsl_vector_view teststat, unitstat;

    smryStat = gsl_matrix_alloc(nParam+1, nVars+1);
    Psmry = gsl_matrix_alloc(nParam+1, nVars+1);
    gsl_matrix_set_zero (Psmry);

    // initialize the design matrix for all hypo tests
    GrpXs = (GrpMat *)malloc((nParam+2)*sizeof(GrpMat));
    GrpXs[0].matrix = gsl_matrix_alloc(nRows, nParam);
    gsl_matrix_memcpy(GrpXs[0].matrix, fit->Xref); // the alt X
    GrpXs[1].matrix = gsl_matrix_alloc(nRows, 1); // overall test
    gsl_matrix_set_all (GrpXs[1].matrix, 1.0);
    for (k=2; k<nParam+2; k++) { // significance tests
         GrpXs[k].matrix = gsl_matrix_alloc(nRows, nParam-1);
         subX2(fit->Xref, k-2, GrpXs[k].matrix);
    }
    GrpOs = (GrpMat *)malloc((nParam+2)*sizeof(GrpMat));
    for (k=0; k<nParam+2; k++) {
        GrpOs[k].matrix = gsl_matrix_alloc(nRows, nVars);
        gsl_matrix_set_zero (GrpOs[k].matrix);
    }
    // calc test statistics
    geeCalc(fit, PtrNull[mtype], smryStat, Rlambda); // for single test only
//    displaymatrix(smryStat, "smryStat");

    // sort id if the unitvaraite test is free step-down
    gsl_permutation **sortid;
    sortid=(gsl_permutation **)malloc((nParam+1)*sizeof(gsl_permutation *));
    for ( k=0; k<nParam+1; k++ ) {
        teststat = gsl_matrix_row (smryStat, k);
        unitstat=gsl_vector_subvector(&teststat.vector, 1, nVars);
        sortid[k] = gsl_permutation_alloc(nVars);
        gsl_sort_vector_index (sortid[k], &unitstat.vector);
        gsl_permutation_reverse(sortid[k]);  // rearrange in descending order
    }
   // ========= Get resampling distribution under H1 ====== //       
    if ( tm->resamp != FREEPERM ) {
        // enable offset for bootstrap    
        gsl_matrix_memcpy (GrpOs[0].matrix, fit->Eta);
        gsl_matrix_view X1=gsl_matrix_submatrix(fit->Xref,0,1,nRows,nParam-1);
        gsl_matrix_view B1=gsl_matrix_submatrix(fit->Beta,1,0,nParam-1,nVars);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,&X1.matrix,&B1.matrix,0,GrpOs[1].matrix);
        for (k=2; k<nParam+2; k++) {
            // offset = X(:,j)*Beta(j, :)
            gsl_vector_view xj=gsl_matrix_column (fit->Xref, k-2);
            gsl_vector_view bj=gsl_matrix_row (fit->Beta, k-2);
            gsl_blas_dger (1.0, &xj.vector, &bj.vector, GrpOs[k].matrix);
    }   }

    unsigned int nP;
    GrpMat *oriGrpXs=NULL, *oriGrpOs=NULL;
    if ( tm->resamp == CASEBOOT ) {
        oriGrpXs=(GrpMat *)malloc((nParam+2)*sizeof(GrpMat));
        oriGrpOs=(GrpMat *)malloc((nParam+2)*sizeof(GrpMat));
        for (k=0; k<nParam+2; k++) {
            nP = GrpXs[k].matrix->size2;
            oriGrpXs[k].matrix = gsl_matrix_alloc(nRows, nP);
            oriGrpOs[k].matrix = gsl_matrix_alloc(nRows, nVars);
            gsl_matrix_memcpy(oriGrpXs[k].matrix, GrpXs[k].matrix);
            gsl_matrix_memcpy(oriGrpOs[k].matrix, GrpOs[k].matrix);
    }   }

    if (tm->resamp==MONTECARLO) {
       if (tm->corr != SHRINK)
           GetR(fit->Res, SHRINK, tm->smry_lambda, 0, Rlambda);
       setMonteCarlo(fit, NULL, Rlambda);
    }

    nSamp=0;
    double *suj, *buj, *puj;
    gsl_matrix *bStat = gsl_matrix_alloc(nParam+1, nVars+1);
    gsl_matrix_set_zero (bStat);
    gsl_matrix *bY = gsl_matrix_alloc(nRows, nVars);
    for ( unsigned int i=0; i<tm->nboot; i++) {
        resampData(fit, bY, oriGrpXs, oriGrpOs, i);
        bAlt[mtype]->regression(bY, GrpXs[0].matrix, GrpOs[0].matrix);
        geeCalc(bAlt[mtype], PtrNull[mtype], bStat, bRlambda);
        for (k=0; k<nParam+1; k++) {
           buj = gsl_matrix_ptr (bStat, k, 0);
           suj = gsl_matrix_ptr (smryStat, k, 0);
           puj = gsl_matrix_ptr (Psmry, k, 0);
           if ( *buj >= *suj ) *puj=*puj+1;
           calcAdjustP(tm->punit, nVars, buj+1, suj+1, puj+1, sortid[k]);
        } // end for j loop
        nSamp++;
    } // end for i loop

    // ========= Get P-values ========= //        
    if ( tm->punit == FREESTEP ) {
       for (k=0; k<nParam+1; k++) {
           puj = gsl_matrix_ptr (Psmry, k, 1);
           reinforceP( puj, nVars, sortid[k] );
    }  }
    // p = (#exceeding observed stat + 1)/(#nboot+1)
    gsl_matrix_add_constant (Psmry, 1.0);
    gsl_matrix_scale (Psmry, (double)1.0/(nSamp+1));
    for (k=0; k<nVars; k++) aic[k]=-fit->ll[k]+2*(nParam+1);

    // === release memory ==== //
    bAlt[mtype]->releaseGlm();
    PtrNull[mtype]->releaseGlm();
    gsl_matrix_free(bStat);
    gsl_matrix_free(bY);

    for (k=0; k<nParam+1; k++) 
       if (sortid[k]!=NULL) gsl_permutation_free(sortid[k]);
    free(sortid);

    if ( GrpXs != NULL ) {
       for ( unsigned int k=0; k<nParam+2; k++ )
           if ( GrpXs[k].matrix != NULL )
              gsl_matrix_free (GrpXs[k].matrix);
       free(GrpXs);
    }
    if ( GrpOs != NULL ) {
       for ( unsigned int k=0; k<nParam+2; k++ )
           if ( GrpOs[k].matrix != NULL )
              gsl_matrix_free (GrpOs[k].matrix);
       free(GrpOs);
    }
    if ( oriGrpXs != NULL ) {
       for ( unsigned int k=0; k<nParam+2; k++ )
           if ( oriGrpXs[k].matrix != NULL )
              gsl_matrix_free (oriGrpXs[k].matrix);
       free(oriGrpXs);
    }
    if ( oriGrpOs != NULL ) {
       for ( unsigned int k=0; k<nParam+2; k++ )
           if ( oriGrpOs[k].matrix != NULL )
              gsl_matrix_free (oriGrpOs[k].matrix);
       free(oriGrpOs);
    }

    return SUCCESS;
}


int GlmTest::anova(gsl_matrix *isXvarIn) 
{
    // Assume the models have been already sorted (in R)
    Xin = isXvarIn;
    nModels = Xin->size1;
    double *rdf = new double [nModels];

    unsigned int nP, i, j, k;
    unsigned int ID0, ID1, nP0, nP1;

    dfDiff = new unsigned int [nModels-1];
    anovaStat = gsl_matrix_alloc(nModels-1, nVars+1);
    Panova = gsl_matrix_alloc(nModels-1, nVars+1);
    gsl_vector *bStat = gsl_vector_alloc(nVars+1);
    gsl_matrix_set_zero (anovaStat);    
    gsl_matrix_set_zero (Panova);
    gsl_vector_set_zero (bStat);

    PoissonGlm pNull(fit->mmRef), pAlt(fit->mmRef);
    LogiGlm lNull(fit->mmRef), lAlt(fit->mmRef);
    NBinGlm nbNull(fit->mmRef), nbAlt(fit->mmRef);
    PoissonGlm pNullb(fit->mmRef), pAltb(fit->mmRef);
    LogiGlm lNullb(fit->mmRef), lAltb(fit->mmRef);
    NBinGlm nbNullb(fit->mmRef), nbAltb(fit->mmRef);
    glm *PtrNull[3] = { &pNull, &nbNull, &lNull };
    glm *PtrAlt[3] = { &pAlt, &nbAlt, &lAlt };
    glm *bNull[3] = { &pNullb, &nbNullb, &lNullb };
    glm *bAlt[3] = { &pAltb, &nbAltb, &lAltb };

    gsl_permutation *sortid=NULL;
    if (tm->punit==FREESTEP) sortid = gsl_permutation_alloc(nVars);

    double *suj, *buj, *puj;
    gsl_vector_view teststat, unitstat, ref1, ref0; 
    gsl_matrix *Xnull=NULL, *L1=NULL, *tmp1=NULL;
    gsl_matrix *bX=NULL, *bXAlt=NULL, *bXnull=NULL, *bO=NULL;
    gsl_matrix *bY=gsl_matrix_alloc(nRows, nVars);
    gsl_matrix_view X1, X0;

    unsigned int mtype = fit->mmRef->model-1;
    // ======= Fit the (first) Alt model =========//
    for (i=0; i<nModels; i++) {
        nP = 0;
        for (k=0; k<nParam; k++) 
	     if (gsl_matrix_get(Xin,i,k)!=FALSE) nP++;   
        rdf[i] = nRows-nP;
    }
//    nP1 = nRows-(unsigned int)rdf[0];
//    X1=gsl_matrix_submatrix(fit->Xref,0,0,nRows,nP1);
//    PtrAlt[mtype]->EstIRLS(fit->Yref, &X1.matrix, NULL, fit->phi); 
    
     PtrAlt[mtype]->copyGlm(fit);

     for (i=1; i<nModels; i++) {       
        // ======= Fit the Null model =========//
        ID0 = i; ID1 = i-1;
        nP0 = nRows - (unsigned int)rdf[ID0];
        nP1 = nRows - (unsigned int)rdf[ID1];

        ref1=gsl_matrix_row(Xin, ID1);
        ref0=gsl_matrix_row(Xin, ID0);
//      Xnull = gsl_matrix_alloc(nRows, nP0);
//      subX(fit->Xref, &ref0.vector, Xnull);
//	PtrNull[mtype]->EstIRLS(fit->Yref, &X0.matrix, NULL, PtrAlt[mtype]->phi);   
        X0=gsl_matrix_submatrix(fit->Xref,0,0,nRows,nP0);
	PtrNull[mtype]->regression(fit->Yref, &X0.matrix, NULL);   

        GetR(PtrAlt[mtype]->Res,tm->corr,tm->anova_lambda,ID1,Rlambda); 

	// ======= Get multivariate test statistics =======//
        // Estimate shrinkage parametr only once under H1 
        // See "FW: Doubts R package "mvabund" (12/14/11)
        teststat = gsl_matrix_row(anovaStat, i-1);
        // L1 encodes info between X0 and X1
        L1 = gsl_matrix_alloc (nP1-nP0, nP1);
        tmp1 = gsl_matrix_alloc (nParam, nP1);
        subX(L, &ref1.vector, tmp1);
        subXrow1(tmp1, &ref0.vector, &ref1.vector, L1);
      
        if (tm->test == WALD) 
           subGeeWald(PtrAlt[mtype], L1, &teststat.vector, Rlambda);
        else if (tm->test == SCORE) 
           subGeeScore(PtrAlt[mtype],PtrNull[mtype],&teststat.vector,Rlambda);
        else subGeeLR(PtrAlt[mtype], PtrNull[mtype], &teststat.vector);	
// Rprintf("LR=%.2f\n", gsl_vector_get(&teststat.vector, 0));

	// ======= Get univariate test statistics =======//
        if (tm->punit == FREESTEP) {	
              unitstat=gsl_vector_subvector(&teststat.vector,1,nVars);
              gsl_sort_vector_index (sortid, &unitstat.vector);
              gsl_permutation_reverse(sortid);        
        }
        // ======= Prepare for resampling methods ======//
	if ( tm->resamp == CASEBOOT ) {
            bX = gsl_matrix_alloc(nRows, nParam);
            bXAlt = gsl_matrix_alloc(nRows, nP1); 
            bXnull = gsl_matrix_alloc(nRows, nP0); 
            bO = gsl_matrix_alloc(nRows, nVars);
	}
        else bX = PtrAlt[mtype]->Xref;
       
	if (tm->resamp == MONTECARLO)  {
           if (tm->corr != SHRINK)  
              GetR(PtrAlt[mtype]->Res, SHRINK, tm->anova_lambda, ID1, Rlambda);
           setMonteCarlo (PtrNull[mtype], NULL, Rlambda);
        }

        // ======= Get resampling distribution under H0 ===== //
	nSamp=0;
//        gsl_matrix_view Beta0;
        for (j=0; j<tm->nboot; j++) {	
//        for (j=0; j<1; j++) {	
//            printf("j=%d \n ", j);
	    gsl_vector_set_zero (bStat);
	    if ( tm->resamp == CASEBOOT ) {
                resampAnovaCase(fit, bY, bX, j);
                subX(bX, &ref1.vector, bXAlt);
                gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, bXAlt, PtrAlt[mtype]->Beta, 0.0, bO);
	        bAlt[mtype]->regression(bY, bXAlt, bO);
//              bAlt[mtype]->display();
                if (tm->test != WALD){
                   subX(bX, &ref0.vector, bXnull);
//                   Beta0 = gsl_matrix_submatrix(PtrAlt[mtype]->Beta, 0, 0, nP0, nVars);
                   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, bXnull, PtrNull[mtype]->Beta, 0.0, bO);
                   bNull[mtype]->EstIRLS(bY,bXnull,bO,bAlt[mtype]->phi);
//                 bNull[mtype]->display()
                }    
           }		
           else {
		resampNonCase(PtrNull[mtype], bY, j);
                bAlt[mtype]->regression(bY, PtrAlt[mtype]->Xref, NULL ); 
//                bAlt[mtype]->EstIRLS(bY,PtrAlt[mtype]->Xref,NULL,PtrAlt[mtype]->phi);
                if (tm->test != WALD) 
                    bNull[mtype]->EstIRLS(bY,PtrNull[mtype]->Xref,NULL,bAlt[mtype]->phi);
                    bNull[mtype]->regression(bY, PtrNull[mtype]->Xref, NULL ); 
	   }

           GetR(bAlt[mtype]->Res,tm->corr,tm->anova_lambda,ID1,bRlambda);   

           if ( tm->test == WALD ) 
               subGeeWald(bAlt[mtype], L1, bStat, bRlambda);
           else if ( tm->test == SCORE ) 
               subGeeScore(bAlt[mtype], bNull[mtype], bStat, bRlambda);
           else 
               subGeeLR(bAlt[mtype], bNull[mtype], bStat);	
        
          // ----- get multivariate counts ------- //	  
          buj = gsl_vector_ptr (bStat,0);
          suj = gsl_matrix_ptr (anovaStat, i-1, 0);
          puj = gsl_matrix_ptr (Panova, i-1, 0);
          if ( *(buj) >= *(suj) ) *puj=*puj+1;
          // ------ get univariate counts ---------//	   	   
          calcAdjustP(tm->punit,nVars,buj+1,suj+1,puj+1,sortid);
	  nSamp++;
      } // end j for loop

      // ========= get p-values ======== //
      if ( tm->punit == FREESTEP) {
         puj = gsl_matrix_ptr (Panova, i-1, 1);
         reinforceP(puj, nVars, sortid);
      }

      // Degrees of freedom
      rdf[i] = PtrNull[mtype]->rdf;
      dfDiff[i-1] = PtrNull[mtype]->rdf - PtrAlt[mtype]->rdf;

      PtrAlt[mtype]->copyGlm(PtrNull[mtype]);

      if ( tm->resamp == CASEBOOT ) {
         if (bX!=NULL) gsl_matrix_free(bX);   
         if (bXAlt!=NULL) gsl_matrix_free(bXAlt);   
         if (bXnull!=NULL) gsl_matrix_free(bXnull);   
         if (bO!=NULL) gsl_matrix_free(bO);   
      }
      if (Xnull!=NULL) gsl_matrix_free(Xnull);   
      if (L1!=NULL) gsl_matrix_free(L1);
      if (tmp1!=NULL) gsl_matrix_free(tmp1);

    } // end i for loop  and test for loop

    // p = (#exceeding observed stat + 1)/(#nboot+1)
    gsl_matrix_add_constant (Panova, 1.0);
    gsl_matrix_scale (Panova, (double)1/(nSamp+1.0));

    PtrAlt[mtype]->releaseGlm();
    PtrNull[mtype]->releaseGlm();
    bAlt[mtype]->releaseGlm();
    if ( tm->test != WALD )
        bNull[mtype]->releaseGlm();
    delete []rdf;
    if (sortid != NULL )
        gsl_permutation_free(sortid);
    gsl_vector_free(bStat);
    gsl_matrix_free(bY);   
    
    return SUCCESS;
}
/*
void GlmTest::displaySmry(void)
{
    unsigned int i, j, k, nk;
    const char *testname[3] // string array, only pointers stored
               ={ "sqrt(WALD)", "SCORE", "LR" }; // 2, 3, 4

    printf("\nSummary of fitting (resampling under H1):\n");
    printf("\n - Regression performance: \n");
    printf("\t\t phi\t AIC\t log-like\t");
    for (j=0; j<nVars; j++){
        printf("\n[Response %d]:\t", (unsigned int)j+1);
	printf("%.4f\t ", fit->phi[j]);
	printf("%.2f\t ", aic[j]);
	printf("%.2f", fit->ll[j]);
    }
       
    if ( tm->corr == SHRINK ) 
       displayvector(tm->smry_lambda, "\n Est. shrink.param in summary\n");

    // significance test
    nk = 1;
    k = tm->test - 2;
    printf("\n - Significance test (Pr>=%s):\n", testname[k]);
    if (tm->punit==FREESTEP) printf("\t (FREESTEP adjusted)\n");
    while ( nk < nParam+1 ) {       
         printf("\t");
         for (i=nk; i<MIN(nk+4, nParam+1); i++) 
             printf("\t [Explain %d] ", (unsigned int)i);
         printf("\n\t ");    
         for (i=nk; i<MIN(nk+4, nParam+1); i++) 
             printf(" %.3f (%.3f) \t", gsl_matrix_get(smryStat,i,0), gsl_matrix_get(Psmry,i,0));
         printf("\n\n");       
         // Significance univariate tests
         if (tm->punit>NONE) {
            for (j=1; j<nVars+1; j++) {
                printf("[Response %d]:\t", (int)j);
                for (i=nk; i<MIN(nk+4, nParam+1); i++) 
                    printf("%.3f (%.3f)\t", gsl_matrix_get(smryStat,i,j), gsl_matrix_get(Psmry,i,j));
                printf("\n");               
            }
         }
         nk = i;
         printf("\n");
     }
     // Overall statistics
     printf("\n - Multivariate test (Pr>=%s): %.3f (%.3f)",testname[k], gsl_matrix_get(smryStat, 0, 0), gsl_matrix_get(Psmry, 0, 0));
     if (tm->punit == FREESTEP) {
        printf("\t (FREESTEP adjusted)\n");
        for (j=1; j<nVars+1; j++) printf( "[Response %d]:\t%.3f (%.3f)\n", (int)j, gsl_matrix_get(smryStat, 0, j), gsl_matrix_get(Psmry,0,j));             
     }
     printf("\n ========================= \n");
    
}

void GlmTest::displayAnova(void)
{
    if (nModels==1) displaySmry();
    else {   
       unsigned int i, j;
       const char *testname[3] // string array, only pointers stored
               ={ "sqrt(WALD)", "SCORE", "LR" }; // 2, 3, 4
    
       printf("\n ========================= \n");
       printf("\nAnova Table (resampling under ");
       if (tm->resamp==CASEBOOT) printf("H1):\n");
       else printf("H0):\n");

       if ( tm->corr == SHRINK ) 
          displayvector(tm->anova_lambda, "Est. shrink.param in anova");

       unsigned int test=tm->test-2;
       printf("Hypo\t Alter\t dff\t %s\t  P-value \n", testname[test]);
       for ( i=0; i<nModels-1; i++ )
           printf("M%d\t M%d\t %d\t %.3f   %.3f\t\t \n",(int)i+1,(int)i,dfDiff[i],gsl_matrix_get(anovaStat, i, 0),gsl_matrix_get(Panova, i, 0));

       if (tm->punit != NONE) {
           if (tm->punit == FREESTEP)
              printf("\nUnivariate Tests (FREESTEP adjusted):\n\t\t");
           else printf("\nUnivariate Tests:\n\t\t");
           for (i=0; i<nModels-1; i++) printf("\tM%d v. M%d\t", (unsigned int)i+1, (unsigned int)i);
           printf("\n");

           for (j=1; j<nVars+1; j++) {
               printf("[Response %d]:", (unsigned int)j);
               for (i=0; i<nModels-1; i++)
                   printf("\t%.3f (%.3f)", gsl_matrix_get(anovaStat,i,j), gsl_matrix_get(Panova,i,j));
               printf("\n");
           }
           printf("\n");
       }
    }
}
*/
int GlmTest::geeCalc(glm *PtrAlt, glm *PtrNull, gsl_matrix *Stats, gsl_matrix *R)
{
    if ( tm->test == WALD ) geeWald(PtrAlt, Stats, R);
    else if (tm->test == SCORE) geeScore(PtrAlt, PtrNull, Stats, R);
    else geeLR(PtrAlt, PtrNull, Stats);
    return SUCCESS;
}



int GlmTest::geeWald(glm *PtrAlt, gsl_matrix *Stats, gsl_matrix *R)
{
   GetR(PtrAlt->Res, tm->corr, tm->smry_lambda, 0, R);  
   // The overall test    
   gsl_matrix_view L1 = gsl_matrix_submatrix(L, 1, 0, nParam-1, nParam);
   gsl_vector_view teststat=gsl_matrix_row(Stats, 0);	
   subGeeWald(PtrAlt, &L1.matrix, &teststat.vector, R);
   // The significance test    
   for (unsigned int k=1; k<nParam+1; k++) {    
       L1 = gsl_matrix_submatrix(L, k-1, 0, 1, nParam);
       teststat=gsl_matrix_row(Stats, k);
       subGeeWald(PtrAlt, &L1.matrix, &teststat.vector, R);
   }  
   return SUCCESS;
}   

int GlmTest::geeScore(glm *PtrAlt, glm *PtrNull, gsl_matrix *Stats, gsl_matrix *R)
{
    for (unsigned int k=1; k<nParam+2; k++) {
        PtrNull->EstIRLS(PtrAlt->Yref,GrpXs[k].matrix,GrpOs[k].matrix,PtrAlt->phi);
        GetR(PtrNull->Res, tm->corr, tm->smry_lambda, k, R);  
        gsl_vector_view teststat=gsl_matrix_row(Stats, k-1);
        subGeeScore(PtrAlt, PtrNull, &teststat.vector, R);
    }     	        
   return SUCCESS;
}

int GlmTest::geeLR(glm *PtrAlt, glm *PtrNull, gsl_matrix *Stats)
{
    for (unsigned int k=0; k<nParam+1; k++) {
        PtrNull->EstIRLS(PtrAlt->Yref,GrpXs[k+1].matrix,GrpOs[k+1].matrix,PtrAlt->phi);
        gsl_vector_view teststat=gsl_matrix_row(Stats, k);
        subGeeLR(PtrAlt, PtrNull, &teststat.vector);
    }     	        
    return SUCCESS;
}

int GlmTest::subGeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat)
{
    double val, result=0;
    for ( unsigned int j=0; j<nVars; j++ ) { // univariates
        val = PtrAlt->ll[j] - PtrNull->ll[j];
        gsl_vector_set(teststat, j+1, val);
        result = result+val;	     
    }
    gsl_vector_set(teststat, 0, result); // multivariate
    return SUCCESS;

}

int GlmTest::subGeeScore(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat, gsl_matrix *R)
{
    double result, alpha;
    unsigned int j, l, nP = PtrAlt->Xref->size2;

    gsl_vector *U = gsl_vector_alloc(nVars*nP);
    gsl_matrix *kRlNull = gsl_matrix_alloc(nVars*nP, nVars*nP);
    gsl_matrix_set_zero (kRlNull);
    gsl_matrix *XwX = gsl_matrix_alloc(nP, nP);
    gsl_vector *tmp=gsl_vector_alloc(nVars*nP);
 
    GrpMat *Z = (GrpMat*)malloc(nVars*sizeof(GrpMat)); 
    for (j=0; j<nVars; j++) {
        Z[j].matrix = gsl_matrix_alloc(nRows, nP);
        gsl_matrix_memcpy(Z[j].matrix, PtrAlt->Xref);	 
//        gsl_matrix_memcpy(Z[j].matrix, PtrNull->Xref);	 

        gsl_matrix_set_zero(Wj);
        gsl_vector_view dw=gsl_matrix_diagonal(Wj);
        gsl_vector_view wj=gsl_matrix_column(PtrNull->wHalf, j);
        gsl_vector_memcpy(&dw.vector, &wj.vector);
        gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1,Wj,Z[j].matrix); // W^1/2 * X
        gsl_vector_view uj=gsl_vector_subvector(U, j*nP, nP); 
        gsl_vector_view rj=gsl_matrix_column(PtrNull->Res, j);
        gsl_blas_dgemv(CblasTrans, 1, Z[j].matrix, &rj.vector, 0, &uj.vector);

        // Get X^T*W*X
        if ( (tm->punit>0) || (tm->corr==IDENTITY) ) {
            gsl_matrix_set_zero(XwX);
            gsl_blas_dsyrk(CblasLower, CblasTrans, 1, Z[j].matrix, 0, XwX);
        }
        if ( tm->corr==IDENTITY ) { // kappaRlNull = X^T*W*X 
           gsl_matrix_view Rl=gsl_matrix_submatrix(kRlNull,j*nP,j*nP,nP,nP);
           gsl_matrix_memcpy (&Rl.matrix, XwX); 
        }
        else { // kappaRlNull = R(j, l)*W0X{j}^T*W0X{l}                 
             for (l=0; l<=j; l++) { // lower half
                 alpha = gsl_matrix_get(R, j, l);
                 gsl_matrix_view Rl=gsl_matrix_submatrix(kRlNull,j*nP,l*nP,nP,nP);
                 gsl_blas_dgemm(CblasTrans, CblasNoTrans, alpha, Z[j].matrix, Z[l].matrix, 0, &Rl.matrix);  
            }   
        }      
        // Get univariate teststat
        if ( (tm->punit>0) ) { // inv(vNaive)
	   gsl_vector_view tmp2=gsl_vector_subvector(tmp, 0, nP);             
           gsl_linalg_cholesky_decomp(XwX); 
           gsl_linalg_cholesky_solve(XwX, &uj.vector, &tmp2.vector); 
           gsl_blas_ddot(&uj.vector, &tmp2.vector, &result);
           gsl_vector_set(teststat, j+1, result);
	}
    } // end for j=1:nVars
      
    // multivariate test stat	
    // kappaRl \ Us(:)
    // Note: kRlNull may be singular because of the correlation matrix R 
    // GSL cholesky can only handle non-singular matrix, so use QR instead
        gsl_linalg_cholesky_decomp (kRlNull); 
        gsl_linalg_cholesky_solve (kRlNull, U, tmp);  
//    invLSQ(kRlNull, U, tmp);

    gsl_blas_ddot (U, tmp, &result);
    gsl_vector_set(teststat, 0, result);
          
    gsl_vector_free(U);
    gsl_vector_free(tmp);
    gsl_matrix_free(XwX);
    gsl_matrix_free(kRlNull);
    for (j=0; j<nVars; j++) gsl_matrix_free(Z[j].matrix);
    free(Z);    

    return SUCCESS;

}

// Wald Test used in both summary and anova (polymophism)
int GlmTest::subGeeWald(glm *Alt, gsl_matrix *LL, gsl_vector *teststat, gsl_matrix *R)
{
    unsigned int i, j, l;
    double alpha, result;
    unsigned int nP = Alt->nParams;
    unsigned int nDF = LL->size1;    

    gsl_vector *LBeta = gsl_vector_alloc(nVars*nDF);	 
    gsl_vector_set_zero(LBeta);
    gsl_matrix *w1jX1=gsl_matrix_alloc(nRows, nP);
    gsl_matrix *XwX=gsl_matrix_alloc(nP, nP);
    gsl_matrix *Rl2 = gsl_matrix_alloc(nP, nDF);
    gsl_matrix *IinvN = gsl_matrix_alloc(nDF, nDF);
    gsl_matrix *IinvRl = gsl_matrix_alloc(nVars*nDF, nVars*nDF);
    gsl_vector *tmp = gsl_vector_alloc(nVars*nDF);
    gsl_vector_view tmp2, dw, wj, wXj, zj, LBj, bj, li, rl;

    gsl_matrix_set_zero(IinvRl);
    GrpMat *Z = (GrpMat*)malloc(nVars*sizeof(GrpMat)); 
    for (j=0; j<nVars; j++){
       Z[j].matrix = gsl_matrix_alloc(nP, nRows);
       gsl_matrix_memcpy(w1jX1, Alt->Xref);	 
       gsl_matrix_set_zero(Wj);
       dw=gsl_matrix_diagonal(Wj);
       wj=gsl_matrix_column(Alt->wHalf, j);
       gsl_vector_memcpy(&dw.vector, &wj.vector);
       // W^1/2 * X
       gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1,Wj,w1jX1);
       // X^T WX
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, w1jX1, 0.0, XwX);
       gsl_linalg_cholesky_decomp (XwX);  // note XwX has to be positive def
       // Z = (X^T W X)^-1 * X^T W^1/2. 
       for (i=0; i<nRows; i++) {
           wXj=gsl_matrix_row (w1jX1, i);
           zj=gsl_matrix_column(Z[j].matrix, i );
           gsl_linalg_cholesky_solve (XwX, &wXj.vector, &zj.vector);
       }	        
       // LBeta = L*Beta       
       LBj=gsl_vector_subvector (LBeta, j*nDF, nDF);
       bj=gsl_matrix_column (Alt->Beta, j);
       gsl_blas_dgemv(CblasNoTrans,1,LL,&bj.vector,0,&LBj.vector);

       // Get Naive estimator for IinvN = L*inv(X'*W*X)*L'
       if ( (tm->punit>0) || (tm->corr==IDENTITY) ){
          gsl_matrix_set_zero(IinvN);
          gsl_matrix_set_zero(Rl2);
          // Rl2 = (X^T W X)^-1 * L^T
          for (i=0; i<nDF; i++){
              li=gsl_matrix_row (LL, i);
              rl=gsl_matrix_column (Rl2, i);
              gsl_linalg_cholesky_solve(XwX, &li.vector, &rl.vector);
          }
          // IinvN = L*inv(vNaive)*L^T = L * Rl2 
          gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LL,Rl2,0.0,IinvN);
       }

       // Get multivariate teststat
       if ( tm->corr==IDENTITY ) { // Diag(IinvRl) = IinvN
          gsl_matrix_view Rl=gsl_matrix_submatrix(IinvRl,j*nDF,j*nDF,nDF,nDF);
          gsl_matrix_memcpy(&Rl.matrix, IinvN);
       }
       else { // Get IinvRl=L*vSandRl*L^T 
          for (l=0; l<=j; l++) {
              gsl_matrix_set_zero(XwX);
              gsl_matrix_set_zero(Rl2);
              alpha = gsl_matrix_get(R, j, l);
              gsl_matrix_view Rl=gsl_matrix_submatrix(IinvRl,j*nDF,l*nDF,nDF,nDF);
              gsl_blas_dgemm(CblasNoTrans,CblasTrans,alpha,Z[j].matrix,Z[l].matrix, 0.0, XwX);  // borrow XwX to accommodate vSandRl
//	      displaymatrix(XwX, "vSandRl");
              gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,XwX,LL,0.0,Rl2); // borrow Rl2 space to accommodate vSandRl*L^T 
              gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LL,Rl2,0.0,&Rl.matrix); // L*vSandRl*L^T
          } // end l
       }  // end if (tm->corr) 

       // Get univariate teststat: statj=LBeta'*inv(IinvN)*LBeta
       if ( tm->punit>0 ){
          tmp2=gsl_vector_subvector(tmp, 0, nDF);
          invLSQ(IinvN, &LBj.vector, &tmp2.vector);
          gsl_blas_ddot (&LBj.vector, &tmp2.vector, &result);
          gsl_vector_set(teststat, j+1, sqrt(result));
       }
    } // end for j=1:nVars	 

    // Get multivariate test stat = LBeta^T * inv(IinvRl) * LBeta
// Note: kRlNull may be singular because of the correlation matrix R 
// // GSL cholesky can only handle non-singular matrix, so use QR instead
    gsl_linalg_cholesky_decomp (IinvRl); 
    gsl_linalg_cholesky_solve (IinvRl, LBeta, tmp); 
//    invLSQ(IinvRl, LBeta, tmp);

    gsl_blas_ddot (LBeta, tmp, &result);
    gsl_vector_set(teststat, 0, sqrt(result));	  
   
    for (j=0; j<nVars; j++) gsl_matrix_free(Z[j].matrix);
    free(Z);    
    gsl_vector_free(LBeta);
    gsl_matrix_free(w1jX1);
    gsl_matrix_free(XwX);
    gsl_matrix_free(Rl2);
    gsl_matrix_free(IinvN);
    gsl_matrix_free(IinvRl);
    gsl_vector_free(tmp);

    return SUCCESS;
}

int GlmTest::resampData(glm *model, gsl_matrix *bT, GrpMat *oriXs, GrpMat *oriOs, unsigned int i )
{
    unsigned int j, k, id;
    gsl_vector_view yj, oj;
//  gsl_vector_view xj;
    
    if ( tm->resamp == CASEBOOT ){        
        for (j=0; j<nRows; j++) {
           if (bootID == NULL) id = gsl_rng_uniform_int(rnd, nRows);
	   else id = (unsigned int) gsl_matrix_get(bootID, i, j);
	   // resample Y, X, offsets accordingly
           yj = gsl_matrix_row(model->Yref, id);
           gsl_matrix_set_row(bT, j, &yj.vector);
	   for (k=0; k<nParam+2; k++) {
//               xj = gsl_matrix_row(oriXs[k].matrix, j);
//               gsl_matrix_set_row(GrpXs[k].matrix, j, &xj.vector);          
               oj = gsl_matrix_row(oriOs[k].matrix, id);
 	       gsl_matrix_set_row(GrpOs[k].matrix, j, &oj.vector);
    	   }   
        }
    }
    else resampNonCase(model, bT, i);       

    return SUCCESS;
}

int GlmTest::resampAnovaCase(glm *model, gsl_matrix *bT, gsl_matrix *bX, unsigned int i)
{
    unsigned int j, id, isSingular;
    gsl_vector_view yj;
    gsl_vector_view xj;
    double det;
    unsigned int nParams = model->Xref->size2;
    gsl_matrix *txX = gsl_matrix_alloc(nParams, nParams);
    gsl_matrix_set_zero(txX);

    if (bootID == NULL) {
       isSingular=TRUE;
       while (isSingular==TRUE) {
            for (j=0; j<nRows; j++) {   
                id = gsl_rng_uniform_int(rnd, nRows);
                // resample Y and X and offset
                yj=gsl_matrix_row(model->Yref, id);
                gsl_matrix_set_row (bT, j, &yj.vector);
                xj = gsl_matrix_row(model->Xref, id);
                gsl_matrix_set_row(bX, j, &xj.vector);
             }
             gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, bX, 0.0, txX);
             det = calcDet(txX); 
       //      printf("det=%.4f\n", det);
             if (det>1e-8) isSingular=FALSE;
       } 
   }   		    	
   else {
       for (j=0; j<nRows; j++) {   
          id = (unsigned int) gsl_matrix_get(bootID, i, j);
          // resample Y and X and offset
          yj=gsl_matrix_row(model->Yref, id);
          gsl_matrix_set_row (bT, j, &yj.vector);
          xj = gsl_matrix_row(model->Xref, id);
          gsl_matrix_set_row(bX, j, &xj.vector);
       }
   }

   gsl_matrix_free(txX);

   return SUCCESS;
} 


int GlmTest::resampNonCase(glm *model, gsl_matrix *bT, unsigned int i)
{
   unsigned int j, k, id;
   double bt, score, yij, eij;
   gsl_vector_view yj;

   // note that residuals have got means subtracted
   switch (tm->resamp) {
   case RESIBOOT: 
       for (j=0; j<nRows; j++) {
           if (bootID==NULL) id = gsl_rng_uniform_int(rnd, nRows);
           else id = (unsigned int) gsl_matrix_get(bootID, i, j);
           // bY = mu+(bootr*sqrt(variance))
           for (k=0; k<nVars; k++) { 
               bt=gsl_matrix_get(model->Mu,j,k)+sqrt(gsl_matrix_get(model->Var,j,k))*gsl_matrix_get(model->Res, id, k)*gsl_matrix_get(model->sqrt1_Hii, j, k);  
               bt = MAX(bt, 0.0);
               if (model->mmRef->model == LOGIT) bt = MIN(bt, 1.0);
               gsl_matrix_set(bT, j, k, bt);
        }   }   	  	
        break;
   case SCOREBOOT: 
        for (j=0; j<nRows; j++) {
            if (bootID==NULL) score = gsl_ran_ugaussian (rnd);
            else score = (double)gsl_matrix_get(bootID, i, j);
            // bY = mu + score*sqrt(variance)  
	    for (k=0; k<nVars; k++){
                bt=gsl_matrix_get(model->Mu, j, k)+sqrt(gsl_matrix_get(model->Var, j, k))*gsl_matrix_get(model->Res, j, k)*score*gsl_matrix_get(model->sqrt1_Hii, j, k);
                bt = MAX(bt, 0.0);
                if (model->mmRef->model == LOGIT) bt = MIN(bt, 1.0);
                gsl_matrix_set(bT, j, k, bt);
        }   }	    
	break;
   case PERMUTE: 
        if (bootID==NULL) 
            gsl_ran_shuffle(rnd,permid,nRows,sizeof(unsigned int));
        for (j=0; j<nRows; j++) {
            if (bootID==NULL) id = permid[j];
            else id = (unsigned int) gsl_matrix_get(bootID, i, j);
	    // bY = mu + bootr * sqrt(var)
	    for (k=0; k<nVars; k++) {
                bt=gsl_matrix_get(model->Mu,j,k)+sqrt(gsl_matrix_get(model->Var,j,k))*gsl_matrix_get(model->Res, id, k)*gsl_matrix_get(model->sqrt1_Hii, j, k);
            bt = MAX(bt, 0.0);
            if (model->mmRef->model == LOGIT) bt = MIN(bt, 1.0);
            gsl_matrix_set(bT, j, k, bt);
        }   }
        break;
   case FREEPERM:
         if (bootID==NULL) 
             gsl_ran_shuffle(rnd,permid,nRows,sizeof(unsigned int));
         for (j=0; j<nRows; j++) {
              if (bootID==NULL)  id = permid[j];
              else id = (unsigned int) gsl_matrix_get(bootID, i, j);
              yj=gsl_matrix_row(model->Yref, id);
              gsl_matrix_set_row (bT, j, &yj.vector);
 	 }
	 break;
   case MONTECARLO:
        if (model->mmRef->model == POISSON) {
           for (j=0; j<nRows; j++) 
           for (k=0; k<nVars; k++) {
                eij = gsl_matrix_get(XBeta, j, k);
                yij = gsl_ran_poisson(rnd, exp(eij));
                gsl_matrix_set(bT, j, k, yij);
        }   }
        else if (model->mmRef->model == LOGIT) {
//           displaymatrix(Sigma, "Sigma");
//           displaymatrix(XBeta, "Xbeta");
           for (j=0; j<nRows; j++) {
               yj = gsl_matrix_row(bT, j);
               // get random effect ej*         
            //   rmvnorm(rnd, nVars, Sigma, &yj.vector);
               semirmvnorm(rnd, nVars, Sigma, &yj.vector);
               for (k=0; k<nVars; k++) {
                   eij = gsl_matrix_get (XBeta, j, k); // logit(m)
                   eij = eij + gsl_vector_get(&yj.vector, k); 
                   yij = gsl_ran_bernoulli(rnd, exp(eij)/(1+exp(eij)));
                   gsl_matrix_set(bT, j, k, yij);
        }   }   }
        else if (model->mmRef->model == NB){  // Poisson log-normal
           for (j=0; j<nRows; j++) {
               yj = gsl_matrix_row(bT, j);
               // get random effect ej*         
               rmvnorm(rnd, nVars, Sigma, &yj.vector);
               for (k=0; k<nVars; k++) {
                   eij=gsl_matrix_get (XBeta, j, k);
                   // m_j = X_j * Beta_j + random_effect
                   if ( model->phi[k]>0 ) // add random effect
                      eij = eij + gsl_vector_get(&yj.vector, k);
                      // printf("%.2f\t", exp(mij));
                      yij = gsl_ran_poisson(rnd, exp(eij));
                      gsl_matrix_set(bT, j, k, yij);
         }   }   }   //  printf("\n");                     
        else 
            GSL_ERROR("The model type is not supported", GSL_ERANGE); 
        break;
    default: GSL_ERROR("The resampling method is not supported", GSL_ERANGE); break;
    }

    return SUCCESS;
} 

int GlmTest::setMonteCarlo(glm *model, gsl_matrix *Os, gsl_matrix *R)
{
   unsigned int j;
   double vij, sd, scale;
   double k = 16*sqrt(3)/15/M_PI;

   gsl_matrix *Sd = gsl_matrix_alloc (nVars, nVars);
   gsl_vector  *s = gsl_vector_alloc (nVars);
   gsl_vector  *var = gsl_vector_alloc (nVars);
   gsl_vector_view mj, d;

   if (model->mmRef->model == POISSON) {
       // Assuming no random effects, i.e. e*=0
       // M = X * Beta
       gsl_matrix_memcpy (XBeta, model->Eta);
       gsl_matrix_set_identity (Sigma);
    }          
   else if (model->mmRef->model == LOGIT) { 
       // Adjusting all betas and assuming var=1 
       // logit(M) = X * sqrt(1+0.346 var) beta = 1.1601 Eta 
       // See MATH5885 LDA lecture notes W9-11, Section 6.6.1
       sd = 1;
       scale = sqrt(1 + gsl_pow_2(k)*gsl_pow_2(sd));
//     printf("scale=%.4f\n", scale);
       gsl_matrix_memcpy (XBeta, model->Eta);
       gsl_matrix_scale (XBeta, scale);
       gsl_matrix_memcpy (Sigma, R);
   }
   else if (model->mmRef->model == NB) {
       // Adjusting the intercept to account for random effects
       //  i.e., M = X * beta - 0.5 * var 
       // var = log(1+phi)    
       for ( j=0; j<nVars; j++) {
           vij = log(1+model->phi[j]);
           gsl_vector_set(var, j, vij);
           gsl_vector_set(s, j, sqrt(vij));
          // adjust E(mj) = X*beta for the random effects
           mj=gsl_matrix_column (XBeta, j);
           gsl_vector_set_all(&mj.vector, -0.5*vij);         
           //  printf("%.2f ", model->phi[j]);
       }
       gsl_matrix_set_zero (Sd);
       gsl_blas_dger (1.0, s, s, Sd);    

       // if phi=0, then vij=0, i.e., no random effects (independence)
       // So it has zero impact / correlation on other variables 
       d = gsl_matrix_diagonal(Sd);
       for (j=0; j<nVars; j++) {
           if (model->phi[j]==0) 
               gsl_vector_set(&d.vector, j, 1.0);
       }
       // displaymatrix(Sd, "Sd");
       
       // Sigma = diag(var)*R*diag(var)
       gsl_matrix_memcpy(Sigma, R);
       gsl_matrix_mul_elements(Sigma, Sd);
       // displaymatrix(Sigma, "log-normal Sigma");
        
       gsl_matrix_add (XBeta, model->Eta);

       // free memory
        gsl_matrix_free(Sd);
        gsl_vector_free(s);
        gsl_vector_free(var);
   }
   else
       GSL_ERROR("The model type is not supported", GSL_ERANGE); 
   
    // subtract the offset
    if ( Os != NULL ) gsl_matrix_sub(XBeta, Os);

    return SUCCESS;
}

