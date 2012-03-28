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
    double oij, lTol=fit->lTol;

    unsigned int nP, i, j, k, l;
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
    gsl_matrix *bXAlt=NULL, *bXnull=NULL, *bO=NULL;
//    gsl_matrix *XO=NULL, *BetaO=NULL;
//    gsl_matrix *Onull=NULL, *bOnull=NULL;
    gsl_matrix *bY=gsl_matrix_alloc(nRows, nVars);
    gsl_matrix_view X0;

    unsigned int mtype = fit->mmRef->model-1;
    // ======= Fit the (first) Alt model =========//
    for (i=0; i<nModels; i++) {
        nP = 0;
        for (k=0; k<nParam; k++) 
	     if (gsl_matrix_get(Xin,i,k)!=FALSE) nP++;   
        rdf[i] = nRows-nP;
    }
    PtrAlt[mtype]->copyGlm(fit);

     for (i=1; i<nModels; i++) {       
        // ======= Fit the Null model =========//
        ID0 = i; ID1 = i-1;
        nP0 = nRows - (unsigned int)rdf[ID0];
        nP1 = nRows - (unsigned int)rdf[ID1];

        ref1=gsl_matrix_row(Xin, ID1);
        ref0=gsl_matrix_row(Xin, ID0);
        X0=gsl_matrix_submatrix(fit->Xref,0,0,nRows,nP0);
        PtrNull[mtype]->EstIRLS(fit->Yref, &X0.matrix, NULL, PtrAlt[mtype]->phi); 
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
        else 
           subGeeLR(PtrAlt[mtype], PtrNull[mtype], &teststat.vector);	

	// ======= Get univariate test statistics =======//
        if (tm->punit == FREESTEP) {	
              unitstat=gsl_vector_subvector(&teststat.vector,1,nVars);
              gsl_sort_vector_index (sortid, &unitstat.vector);
              gsl_permutation_reverse(sortid);        
        }
        // ======= Prepare for resampling methods ======//
	if ( tm->resamp == CASEBOOT ) {
            bXAlt = gsl_matrix_alloc(nRows, nP1); 
            bXnull = gsl_matrix_alloc(nRows, nP0); 
            bO = gsl_matrix_alloc(nRows, nVars);
//            XO = gsl_matrix_alloc(nRows, nP1-nP0); 
//            BetaO = gsl_matrix_alloc(nP1-nP0, nVars); 
//            BetaO = gsl_matrix_alloc(nP0, nVars); 
//            Onull = gsl_matrix_alloc(nRows, nVars);
//            bOnull = gsl_matrix_alloc(nRows, nVars);
            PtrAlt[mtype]->Oref = gsl_matrix_alloc(nRows, nVars);

            gsl_matrix_memcpy(PtrAlt[mtype]->Oref, PtrAlt[mtype]->Eta);
            ref0=gsl_matrix_subrow(Xin, ID0, 0, nP1);            
//            subX1(PtrAlt[mtype]->Xref, &ref0.vector, XO);
//            subXrow2(PtrAlt[mtype]->Beta, &ref0.vector, BetaO);
//            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,XO,BetaO,0.0,Onull);
//            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,PtrNull[mtype]->Xref,BetaO,0.0,Onull);
//            for (k=0; k<nRows; k++)
//            for (l=0; l<nVars; l++){
//                oij=gsl_matrix_get(Onull, k, l);
//                oij=(oij<-lTol)?-lTol:((oij>lTol)?lTol:oij);
//                gsl_matrix_set(Onull, k, l, oij);
//            }
	}
       
	if (tm->resamp == MONTECARLO)  {
           if (tm->corr != SHRINK)  
              GetR(PtrAlt[mtype]->Res, SHRINK, tm->anova_lambda, ID1, Rlambda);
           setMonteCarlo (PtrNull[mtype], NULL, Rlambda);
        }

        // ======= Get resampling distribution under H0 ===== //
	nSamp=0;
        for (j=0; j<tm->nboot; j++) {	
	    gsl_vector_set_zero (bStat);
	    if ( tm->resamp == CASEBOOT ) {
//                resampAnovaCase(PtrAlt[mtype], Onull, bY, bXAlt, bO, bOnull, j);
                resampAnovaCase(PtrAlt[mtype], bY, bXAlt, bO, j);
	        bAlt[mtype]->regression(bY, bXAlt, bO);
                if (tm->test != WALD){
                   subX(bXAlt, &ref0.vector, bXnull);
//                   gsl_matrix_sub (bO, bOnull);
                   bNull[mtype]->EstIRLS(bY,bXnull,bO,bAlt[mtype]->phi);
//                   bNull[mtype]->EstIRLS(bY,bXnull,bOnull,bAlt[mtype]->phi);
                }    
           }		
           else {
		resampNonCase(PtrNull[mtype], bY, j);
                bAlt[mtype]->regression(bY,PtrAlt[mtype]->Xref,NULL); 
                if (tm->test != WALD) 
                    bNull[mtype]->EstIRLS(bY,PtrNull[mtype]->Xref,NULL,bAlt[mtype]->phi);                
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
         if (bXAlt!=NULL) gsl_matrix_free(bXAlt);   
         if (bXnull!=NULL) gsl_matrix_free(bXnull);   
         if (bO!=NULL) gsl_matrix_free(bO);   
//         if (XO!=NULL) gsl_matrix_free(XO);   
//         if (BetaO!=NULL) gsl_matrix_free(BetaO);
//         if (Onull!=NULL) gsl_matrix_free(Onull);
//         if (bOnull!=NULL) gsl_matrix_free(bOnull);
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
        if ( tm->resamp == CASEBOOT )
           PtrNull->EstIRLS(PtrAlt->Yref,GrpXs[k].matrix,GrpOs[0].matrix,PtrAlt->phi); // always under H1
        else 
           PtrNull->EstIRLS(PtrAlt->Yref,GrpXs[k].matrix,GrpOs[k].matrix,PtrAlt->phi); // resample under H0
        GetR(PtrNull->Res, tm->corr, tm->smry_lambda, k, R);  
        gsl_vector_view teststat=gsl_matrix_row(Stats, k-1);
        subGeeScore(PtrAlt, PtrNull, &teststat.vector, R);
    }     	        
   return SUCCESS;
}

int GlmTest::geeLR(glm *PtrAlt, glm *PtrNull, gsl_matrix *Stats)
{
    for (unsigned int k=1; k<nParam+2; k++) {
        if ( tm->resamp == CASEBOOT )
           PtrNull->EstIRLS(PtrAlt->Yref,GrpXs[k].matrix,GrpOs[0].matrix,PtrAlt->phi); // always under H1
        else 
           PtrNull->EstIRLS(PtrAlt->Yref,GrpXs[k].matrix,GrpOs[k].matrix,PtrAlt->phi); // resample under H0
        gsl_vector_view teststat=gsl_matrix_row(Stats, k-1);
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
    double result, alpha, sum=0;
    unsigned int i, j, l, nP = PtrAlt->Xref->size2;

    gsl_vector *U = gsl_vector_alloc(nVars*nP);
    gsl_matrix *kRlNull = gsl_matrix_alloc(nVars*nP, nVars*nP);
    gsl_matrix_set_zero (kRlNull);
    gsl_matrix *XwX = gsl_matrix_alloc(nP, nP);
    gsl_vector *tmp=gsl_vector_alloc(nVars*nP);
    gsl_vector_view wj, uj, rj, tmp2;
    gsl_matrix_view Rl;

    GrpMat *Z = (GrpMat*)malloc(nVars*sizeof(GrpMat));
    for (j=0; j<nVars; j++) {
        Z[j].matrix = gsl_matrix_alloc(nRows, nP);
        // get W^1/2 * X
        wj = gsl_matrix_column (PtrNull->wHalf, j);
        for (i=0; i<nP; i++)
            gsl_matrix_set_col (Z[j].matrix, i, &wj.vector);
        gsl_matrix_mul_elements (Z[j].matrix, PtrAlt->Xref);

        uj=gsl_vector_subvector(U, j*nP, nP);
        rj=gsl_matrix_column(PtrNull->Res, j);
        gsl_blas_dgemv(CblasTrans, 1, Z[j].matrix, &rj.vector, 0, &uj.vector);

        if ( (tm->punit>0) || (tm->corr==IDENTITY) ) {
           gsl_matrix_set_zero(XwX);
           gsl_blas_dsyrk(CblasLower, CblasTrans, 1, Z[j].matrix, 0, XwX);
           // univariate test = U^T*(XwX)^-1*U
           tmp2=gsl_vector_subvector(tmp, 0, nP);
           gsl_linalg_cholesky_decomp(XwX);
           gsl_linalg_cholesky_solve(XwX, &uj.vector, &tmp2.vector);
           gsl_blas_ddot(&uj.vector, &tmp2.vector, &result);
           gsl_vector_set(teststat, j+1, result);
           sum = sum+result;           
        }

        if ( tm->corr!=IDENTITY) {
            for (l=0; l<=j; l++) { // lower half
                alpha = gsl_matrix_get(R, j, l);
                Rl=gsl_matrix_submatrix(kRlNull,j*nP,l*nP,nP,nP);
                gsl_blas_dgemm(CblasTrans, CblasNoTrans, alpha, Z[j].matrix, Z[l].matrix, 0, &Rl.matrix);
            }
        }
    } // end for j=1:nVars

    // multivariate test stat   
    if ( tm->corr==IDENTITY )  
        gsl_vector_set(teststat, 0, sum);
    else {    
        gsl_linalg_cholesky_decomp (kRlNull);
        gsl_linalg_cholesky_solve (kRlNull, U, tmp);
        gsl_blas_ddot (U, tmp, &result);
        gsl_vector_set(teststat, 0, result);
    }

   // clear memory
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
    double alpha, result, sum=0;
    unsigned int nP = Alt->nParams;
    unsigned int nDF = LL->size1;

    gsl_vector *LBeta = gsl_vector_alloc(nVars*nDF);
    gsl_vector_set_zero(LBeta);
    gsl_matrix *w1jX1=gsl_matrix_alloc(nRows, nP);
    gsl_matrix *XwX=gsl_matrix_alloc(nP, nP);
    gsl_matrix *Rl2 = gsl_matrix_alloc(nDF, nP);
    gsl_matrix *IinvN = gsl_matrix_alloc(nDF, nDF);
    gsl_matrix *IinvRl = gsl_matrix_alloc(nVars*nDF, nVars*nDF);
    gsl_vector *tmp = gsl_vector_alloc(nVars*nDF);
    gsl_vector_view tmp2, wj, LBj, bj; 
    gsl_matrix_view Rl;

    gsl_matrix_set_zero(IinvRl);
    GrpMat *Z = (GrpMat*)malloc(nVars*sizeof(GrpMat));
    for (j=0; j<nVars; j++){
       Z[j].matrix = gsl_matrix_alloc(nP, nRows);
       // Get X * W^1/2
       wj=gsl_matrix_column(Alt->wHalf, j);
       for (i=0; i<nP; i++)
           gsl_matrix_set_col (w1jX1, i, &wj.vector);
       gsl_matrix_mul_elements (w1jX1, Alt->Xref);

       // Z = (X^T W X)^-1 * X^T W^1/2. 
       gsl_matrix_transpose_memcpy (Z[j].matrix, w1jX1); // X^T W^1/2
       gsl_matrix_set_zero (XwX);
       gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, w1jX1, 0.0, XwX);
       gsl_linalg_cholesky_decomp (XwX);
       gsl_blas_dtrsm (CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, XwX, Z[j].matrix);

       // LBeta = L*Beta       
       LBj=gsl_vector_subvector (LBeta, j*nDF, nDF);
       bj=gsl_matrix_column (Alt->Beta, j);
       gsl_blas_dgemv(CblasNoTrans,1,LL,&bj.vector,0,&LBj.vector);

       if ( (tm->punit>0) || (tm->corr==IDENTITY) ){
          // statj=LBeta'*(IinvN)^-1*LBeta
          gsl_linalg_cholesky_invert (XwX); // inv(X^T*W*X)
          gsl_linalg_cholesky_decomp (XwX); // A*A^T = inv(X^T*W*X)
          gsl_matrix_memcpy(Rl2, LL);
          gsl_blas_dtrmm (CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, XwX, Rl2); // L*A
          // IinvN = L*inv(X^T*W*X)*L^T = (L*A)*(L*A)^T
          gsl_matrix_set_zero(IinvN);
          gsl_blas_dsyrk (CblasLower, CblasNoTrans, 1.0, Rl2, 0.0, IinvN);
          gsl_linalg_cholesky_decomp (IinvN);
          tmp2=gsl_vector_subvector(tmp, 0, nDF);
          gsl_vector_memcpy(&tmp2.vector, &LBj.vector);
          gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasNonUnit, IinvN, &tmp2.vector); // inv(IinvN)*LBeta
          gsl_blas_ddot (&LBj.vector, &tmp2.vector, &result);
          gsl_vector_set(teststat, j+1, sqrt(result));
          sum = sum + result;
       }
       if (tm->corr!=IDENTITY) {
          // IinvRl=L*vSandRl*L^T 
          for (l=0; l<=j; l++) {
              Rl=gsl_matrix_submatrix(IinvRl,j*nDF,l*nDF,nDF,nDF);
              alpha = gsl_matrix_get(R, j, l);
              gsl_blas_dgemm(CblasNoTrans,CblasTrans,alpha,Z[j].matrix,Z[l].matrix, 0.0, XwX);  // borrow XwX space to store vSandRl 
              // Rl2 = L*vSandRl
              gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,LL,XwX,0.0,Rl2);
              // Rl = L*vSandRl*L^T = (L*A)*(L*A)^T
              gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,Rl2,LL,0.0,&Rl.matrix);
          } // end l
       }  // end if (tm->corr) 
    } // end for j=1:nVars       

    // Get multivariate test stat = LBeta^T * inv(IinvRl) * LBeta
    if ( tm->corr==IDENTITY ) 
        gsl_vector_set(teststat, 0, sqrt(sum));
    else {
        gsl_linalg_cholesky_decomp (IinvRl);
        gsl_linalg_cholesky_solve (IinvRl, LBeta, tmp);
        gsl_blas_ddot (LBeta, tmp, &result);
        gsl_vector_set(teststat, 0, sqrt(result));
    }

    // free memory
    for (j=0; j<nVars; j++) 
        gsl_matrix_free(Z[j].matrix);
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
    double det;
    unsigned int j, k, id, nP, isSingular;
    gsl_vector_view yj, oj, xj;
    GrpMat *GrptX = NULL;
    
    if ( tm->resamp == CASEBOOT ){
        if (bootID == NULL) {
            GrptX = (GrpMat *)malloc((nParam+2)*sizeof(GrpMat));
            for (k=0; k<nParam+2; k++){
               nP = GrpXs[k].matrix->size2;
               GrptX[k].matrix = gsl_matrix_alloc(nP, nP);
            }
            isSingular=TRUE;
            while (isSingular==TRUE) {
                for (j=0; j<nRows; j++) {
                    id = gsl_rng_uniform_int(rnd, nRows);
	            // resample Y, X, offsets accordingly
                    yj = gsl_matrix_row(model->Yref, id);
                    gsl_matrix_set_row(bT, j, &yj.vector);
               	    for (k=0; k<nParam+2; k++) {
                        xj = gsl_matrix_row(oriXs[k].matrix, id);
                        oj = gsl_matrix_row(oriOs[k].matrix, id);
                        gsl_matrix_set_row(GrpXs[k].matrix, j, &xj.vector); 
 	                gsl_matrix_set_row(GrpOs[k].matrix, j, &oj.vector);
                    }
                 }
                 for (k=0; k<nParam+2; k++) {
                     gsl_matrix_set_zero(GrptX[k].matrix);             
                     gsl_blas_dsyrk (CblasLower,CblasTrans,1.0,GrpXs[k].matrix,0.0,GrptX[k].matrix);
                     det = calcDet(GrptX[k].matrix); 
                     if (det>1e-8) {
                        isSingular=FALSE;
                        break;
            }     }    }
        }
        else {
             for (j=0; j<nRows; j++) {   
                 id = (unsigned int) gsl_matrix_get(bootID, i, j);
                 // resample Y and X and offset
                 yj=gsl_matrix_row(model->Yref, id);
                 gsl_matrix_set_row (bT, j, &yj.vector);
          	 for (k=0; k<nParam+2; k++) {
                     xj = gsl_matrix_row(oriXs[k].matrix, id);
                     oj = gsl_matrix_row(oriOs[k].matrix, id);
                     gsl_matrix_set_row(GrpXs[k].matrix, j, &xj.vector); 
 	             gsl_matrix_set_row(GrpOs[k].matrix, j, &oj.vector);
       }    }    }
    }
    else resampNonCase(model, bT, i);  

    if (GrptX!=NULL) {
       for (k=0; k<nParam+2; k++) 
           gsl_matrix_free(GrptX[k].matrix);
       free(GrptX);
    }    

    return SUCCESS;
}

//int GlmTest::resampAnovaCase(glm *model, gsl_matrix *Onull, gsl_matrix *bT, gsl_matrix *bX, gsl_matrix *bO, gsl_matrix *bOnull, unsigned int i)
int GlmTest::resampAnovaCase(glm *model, gsl_matrix *bT, gsl_matrix *bX, gsl_matrix *bO, unsigned int i)
{
    double det;
    unsigned int j, id, isSingular, nP;
    gsl_vector_view yj, xj, oj; //, o0j;
    nP = model->Xref->size2;
    gsl_matrix *txX = gsl_matrix_alloc(nP, nP);
    gsl_matrix_set_zero(txX);

    if (bootID == NULL) {
       isSingular=TRUE;
       while (isSingular==TRUE) {
            for (j=0; j<nRows; j++) {   
                id = gsl_rng_uniform_int(rnd, nRows);
                // resample Y and X and offset
                yj=gsl_matrix_row(model->Yref, id);
                xj = gsl_matrix_row(model->Xref, id);
                oj = gsl_matrix_row(model->Oref, id);
//                o0j = gsl_matrix_row(Onull, id);
                gsl_matrix_set_row (bT, j, &yj.vector);
                gsl_matrix_set_row(bX, j, &xj.vector);
                gsl_matrix_set_row(bO, j, &oj.vector);
//                gsl_matrix_set_row(bOnull, j, &o0j.vector);
             }
             gsl_blas_dsyrk (CblasLower, CblasTrans, 1.0, bX, 0.0, txX);
             det = calcDet(txX); 
             if (det>1e-8) isSingular=FALSE;
       } 
   }   		    	
   else {
       for (j=0; j<nRows; j++) {   
          id = (unsigned int) gsl_matrix_get(bootID, i, j);
          // resample Y and X and offset
          yj=gsl_matrix_row(model->Yref, id);
          xj = gsl_matrix_row(model->Xref, id);
          oj = gsl_matrix_row(model->Oref, id);
//          o0j = gsl_matrix_row(Onull, id);
          gsl_matrix_set_row (bT, j, &yj.vector);
          gsl_matrix_set_row(bX, j, &xj.vector);
          gsl_matrix_set_row(bO, j, &oj.vector);
//          gsl_matrix_set_row(bOnull, j, &o0j.vector);
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
                      yij = gsl_ran_poisson(rnd, exp(eij));
                      gsl_matrix_set(bT, j, k, yij);
         }   }   }                        
        else GSL_ERROR("The model type is not supported", GSL_ERANGE); 
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
   else GSL_ERROR("The model type is not supported", GSL_ERANGE); 
   
   // subtract the offset
   if ( Os != NULL ) gsl_matrix_sub(XBeta, Os);

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
