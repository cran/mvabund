#include "mvRcpp.h"
extern "C"{
#include "resampTest.h"
#include "time.h"
}

// params for methods, Matrices X, Y, Res, Coef are from the lm model
// nlist is numerical list including model.rdf, testStats, etc.
// df_Hats is a dataframe with concatenated matrices 
// Mat_isXvarIn is used for anova
RcppExport SEXP RtoAnovaCpp(SEXP params, SEXP Ysexp, SEXP Xsexp,  
                             SEXP INsexp, SEXP bIDsexp )//, 
//			     SEXP nl_stats)
//			     SEXP df_Hats, SEXP Mat_isXvarIn)
{
    SEXP rl=R_NilValue; // Use this when there is nothing to be returned.
    char *exceptionMesg=NULL;

    try {
	// Get parameters in params.
	RcppParams rparam(params);

        // pass parameters
        mv_Method mm;	
	mm.tol = rparam.getDoubleValue("tol");
	mm.nboot = rparam.getIntValue("nboot");
	mm.corr = rparam.getIntValue("cor_type");
	mm.shrink_param = rparam.getDoubleValue("shrink_param");
        mm.test = rparam.getIntValue("test_type");
	mm.resamp = rparam.getIntValue("resamp");
	mm.reprand = rparam.getIntValue("reprand");
	mm.student = rparam.getIntValue("studentize");
	mm.punit = rparam.getIntValue("punit");
	mm.rsquare = rparam.getIntValue("rsquare");

	// for debug
//	Rprintf("Input param arguments:\n tol=%g, nboot=%d, cor_type=%d, shrink_param=%g, test_type=%d, resamp=%d, reprand=%d\n",mm.tol, mm.nboot, mm.corr, mm.shrink_param, mm.test, mm.resamp, mm.reprand);


        RcppMatrixView<int> Yr(Ysexp);
	RcppMatrixView<double> Xr(Xsexp);
	RcppMatrixView<int> INr(INsexp);
//	Rprintf("nModels=%d\n", nModels);
        int nRows = Yr.dim1();
	int nVars = Yr.dim2();
	int nParam = Xr.dim2();
	int nModels = INr.dim1();
	// for debug
//	Rprintf("nRows=%d, nVars=%d, nParam=%d\n", nRows, nVars, nParam);

        // Rcpp -> gsl
	int i, j, k;
	gsl_matrix *X = gsl_matrix_alloc(nRows, nParam);        
	gsl_matrix *Y = gsl_matrix_alloc(nRows, nVars);
	gsl_matrix *isXvarIn = gsl_matrix_alloc(nModels, nParam);
        for (i=0; i<nRows; i++){
	    for (j=0; j<nVars; j++) {
	        gsl_matrix_set(Y, i, j, (int)Yr(i, j));
//		Rprintf("%d ", (int)gsl_matrix_get(Y, i, j));
            }
//	    Rprintf("\t");
            for (k=0; k<nParam; k++){
	        gsl_matrix_set(X, i, k, (double)Xr(i, k));
//	        Rprintf("%.2f ", gsl_matrix_get(X, i, k));
	    }
//	    Rprintf("\n");
         }
       
	for (i=0; i<nModels; i++){
	for (j=0; j<nParam; j++){
            gsl_matrix_set(isXvarIn, i, j, (int)INr(i, j));
//	    Rprintf("%d ", (int)gsl_matrix_get(isXvarIn, i, j));
	}
//	Rprintf("\n");
	}

        // do stuff	
	clock_t clk_start, clk_end;
	clk_start = clock();

        // initialize anova class
//	Rprintf("Initialize anova test\n");
        AnovaTest anova(&mm, Y, X, isXvarIn);
	
	// Resampling indices
        if ( !Rf_isNumeric(bIDsexp) || !Rf_isMatrix(bIDsexp) ) {
//	   Rprintf("Calc bootID on the fly.\n");
	 }
	else {
	   if ( mm.resamp == SCOREBOOT ) {
	      RcppMatrixView<double> bIDr(bIDsexp);
              mm.nboot = bIDr.dim1();	   
	      anova.bootID = gsl_matrix_alloc(mm.nboot, nRows);
	      for (i=0; i<mm.nboot; i++)
	      for (j=0; j<nRows; j++)
                  gsl_matrix_set(anova.bootID, i, j, bIDr(i, j));
	   }
	   else{
	      RcppMatrixView<int> bIDr(bIDsexp);
              mm.nboot = bIDr.dim1();	   
	//      Rprintf("nboot = %d, nrow=%d\n", mm.nboot, bIDr.dim2());
	      anova.bootID = gsl_matrix_alloc(mm.nboot, nRows);
	      for (i=0; i<mm.nboot; i++)
              for (j=0; j<nRows; j++)
                  gsl_matrix_set(anova.bootID, i, j, bIDr(i, j)-1);
        }  } 

       // resampling test
  //     Rprintf("Start resampling test\n");
	anova.resampTest();
//	anova.display();

	clk_end = clock();
	float dif = (float)(clk_end - clk_start)/(float)(CLOCKS_PER_SEC);
//	Rprintf("Time elapsed: %.4f seconds\n", dif);

        // gsl -> Rcpp
	RcppVector<double> Vec_mul(nModels-1);
	RcppVector<double> Vec_Pm(nModels-1);
	RcppVector<double> Vec_df(nModels-1);

	RcppMatrix<double> Mat_statj(nModels-1, nVars);
	RcppMatrix<double> Mat_Pstatj(nModels-1, nVars);

        for (i=0; i<nModels-1; i++){
	    Vec_mul(i) = anova.multstat[i];
	    Vec_Pm(i) = anova.Pmultstat[i];
	    Vec_df(i) = anova.dfDiff[i];

            for (j=0; j<nVars; j++){
                Mat_statj(i, j) = gsl_matrix_get(anova.statj, i, j);
		Mat_Pstatj(i, j) = gsl_matrix_get(anova.Pstatj, i, j);
        }   }

        int nbootsrun = anova.nSamp;

        // clear objects
        anova.releaseTest(); 
        gsl_matrix_free(Y);
	gsl_matrix_free(X);
	gsl_matrix_free(isXvarIn);

	// Rcpp -> R
        RcppResultSet rs;
	rs.add("multstat", Vec_mul);
        rs.add("Pmultstat", Vec_Pm);
	rs.add("dfDiff", Vec_df);
	rs.add("statj", Mat_statj);
	rs.add("Pstatj", Mat_Pstatj);
	rs.add("nSamp", nbootsrun);
	rl = rs.getReturnList();
	
    } catch(std::exception& ex) {
	exceptionMesg = copyMessageToR(ex.what());
    } catch(...) {
	exceptionMesg = copyMessageToR("unknown reason");
    }
    
    if(exceptionMesg != NULL)
	Rf_error(exceptionMesg);

    return rl;
}

