// Interface between R and glmtest.cpp (ANOVA) (Rcpp API >= 0.8.6)
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 20-April-2011

#include <Rcpp.h>
extern "C"{
#include "resampTest.h"
#include "time.h"
}

RcppExport SEXP RtoGlmAnova(SEXP mpar, SEXP tpar, SEXP Ysexp, SEXP Xsexp, SEXP Osexp, 
                            SEXP INsexp, SEXP bIDsexp, SEXP LamSexp) 
{
    using namespace Rcpp;

    // Get parameters in params.
    List sparam(mpar);
    // pass parameters
    reg_Method mm;
    mm.tol = as<double>(sparam["tol"]);
    mm.model = as<unsigned int>(sparam["regression"]);
    mm.estiMethod = as<unsigned int>(sparam["estimation"]);
    mm.varStab = as<unsigned int>(sparam["stablizer"]);   
    mm.n = as<unsigned int>(sparam["n"]);   
    mm.maxiter = as<unsigned int>(sparam["maxiter"]);   
    mm.maxiter2 = as<unsigned int>(sparam["maxiter2"]);   
    mm.warning = as<unsigned int>(sparam["warning"]);   

    List rparam(tpar);
    // pass parameters
    mv_Method tm;	
    tm.nboot = as<unsigned int>(rparam["nboot"]);
    tm.corr = as<unsigned int>(rparam["cor_type"]);
    tm.test = as<unsigned int>(rparam["test_type"]);
    tm.resamp = as<unsigned int>(rparam["resamp"]);
    tm.reprand = as<unsigned int>(rparam["reprand"]);
    tm.punit = as<unsigned int>(rparam["punit"]);
    tm.showtime = as<unsigned int>(rparam["showtime"]);
    tm.warning = as<unsigned int>(rparam["warning"]);   

//  for debug
//    Rprintf("Input param arguments:\n tol=%g, mm.model=%d, nboot=%d, cor_type=%d, shrink_param=%g, test_type=%d, resamp=%d, n=%d, showtime=%d\n", mm.tol, mm.model, tm.nboot, tm.corr, tm.shrink_param, tm.test, tm.resamp, mm.n, tm.showtime);

    NumericMatrix Yr(Ysexp);
    NumericMatrix Xr(Xsexp);
    NumericMatrix Or(Osexp);
    IntegerMatrix INr(INsexp);
    NumericVector lambda(LamSexp);

    unsigned int nRows = Yr.nrow();
    unsigned int nVars = Yr.ncol();
    unsigned int nParam = Xr.ncol();
    unsigned int nModels = INr.nrow();
    unsigned int nLambda = lambda.size();
    tm.nRows = nRows;
    tm.nVars = nVars;
    tm.nParam = nParam;
// for debug
//    Rprintf("nRows=%d, nVars=%d, nParam=%d\n", nRows, nVars, nParam);

    // Rcpp -> gsl
    unsigned int i, j, k;
    tm.anova_lambda = gsl_vector_alloc(nLambda);
    for (i=0; i<nLambda; i++) {
        gsl_vector_set(tm.anova_lambda, i, lambda(i));
//    	Rprintf("%.2f ", gsl_vector_get(tm.anova_lambda, i));
    }
//    Rprintf("\n");

    gsl_matrix *X = gsl_matrix_alloc(nRows, nParam);        
    gsl_matrix *Y = gsl_matrix_alloc(nRows, nVars);
    gsl_matrix *O = gsl_matrix_alloc(nRows, nVars);
    gsl_matrix *isXvarIn = gsl_matrix_alloc(nModels, nParam);
//  Must be careful using std::copy for matrix. The following direct
//  use is not right as row elements are copied column-wise.  
//    std::copy( Yr.begin(), Yr.end(), Y->data );
//    std::copy( Xr.begin(), Xr.end(), X->data );
//    std::copy( INr.begin(), INr.end(), isXvarIn->data );

    for (i=0; i<nRows; i++) {
        for (j=0; j<nVars; j++) {
            gsl_matrix_set(Y, i, j, Yr(i, j));
            gsl_matrix_set(O, i, j, Or(i, j));
//            Rprintf("%d ", (int)gsl_matrix_get(Y, i, j));
        }
//        Rprintf("\t");
        
        for (k=0; k<nParam; k++){
            gsl_matrix_set(X, i, k, Xr(i, k));
//            Rprintf("%.2f ", gsl_matrix_get(X, i, k));
        }
//        Rprintf("\n");    
    }

    for (i=0; i<nModels; i++) {
        for (j=0; j<nParam; j++){
            gsl_matrix_set(isXvarIn, i, j, INr(i, j));
//            Rprintf("%d ", (int)gsl_matrix_get(isXvarIn, i, j));
        }
//        Rprintf("\n");
    }

    // do stuff	
  clock_t clk_start, clk_end;
  clk_start = clock();

    // glmfit
    PoissonGlm pfit(&mm);
    NBinGlm nbfit(&mm);
    BinGlm binfit(&mm);
    glm *glmPtr[3] = { &pfit, &nbfit, &binfit };
    unsigned int mtype = mm.model-1;
    glmPtr[mtype]->regression(Y, X, O, NULL);
//    glmPtr[mtype]->display();

    GlmTest myTest(&tm);
    // Resampling indices
    if ( !Rf_isNumeric(bIDsexp) || !Rf_isMatrix(bIDsexp) ) {
//       Rprintf("Calc bootID on the fly.\n");
     }
    else {
        if ( tm.resamp == SCOREBOOT ) {
            NumericMatrix bIDr(bIDsexp);
            tm.nboot = bIDr.nrow();	   
            myTest.bootID = gsl_matrix_alloc(tm.nboot, nRows);
//	    std::copy ( bIDr.begin(), bIDr.end(), anova.bootID->data);
	    for (i=0; i<tm.nboot; i++)
	    for (j=0; j<nRows; j++)
                gsl_matrix_set(myTest.bootID, i, j, bIDr(i, j));
	 }
	else{
	    IntegerMatrix bIDr(bIDsexp);
            tm.nboot = bIDr.nrow();	   
	    myTest.bootID = gsl_matrix_alloc(tm.nboot, nRows);
	    // integer -> double
	    for (i=0; i<tm.nboot; i++)
            for (j=0; j<nRows; j++)
                gsl_matrix_set(myTest.bootID, i, j, bIDr(i, j)-1);
    }  } 

    // resampling test
    myTest.anova(glmPtr[mtype],isXvarIn);
//    myTest.displayAnova();

    clk_end = clock();
    unsigned long int dif = floor((double)(clk_end - clk_start)/(double)(CLOCKS_PER_SEC));  
    unsigned int hours = floor((double)(dif/(double)3600));
    unsigned int min = floor((double)(dif%3600)/(double)60);
    unsigned int sec = dif%60;   
    if (tm.showtime>=TRUE)
        Rprintf("Time elapsed: %d hr %d min %d sec\n", hours, min, sec);

    // Wrap the gsl objects with Rcpp 
    NumericVector Vec_df(myTest.dfDiff, myTest.dfDiff+nModels-1);
    NumericVector Vec_mul(nModels-1);
    NumericVector Vec_Pmul(nModels-1);
   
    NumericMatrix Mat_statj(nModels-1, nVars);
    NumericMatrix Mat_Pstatj(nModels-1, nVars);
    for (i=0; i<nModels-1; i++) {
        Vec_mul(i) = gsl_matrix_get(myTest.anovaStat, i, 0);
        Vec_Pmul(i) = gsl_matrix_get(myTest.Panova, i, 0);
        for (j=0; j<nVars; j++){
            Mat_statj(i, j) = gsl_matrix_get(myTest.anovaStat, i, j+1);
	    Mat_Pstatj(i, j) = gsl_matrix_get(myTest.Panova, i, j+1);        
	}    
    }
//    double *uj = gsl_matrix_ptr(anova.statj, 0, 0);
//    double *pj = gsl_matrix_ptr(anova.Pstatj, 0, 0);
//    std::copy(uj, uj+nVars*(nModels-1), Mat_statj.begin());
//    std::copy(pj, pj+nVars*(nModels-1), Mat_Pstatj.begin());
    
    unsigned int nSamp = myTest.nSamp;

    // Rcpp -> R
    List rs = List::create(
         _["multstat" ] = Vec_mul,
         _["Pmultstat"] = Vec_Pmul,
	 _["dfDiff"   ] = Vec_df,
	 _["statj"    ] = Mat_statj,
	 _["Pstatj"   ] = Mat_Pstatj,
	 _["nSamp"    ] = nSamp
    );

    // clear objects
    myTest.releaseTest();
    glmPtr[mtype]->releaseGlm();
    gsl_matrix_free(Y);
    gsl_matrix_free(X);
    gsl_matrix_free(O);
    gsl_matrix_free(isXvarIn);
    gsl_vector_free(tm.anova_lambda);

    return rs;
}

