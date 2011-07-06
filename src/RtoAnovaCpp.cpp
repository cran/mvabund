// Interface between R and anova.cpp (Rcpp API >= 0.7.11)
//
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// Last modified: 20-April-2010

#include <Rcpp.h>
extern "C"{
#include "resampTest.h"
#include "time.h"
}

RcppExport SEXP RtoAnovaCpp(SEXP params, SEXP Ysexp, SEXP Xsexp,  
                             SEXP INsexp, SEXP bIDsexp ) 
{
    using namespace Rcpp;

    // Get parameters in params.
    List rparam(params);
    // pass parameters
    mv_Method mm;	
    mm.tol = as<double>(rparam["tol"]);
    mm.nboot = as<unsigned int>(rparam["nboot"]);
    mm.corr = as<unsigned int>(rparam["cor_type"]);
    mm.shrink_param = as<double>(rparam["shrink_param"]);
    mm.test = as<unsigned int>(rparam["test_type"]);
    mm.resamp = as<unsigned int>(rparam["resamp"]);
    mm.reprand = as<unsigned int>(rparam["reprand"]);
    mm.student = as<unsigned int>(rparam["studentize"]);
    mm.punit = as<unsigned int>(rparam["punit"]);
    mm.rsquare = as<unsigned int>(rparam["rsquare"]);

// for debug
//    Rprintf("Input param arguments:\n tol=%g, nboot=%d, cor_type=%d, shrink_param=%g, test_type=%d, resamp=%d, reprand=%d\n",mm.tol, mm.nboot, mm.corr, mm.shrink_param, mm.test, mm.resamp, mm.reprand);

    IntegerMatrix Yr(Ysexp);
    NumericMatrix Xr(Xsexp);
    IntegerMatrix INr(INsexp);
    unsigned int nRows = Yr.nrow();
    unsigned int nVars = Yr.ncol();
    unsigned int nParam = Xr.ncol();
    unsigned int nModels = INr.nrow();
// for debug
//    Rprintf("nRows=%d, nVars=%d, nParam=%d\n", nRows, nVars, nParam);

    // Rcpp -> gsl
    unsigned int i, j, k;
    gsl_matrix *X = gsl_matrix_alloc(nRows, nParam);        
    gsl_matrix *Y = gsl_matrix_alloc(nRows, nVars);
    gsl_matrix *isXvarIn = gsl_matrix_alloc(nModels, nParam);

//  Must be careful about using std::copy for matrix. The following direct
//  use is not doing right - row elements are copied to columns. Need to
//  fix it later on. 
//    std::copy( Yr.begin(), Yr.end(), Y->data );
//    std::copy( Xr.begin(), Xr.end(), X->data );
//    std::copy( INr.begin(), INr.end(), isXvarIn->data );

    for (i=0; i<nRows; i++){
        for (j=0; j<nVars; j++) {
            gsl_matrix_set(Y, i, j, Yr(i, j));
//            Rprintf("%d ", (int)gsl_matrix_get(Y, i, j));
        }
//        Rprintf("\t");
//
        for (k=0; k<nParam; k++){
	    gsl_matrix_set(X, i, k, Xr(i, k));
//	    Rprintf("%.2f ", gsl_matrix_get(X, i, k));
	}
//	Rprintf("\n");
    }
       
    for (i=0; i<nModels; i++){
        for (j=0; j<nParam; j++){
            gsl_matrix_set(isXvarIn, i, j, INr(i, j));
//            Rprintf("%d ", (int)gsl_matrix_get(isXvarIn, i, j));
        }
//	Rprintf("\n");
    }

    // do stuff	
    clock_t clk_start, clk_end;
    clk_start = clock();

// initialize anova class
    AnovaTest anova(&mm, Y, X, isXvarIn);
	
// Resampling indices
    if ( !Rf_isNumeric(bIDsexp) || !Rf_isMatrix(bIDsexp) ) {
//      Rprintf("Calc bootID on the fly.\n");
     }
    else {
        if ( mm.resamp == SCOREBOOT ) {
            NumericMatrix bIDr(bIDsexp);
            mm.nboot = bIDr.nrow();	   
            anova.bootID = gsl_matrix_alloc(mm.nboot, nRows);
//	    std::copy ( bIDr.begin(), bIDr.end(), anova.bootID->data);
	    for (i=0; i<mm.nboot; i++)
	    for (j=0; j<nRows; j++)
                gsl_matrix_set(anova.bootID, i, j, bIDr(i, j));
	 }
	else{
	    IntegerMatrix bIDr(bIDsexp);
            mm.nboot = bIDr.nrow();	   
	    anova.bootID = gsl_matrix_alloc(mm.nboot, nRows);
	    // integer -> double
	    for (i=0; i<mm.nboot; i++)
            for (j=0; j<nRows; j++)
                gsl_matrix_set(anova.bootID, i, j, bIDr(i, j)-1);
    }  } 

    // resampling test
    anova.resampTest();
//    anova.display();

    clk_end = clock();
    double dif = (double)(clk_end - clk_start)/(double)(CLOCKS_PER_SEC);
    Rprintf("Time elapsed: %d seconds\n", (unsigned int) dif);

    // Wrap the gsl objects with Rcpp 
    NumericVector Vec_mul(anova.multstat, anova.multstat+nModels-1);
    NumericVector Vec_Pm(anova.Pmultstat, anova.Pmultstat+nModels-1);
    NumericVector Vec_df(anova.dfDiff, anova.dfDiff+nModels-1);
    
    NumericMatrix Mat_statj(nModels-1, nVars);
    NumericMatrix Mat_Pstatj(nModels-1, nVars);
    for (i=0; i<nModels-1; i++)
    for (j=0; j<nVars; j++){
        Mat_statj(i, j) = gsl_matrix_get(anova.statj, i, j);
	Mat_Pstatj(i, j) = gsl_matrix_get(anova.Pstatj, i, j);        
    }
//    double *uj = gsl_matrix_ptr(anova.statj, 0, 0);
//    double *pj = gsl_matrix_ptr(anova.Pstatj, 0, 0);
//    std::copy(uj, uj+nVars*(nModels-1), Mat_statj.begin());
//    std::copy(pj, pj+nVars*(nModels-1), Mat_Pstatj.begin());

    unsigned int nSamp = anova.nSamp;

    // Rcpp -> R
    List rs = List::create(
         _["multstat" ] = Vec_mul,
         _["Pmultstat"] = Vec_Pm,
	 _["dfDiff"   ] = Vec_df,
	 _["statj"    ] = Mat_statj,
	 _["Pstatj"   ] = Mat_Pstatj,
	 _["nSamp"    ] = nSamp
    );

    // clear objects
    anova.releaseTest(); 
    gsl_matrix_free(Y);
    gsl_matrix_free(X);
    gsl_matrix_free(isXvarIn);


    return rs;
}

