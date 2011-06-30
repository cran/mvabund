// Interface between R and summary.cpp (Rcpp API >= 0.7.11 )
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// Last modified: 20-April-2010

#include "Rcpp.h"
extern "C"{
#include "resampTest.h"
#include "time.h"
}

RcppExport SEXP RtoSmryCpp(SEXP params, SEXP Ysexp, SEXP Xsexp,  
                             SEXP bIDsexp )
{
    using namespace Rcpp;

    // Get parameters in params.
    List rparam(params);
    // pass parameters
    mv_Method mm;	
    mm.tol = rparam["tol"];
    mm.nboot = rparam["nboot"];
    mm.corr = rparam["cor_type"];
    mm.shrink_param = rparam["shrink_param"];
    mm.test = rparam["test_type"];
    mm.resamp = rparam["resamp"];
    mm.reprand = rparam["reprand"];
    mm.student = rparam["studentize"];
    mm.punit = rparam["punit"];
    mm.rsquare = rparam["rsquare"];

//    // for debug
//   Rprintf("Input param arguments:\n tol=%g, nboot=%d, cor_type=%d, shrink_param=%g, test_type=%d, resamp=%d, reprand=%d\n",mm.tol, mm.nboot, mm.corr, mm.shrink_param, mm.test, mm.resamp, mm.reprand);

    IntegerMatrix Yr(Ysexp);
    NumericMatrix Xr(Xsexp);
    int nRows = Yr.nrow();
    int nVars = Yr.ncol();
    int nParam = Xr.ncol();

    // Rcpp -> gsl
    int i, j, k;
    gsl_matrix *X = gsl_matrix_alloc(nRows, nParam);        
    gsl_matrix *Y = gsl_matrix_alloc(nRows, nVars);    
    
//  Must be careful about using std::copy for matrix. The following direct
//  use is not doing right - row elements are copied to columns. Need to 
//  fix it later on.
//    std::copy( Yr.begin(), Yr.end(), Y->data );
//    std::copy( Xr.begin(), Xr.end(), X->data );

    for (i=0; i<nRows; i++)
    for (j=0; j<nVars; j++){ 
        gsl_matrix_set(Y, i, j, Yr(i, j));
        for (k=0; k<nParam; k++)
            gsl_matrix_set(X, i, k, Xr(i, k));
    }
       
    // do stuff	
//  clock_t clk_start, clk_end;
//  clk_start = clock();

    // initialize summary class
    Summary smry(&mm, Y, X);
	
    // Resampling indices
    if ( !Rf_isNumeric(bIDsexp) || !Rf_isMatrix(bIDsexp) ) {
//      Rprintf("Calc bootID on the fly.\n");
    }	   
    else {
        if ( mm.resamp == SCOREBOOT ) {
            NumericMatrix bIDr(bIDsexp);
            mm.nboot = bIDr.nrow();	   
            smry.bootID = gsl_matrix_alloc(mm.nboot, nRows);
//	    std::copy( bIDr.begin(), bIDr.end(), smry.bootID->data );
            for (i=0; i<mm.nboot; i++)
            for (j=0; j<nRows; j++)
                gsl_matrix_set(smry.bootID, i, j, bIDr(i, j));
	}
        else{
	    IntegerMatrix bIDr(bIDsexp);
            mm.nboot = bIDr.nrow();	   
	    smry.bootID = gsl_matrix_alloc(mm.nboot, nRows);
	    // integer -> double
	    for (i=0; i<mm.nboot; i++)
            for (j=0; j<nRows; j++)
                gsl_matrix_set(smry.bootID, i, j, bIDr(i, j)-1);
    }   } 

// resampling test
    smry.resampTest();
//    smry.display();

//  clk_end = clock();
//  float dif = (float)(clk_end - clk_start)/(float)(CLOCKS_PER_SEC);
//  Rprintf("Time elapsed: %.4f seconds\n", dif);

    // Wrap gsl vectors with Rcpp 
    NumericVector Vec_signific(smry.multstat+1, smry.multstat+nParam+1);
    NumericVector Vec_Psignific(smry.Pmultstat+1, smry.Pmultstat+nParam+1);

    // Copy gsl matrix to Rcpp Matrix
    NumericVector Vec_unitmult(nVars);
    NumericVector Vec_Punitmult(nVars);
    NumericMatrix Mat_unitsign(nParam, nVars);
    NumericMatrix Mat_Punitsign(nParam, nVars);

    double *uj = gsl_matrix_ptr(smry.unitstat, 0, 0);
    double *pj = gsl_matrix_ptr(smry.Punitstat, 0, 0);
    Vec_unitmult = wrap(uj, uj+nVars);
    Vec_Punitmult = wrap(pj, pj+nVars);

    for (i=0; i<nParam; i++){
        uj=gsl_matrix_ptr (smry.unitstat, i+1, 0);
	pj=gsl_matrix_ptr (smry.Punitstat, i+1, 0);	
        for (j=0; j<nVars; j++){
            Mat_unitsign(i, j) = *(uj+j);
	    Mat_Punitsign(i, j) = *(pj+j);
	}    
    }
//    std::copy(ptr_statj+nVars, ptr_statj+nVars*nParam-1, 
//              Mat_unitsign.begin());
//    std::copy(ptr_Pstatj+nVars+1, ptr_Pstatj+1+nVars*nParam, 
//              Mat_Punitsign.begin());
//    Rprintf("Done\n.");

    // clear objects
    smry.releaseSummary(); 
    gsl_matrix_free(Y);
    gsl_matrix_free(X);
    
    // Rcpp -> R
    List rs = List::create(
	_["multstat"] = smry.multstat[0],
	_["Pmultstat"] = smry.Pmultstat[0],
        _["unitmult"] = Vec_unitmult,
        _["Punitmult"]= Vec_Punitmult,
	_["signific"]= Vec_signific,
	_["Psignific"]= Vec_Psignific,
	_["unitsign"]= Mat_unitsign,
	_["Punitsign"]= Mat_Punitsign,
	_["nSamp"]= smry.nSamp,
	_["R2"]= smry.R2
    );
    return rs;
}

