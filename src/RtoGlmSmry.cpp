// Interface between R and glmtest.cpp (summary function) (Rcpp API >= 0.8.6 )
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 20-April-2011

#include <Rcpp.h>
extern "C"{
#include "resampTest.h"
#include "time.h"
//#include "math.h"
}

RcppExport SEXP RtoGlmSmry(SEXP mpar, SEXP tpar, SEXP Ysexp, SEXP Xsexp,  
                             SEXP bIDsexp, SEXP LamSexp )
{
    using namespace Rcpp;

    // Get parameters in params.
    List sparam(mpar);
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
    tm.corr = as<unsigned int>(rparam["cor_type"]);
    tm.test = as<unsigned int>(rparam["test_type"]);
    tm.resamp = as<unsigned int>(rparam["resamp"]);
    tm.reprand = as<unsigned int>(rparam["reprand"]);
    tm.punit = as<unsigned int>(rparam["punit"]);
    tm.nboot = as<unsigned int>(rparam["nboot"]);
    tm.showtime = as<unsigned int>(rparam["showtime"]);
    tm.warning = as<unsigned int>(rparam["warning"]);

//    // for debug
//    Rprintf("Input param arguments:\n tol=%.4f, nboot=%d, cor_type=%d, test_type=%d, resamp=%d, showtime=%d\n",mm.tol, tm.nboot, tm.corr, tm.test, tm.resamp, tm.showtime);

    NumericMatrix Yr(Ysexp);
    NumericMatrix Xr(Xsexp);
    NumericVector lambda(LamSexp);
    unsigned int nRows = Yr.nrow();
    unsigned int nVars = Yr.ncol();
    unsigned int nParam = Xr.ncol();
    unsigned int nLambda = lambda.size();
    tm.nRows = nRows;
    tm.nVars = nVars;
    tm.nParam = nParam;

    // Rcpp -> gsl
    unsigned int i, j, k;
    tm.smry_lambda = gsl_vector_alloc(nLambda);
    for (i=0; i<nLambda; i++)
        gsl_vector_set(tm.smry_lambda, i, lambda(i));	    
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
    clock_t clk_start, clk_end;
    clk_start = clock();

    // Glm fit
    PoissonGlm pfit(&mm);
    BinGlm lfit(&mm);
    NBinGlm nbfit(&mm);
    glm *glmPtr[3] = { &pfit, &nbfit, &lfit };
    unsigned int mtype = mm.model-1;
    glmPtr[mtype]->regression(Y, X, NULL, NULL);
//    glmPtr[mtype]->display();

    GlmTest myTest(&tm);    
    // Resampling indices
    if ( !Rf_isNumeric(bIDsexp) || !Rf_isMatrix(bIDsexp) ) {
//      Rprintf("Calc bootID on the fly.\n");
    }	   
    else {
        if ( tm.resamp == SCOREBOOT ) {
            NumericMatrix bIDr(bIDsexp);
            tm.nboot = bIDr.nrow();	   
            myTest.bootID = gsl_matrix_alloc(tm.nboot, nRows);
//	    std::copy( bIDr.begin(), bIDr.end(), smry.bootID->data );
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
    }   } 

    // resampling test
    myTest.summary(glmPtr[mtype]);
//    myTest.displaySmry();

    clk_end = clock();
    unsigned long int dif = floor((double)(clk_end - clk_start)/(double)(CLOCKS_PER_SEC));
    unsigned int hours = floor((double)dif/(double)3600);
    unsigned int min = floor((double)(dif-hours*3600)/(double)60);
    unsigned int sec = dif - hours*3600 - min*60;
    if (tm.showtime>=TRUE)
       Rprintf("Time elapsed: %d hr %d min %d sec\n", hours, min, sec);

    // Wrap gsl vectors with Rcpp 
    double multstat, Pmultstat;
    multstat = gsl_matrix_get(myTest.smryStat, 0, 0);
    Pmultstat = gsl_matrix_get(myTest.Psmry, 0, 0);
    NumericVector Vec_aic(myTest.aic, myTest.aic+nVars);
    double *uj = gsl_matrix_ptr(myTest.smryStat, 0, 1);
    double *pj = gsl_matrix_ptr(myTest.Psmry, 0, 1);
    NumericVector Vec_unitmult = wrap(uj, uj+nVars);
    NumericVector Vec_Punitmult = wrap(pj, pj+nVars);
    
    NumericVector Vec_signi(nParam);
    NumericVector Vec_Psigni(nParam);
    NumericMatrix Mat_unitsigni(nParam, nVars);
    NumericMatrix Mat_Punitsigni(nParam, nVars);
    for (i=0; i<nParam; i++){        
        uj=gsl_matrix_ptr (myTest.smryStat, i+1, 0);
	pj=gsl_matrix_ptr (myTest.Psmry, i+1, 0);
        Vec_signi(i) = *(uj); 
        Vec_Psigni(i) = *(pj);
	for (j=0; j<nVars; j++){
	    Mat_unitsigni(i, j) = *(uj+j+1);
	    Mat_Punitsigni(i, j) = *(pj+j+1);
	}
    }

//    std::copy(ptr_statj+nVars, ptr_statj+nVars*nParam-1, 
//              Mat_unitsign.begin());
//    std::copy(ptr_Pstatj+nVars+1, ptr_Pstatj+1+nVars*nParam, 
//              Mat_Punitsign.begin());
//    Rprintf("Done\n.");

    unsigned int nSamp = myTest.nSamp;

    // Rcpp -> R
    List rs = List::create(
	_["multstat" ] = multstat,
	_["Pmultstat"] = Pmultstat,
        _["unitmult" ] = Vec_unitmult,
        _["Punitmult"] = Vec_Punitmult,
	_["signific" ] = Vec_signi,
	_["Psignific"] = Vec_Psigni,
	_["unitsign" ] = Mat_unitsigni,
	_["Punitsign"] = Mat_Punitsigni,
	_["nSamp"    ] = nSamp,
	_["aic"      ] = Vec_aic
    );

    // clear objects
    glmPtr[mtype]->releaseGlm();
    myTest.releaseTest();    
    gsl_matrix_free(Y);
    gsl_matrix_free(X);
    gsl_vector_free(tm.smry_lambda);
    
    return rs;
}

