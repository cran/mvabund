// Interface between R and glm.cpp (Rcpp API >= 0.8.6)
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// 20-April-2011

#include <Rcpp.h>
extern "C"{
#include "resampTest.h"
//#include "time.h"
}

RcppExport SEXP RtoGlm(SEXP params, SEXP Ysexp, SEXP Xsexp, SEXP Osexp) 
{
    using namespace Rcpp;

    // Get parameters in params.
    List rparam(params);
    // pass parameters
    reg_Method mm;	
    mm.tol = as<double>(rparam["tol"]);
    mm.model = as<unsigned int>(rparam["regression"]);
    mm.estiMethod = as<unsigned int>(rparam["estimation"]);
    mm.varStab = as<unsigned int>(rparam["stablizer"]);
    mm.n = as<unsigned int>(rparam["n"]);
    mm.maxiter = as<unsigned int>(rparam["maxiter"]);
    mm.maxiter2 = as<unsigned int>(rparam["maxiter2"]);
    mm.warning = as<unsigned int>(rparam["warning"]);
// for debug
//    Rprintf("tol=%g, model=%d, estiMethod=%d, varStab=%d\n", mm.tol, mm.model, mm.estiMethod, mm.varStab);

    NumericMatrix Yr(Ysexp);
    NumericMatrix Xr(Xsexp);
    NumericMatrix Or(Osexp);
    unsigned int nRows = Yr.nrow();
    unsigned int nVars = Yr.ncol();
    unsigned int nParam = Xr.ncol();
// for debug
//    Rprintf("nRows=%d, nVars=%d, nParam=%d\n", nRows, nVars, nParam);

    // Rcpp -> gsl
    unsigned int i, j, k;
    gsl_matrix *X = gsl_matrix_alloc(nRows, nParam);        
    gsl_matrix *Y = gsl_matrix_alloc(nRows, nVars);
    gsl_matrix *O = gsl_matrix_alloc(nRows, nVars);

//  Must be careful about using std::copy for matrix. The following direct
//  use is not doing right - row elements are copied to columns. Need to
//  fix it later on. 
//    std::copy( Yr.begin(), Yr.end(), Y->data );
//    std::copy( Xr.begin(), Xr.end(), X->data );
//    std::copy( INr.begin(), INr.end(), isXvarIn->data );

//    Rprintf("Y passed to C\n");
    for (i=0; i<nRows; i++){
        for (j=0; j<nVars; j++) {
            gsl_matrix_set(Y, i, j, Yr(i, j));
            gsl_matrix_set(O, i, j, Or(i, j));
//            Rprintf("%.2f ", gsl_matrix_get(Y, i, j));
        }
//        Rprintf("\t");
//
        for (k=0; k<nParam; k++){
	    gsl_matrix_set(X, i, k, Xr(i, k));
//	    Rprintf("%.2f ", gsl_matrix_get(X, i, k));
	}
//	Rprintf("\n");
    }
       
    // do stuff	
//    clock_t clk_start, clk_end;
//    clk_start = clock();

    PoissonGlm pfit(&mm);
    BinGlm lfit(&mm);
    NBinGlm nbfit(&mm);
    glm *glmPtr[3] = { &pfit, &nbfit, &lfit };
    unsigned int mtype = mm.model-1;
    glmPtr[mtype]->regression(Y, X, O, NULL);
//    glmPtr[mtype]->display();
	
//    clk_end = clock();
//    unsigned long int dif = floor((double)(clk_end - clk_start)/(double)(CLOCKS_PER_SEC));
//    unsigned int hours = floor((double)(dif/(double)3600));
//    unsigned int min = floor((double)(dif%3600)/(double)60);
//    unsigned int sec = dif%60;
//    Rprintf("Time elapsed: %d hr %d min %d sec (%d seconds)\n", hours, min, sec, dif);

    // Wrap the glm object with Rcpp 
    NumericVector theta(glmPtr[mtype]->theta, glmPtr[mtype]->theta+nVars);
    NumericVector ll(glmPtr[mtype]->ll, glmPtr[mtype]->ll+nVars);
    NumericVector aic(glmPtr[mtype]->aic, glmPtr[mtype]->aic+nVars);
    NumericVector dev(glmPtr[mtype]->dev, glmPtr[mtype]->dev+nVars);
    NumericVector iterconv(glmPtr[mtype]->iterconv, glmPtr[mtype]->iterconv+nVars);
    
    NumericMatrix Beta(nParam, nVars);
    NumericMatrix varBeta(nParam, nVars);
    NumericMatrix Mu(nRows, nVars);
    NumericMatrix Eta(nRows, nVars);
    NumericMatrix Vars(nRows, nVars);
    NumericMatrix wHalf(nRows, nVars);
    NumericMatrix Res(nRows, nVars);
    NumericMatrix PitRes(nRows, nVars);
    NumericMatrix sqrt1_Hii(nRows, nVars);

    for (i=0; i<nParam; i++)
    for (j=0; j<nVars; j++) {
        Beta(i, j) = gsl_matrix_get(glmPtr[mtype]->Beta, i, j);
        varBeta(i, j) = gsl_matrix_get(glmPtr[mtype]->varBeta, i, j);
    }
//    Rprintf("Residuals calculated by C\n");
    for (i=0; i<nRows; i++) //{
    for (j=0; j<nVars; j++){
	Mu(i, j) = gsl_matrix_get(glmPtr[mtype]->Mu, i, j);        
	Eta(i, j) = gsl_matrix_get(glmPtr[mtype]->Eta, i, j);        
	Vars(i, j) = gsl_matrix_get(glmPtr[mtype]->Var, i, j);        
	wHalf(i, j) = gsl_matrix_get(glmPtr[mtype]->wHalf, i, j);        
	Res(i, j) = gsl_matrix_get(glmPtr[mtype]->Res, i, j);        
	PitRes(i, j) = gsl_matrix_get(glmPtr[mtype]->PitRes, i, j);        
	sqrt1_Hii(i, j) = gsl_matrix_get(glmPtr[mtype]->sqrt1_Hii, i, j);        
//        Rprintf("%d ", (int)gsl_matrix_get(Y, i, j));
    }
//        Rprintf("\t");
//    }       

//    double *uj = gsl_matrix_ptr(anova.statj, 0, 0);
//    double *pj = gsl_matrix_ptr(anova.Pstatj, 0, 0);
//    std::copy(uj, uj+nVars*(nModels-1), Mat_statj.begin());
//    std::copy(pj, pj+nVars*(nModels-1), Mat_Pstatj.begin());

    // Rcpp -> R
    List rs = List::create(
         _["coefficients"] = Beta,
         _["var.coefficients"] = varBeta,
         _["fitted.values"] = Mu,
         _["linear.predictor"] = Eta,
	 _["residuals"] = Res,
	 _["PIT.residuals"] = PitRes,
	 _["sqrt.1_Hii"] = sqrt1_Hii,
         _["var.estimator"] = Vars,
	 _["sqrt.weight"] = wHalf,
	 _["theta"] = theta,
	 _["two.loglike"] = ll,
	 _["deviance"] = dev,
	 _["aic"] = aic,
	 _["iter"] = iterconv
    );

    // clear objects
    glmPtr[mtype]->releaseGlm();
    gsl_matrix_free(Y);
    gsl_matrix_free(X);
    gsl_matrix_free(O);

    return rs;
}

