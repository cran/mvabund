// Interface between R and summary.cpp
//
// Author: Yi Wang (yi dot wang at computer dot org)
// 16-Nov-2009

#include "Rcpp.h"
extern "C"{
#include "resampTest.h"
#include "time.h"
}

RcppExport SEXP RtoSmryCpp(SEXP params, SEXP Ysexp, SEXP Xsexp,  
                             SEXP bIDsexp )//, 
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
        int nRows = Yr.dim1();
	int nVars = Yr.dim2();
	int nParam = Xr.dim2();

        // Rcpp -> gsl
	int i, j, k;
	gsl_matrix *X = gsl_matrix_alloc(nRows, nParam);        
	gsl_matrix *Y = gsl_matrix_alloc(nRows, nVars);
        for (i=0; i<nRows; i++){
	    for (j=0; j<nVars; j++) 
	        gsl_matrix_set(Y, i, j, (int)Yr(i, j));
            for (k=0; k<nParam; k++)
	        gsl_matrix_set(X, i, k, (double)Xr(i, k));
         }
       
        // do stuff	
	clock_t clk_start, clk_end;
	clk_start = clock();

	// initialize summary class
	Summary smry(&mm, Y, X);
	
	// Resampling indices
        if ( !Rf_isNumeric(bIDsexp) || !Rf_isMatrix(bIDsexp) ) {
//	   Rprintf("Calc bootID on the fly.\n");
        }	   
	else {
	   if ( mm.resamp == SCOREBOOT ) {
	      RcppMatrixView<double> bIDr(bIDsexp);
              mm.nboot = bIDr.dim1();	   
	      smry.bootID = gsl_matrix_alloc(mm.nboot, nRows);
	      for (i=0; i<mm.nboot; i++)
	      for (j=0; j<nRows; j++)
                  gsl_matrix_set(smry.bootID, i, j, bIDr(i, j));
	   }
	   else{
	      RcppMatrixView<int> bIDr(bIDsexp);
              mm.nboot = bIDr.dim1();	   
	      smry.bootID = gsl_matrix_alloc(mm.nboot, nRows);
	      for (i=0; i<mm.nboot; i++)
              for (j=0; j<nRows; j++)
                  gsl_matrix_set(smry.bootID, i, j, bIDr(i, j)-1);
	}   } 

       // resampling test
       smry.resampTest();
       //smry.display();

	clk_end = clock();
	float dif = (float)(clk_end - clk_start)/(float)(CLOCKS_PER_SEC);
//	Rprintf("Time elapsed: %.4f seconds\n", dif);

        // gsl -> Rcpp
	RcppVector<double> Vec_unitmult(nVars);
	RcppVector<double> Vec_Punitmult(nVars);
	RcppVector<double> Vec_signific(nParam);
	RcppVector<double> Vec_Psignific(nParam);
	RcppMatrix<double> Mat_unitsign(nParam, nVars);
	RcppMatrix<double> Mat_Punitsign(nParam, nVars);

        double multstat = smry.multstat[0];
	double Pmultstat = smry.Pmultstat[0];
	double *uj=gsl_matrix_ptr (smry.unitstat, 0, 0);
	double *pj=gsl_matrix_ptr (smry.Punitstat, 0, 0);
	for (j=0; j<nVars; j++){
	    Vec_unitmult(j)=*(uj+j);
	    Vec_Punitmult(j)=*(pj+j);
	}
        for (i=0; i<nParam; i++){
	    Vec_signific(i) = smry.multstat[i+1];
	    Vec_Psignific(i) = smry.Pmultstat[i+1];
	    uj=gsl_matrix_ptr (smry.unitstat, i+1, 0);
	    pj=gsl_matrix_ptr (smry.Punitstat, i+1, 0);
	    for (j=0; j<nVars; j++){
	        Mat_unitsign(i, j)=*(uj+j);
		Mat_Punitsign(i, j)=*(pj+j);
	    }
	}

        int nbootsrun = smry.nSamp;

        // clear objects
        smry.releaseSummary(); 
        gsl_matrix_free(Y);
	gsl_matrix_free(X);

	// Rcpp -> R
        RcppResultSet rs;
	rs.add("multstat", multstat);
	rs.add("Pmultstat", Pmultstat);
	rs.add("unitmult", Vec_unitmult);
        rs.add("Punitmult", Vec_Punitmult);
	rs.add("signific", Vec_signific);
	rs.add("Psignific", Vec_Psignific);
	rs.add("unitsign", Mat_unitsign);
	rs.add("Punitsign", Mat_Punitsign);
	rs.add("nSamp", nbootsrun);
	rs.add("R2", smry.R2);
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

