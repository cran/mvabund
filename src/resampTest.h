// Main header file
// Author: Yi Wang (yi dot wang at unsw dot edu dot au
// 16-Jun-2011

#ifndef _RESAMPTEST_H
#define _RESAMPTEST_H

// original resampTest.h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h> 
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_errno.h>

// rmv.h
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_gamma.h>

// return status
#define SUCCESS 0
#define FAILED 1
#define CannotOpenFile 2  

// logic
#define TRUE 1
#define FALSE 0
// model
#define LM 0
#define POISSON 1
#define NB 2
#define LOGIT 3
// shrinkage
#define NOSHRINK 0
#define IDENTITY 1
#define SHRINK 2
// test
#define LOGWILK 0
#define HOTELING 1
#define WALD 2
#define SCORE 3
#define LR 4
// estiMethod
#define NEWTON 0
#define CHI2 1
#define FISHER 2
// infoMatrix
#define OIM 0
#define EIM 1
// resampling
#define CASEBOOT 0
#define RESIBOOT 1
#define SCOREBOOT 2
#define PERMUTE 3
#define FREEPERM 4
#define MONTECARLO 5
#define RESIZ 6 
#define SCOREZ 7
// p-value adjustment
#define NONE 0
#define UNADJUST 1
#define FREESTEP 2
#define SINGLESTEP 3
#define STEPUP 4
#define NOMONO 5
// R-squared
#define HOOPER 0
#define VECTOR 1
// others
#define TOL 1e-5 
#define MAXITER 999 
#define LAMBDA 0.8  // no shrinkage 
#define NaN -100000
#define MAX_LINE_LENGTH 65536
#define WRAP 4
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))
#define ABS(a) MAX(a, -a)  
//simulations
#define NSIMU 1000
#define NQ 6
#define RHO 0.5
#define ALFA 0.05
#define MAXLINELEN 256
//
typedef struct MethodStruc {
    // hypo test methods
    int nboot;
    int corr;
    int test;
    int resamp;
    int reprand;
    int student;
    int punit;
    int rsquare;
    // numeric
    double shrink_param;
    gsl_vector *smry_lambda;
    gsl_vector *anova_lambda;
    double tol;
} mv_Method;

// used for manylm only
typedef struct matStruc {
    gsl_matrix *mat;   // hat(X)
    gsl_matrix *SS;
    gsl_matrix *Coef;
    gsl_matrix *Res;
    gsl_matrix *X;
    gsl_matrix *Y;
    double teststat;
} mv_mat;

typedef struct GroupMatrix{
    gsl_matrix *matrix;
} GrpMat;

// ManyGlm related
typedef struct RegressionMethod{
    // regression methods
    int model;
    int varStab;
    int estiMethod;
    double tol;
} reg_Method;

// ManyLM related
// anova.cpp
class AnovaTest {
public: mv_Method *mmRef;
	gsl_matrix *Yref;
	gsl_matrix *Xref;
	gsl_matrix *inRef;
	int nSamp;

	double *multstat;
	double *Pmultstat;
	gsl_matrix *statj;
	gsl_matrix *Pstatj;
	int *dfDiff;
	gsl_matrix *bootID;	 

       // Methods
       AnovaTest(mv_Method *, gsl_matrix *, gsl_matrix *, gsl_matrix *isXvarIn);
       virtual ~AnovaTest();
       int resampTest(void); 
       void releaseTest(void);
       void display(void);

private: mv_mat *Hats;
         gsl_permutation **sortid;
	 gsl_vector *bStatj;
	 double bMultStat;
	 size_t nModels, nRows, nVars, nParam;

         // Methods
//         int getBootID(void); // done in R
	 int anovacase(gsl_matrix *bY, gsl_matrix *bX);
         int anovaresi(gsl_matrix *bY, const size_t p);
};

// summary.cpp
class Summary
{
public: mv_Method *mmRef;
	gsl_matrix *Yref;
	gsl_matrix *Xref;
	int nSamp;

        double R2;	
	double *multstat;
	double *Pmultstat;
	gsl_matrix *unitstat;
	gsl_matrix *Punitstat;
	gsl_matrix *bootID;	 

       // Methods
       Summary(mv_Method *, gsl_matrix *, gsl_matrix *);
       virtual ~Summary();
       int resampTest(void); 
       void releaseSummary(void);
       void display(void);

private: mv_mat *Hats;
	 gsl_permutation **sortid;
	 size_t nRows, nVars, nParam;
	 double *bMultStat;
	 gsl_matrix *bUnitStat;

	 // Methods
         int calcR2(void);
//         int getBootID(void); // done in R
	 int smrycase(gsl_matrix *bY, gsl_matrix *bX);
         int smryresi(gsl_matrix *bY);

};

// glm base class
class glm
{
   public:
           glm(const reg_Method *mm);
	   virtual ~glm();	   
	   void initialGlm(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O);
	   void releaseGlm(void);
	   virtual int regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O)=0;
	   virtual int EstIRLS(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O, double *) = 0;
	   int copyGlm(glm *src);
           void display(void); 	  
	   double getMij( double eij ) 
	          {return invLink(eij);}
	   double getWij( double mij, double phi ) 
	          { return weifunc(mij, phi); }
	   double getVij( double mij, double phi ) 
	          { return varfunc(mij, phi); }
	   double getLij( double yij, double mij, double phi ) 
	          { return llfunc(yij, mij, phi); }
           double getDij( double yij, double mij, double phi)
	          { return devfunc(yij, mij, phi); }

	   void InitResampGlm(gsl_matrix *,gsl_matrix *,gsl_matrix *,glm *,gsl_matrix *);

           // input arguments
           const reg_Method *mmRef;
           gsl_matrix *Yref;
	   gsl_matrix *Xref;
	   gsl_matrix *Oref;
	   // return properties
	   gsl_matrix *Beta;
	   gsl_matrix *Mu;
	   gsl_matrix *Eta;
	   gsl_matrix *Res;
	   gsl_matrix *Var;
	   gsl_matrix *wHalf;
	   gsl_matrix *sqrt1_Hii;

           int rdf;
	   double *phi, *ll, *dev, *aic;
	   int *iterconv;  
	   double mintol, lTol;
           size_t maxiter, nRows, nVars, nParams;
   private: 
  	   // abstract 	
	   virtual double link(double) const=0;	   
	   virtual double invLink(double) const=0;	   
	   virtual double rcpLinkDash(double) const=0;	   
	   virtual double weifunc(double, double) const=0;   
	   virtual double varfunc(double, double) const=0;	   
	   virtual double llfunc(double, double, double) const=0;	   
	   virtual double devfunc(double, double, double) const=0;	   
};

// poisson regression
class PoissonGlm : public glm 
{
    public: // public functions
           PoissonGlm(const reg_Method *);
	   virtual ~PoissonGlm();
	   virtual int regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O) 
	        { return EstIRLS ( Y, X, O, NULL); }
	   int EstIRLS( gsl_matrix *, gsl_matrix *, gsl_matrix *, double * );
	   int betaEst( size_t id, size_t iter, double *tol, double a );
	   double getDisper( size_t id ) const;

    private: 
           // Log-link and property functions
           double link(double mui) const 
                { return log(mui); } 
           double invLink(double etai) const 
	        { return exp(etai); }
	   double rcpLinkDash(double mui) const
                { return mui; }
	   double weifunc(double mui, double a) const
	        { return mui; }
	   double varfunc(double mui, double a) const
	        { return mui; }
	   double llfunc(double yi, double mui, double a) const
	        { if (yi==0) return 2*(-mui); else
		  return 2*(yi*log(mui) - mui - gsl_sf_lngamma(yi+1)); }
	   double devfunc(double yi, double mui, double a) const
	        { return 2*(yi*log(MAX(yi, mintol)/mui)-yi+mui); }
};

// logistic regression 
class LogiGlm : public PoissonGlm
{
    public: // public functions
           LogiGlm(const reg_Method *);
	   virtual ~LogiGlm();

    private: // logit link and property functions 
           double link(double mui) const
	       { return log(mui/(1-mui)); }
           double invLink(double etai) const 
	        { return exp(etai)/(1+exp(etai)); }
	   double rcpLinkDash(double mui) const
                { return mui*(1-mui); }
	   double weifunc(double mui, double a) const
	        { return mui*(1-mui); }
	   double varfunc(double mui, double a) const
	        { return mui*(1-mui); }
	   double llfunc(double yi, double mui, double a) const
	        { return (yi * (log(mui)-log(1-mui)) + log(1-mui)); }
	   double devfunc(double yi, double mui, double a) const
	        { return 2*((yi<=mintol)?(-log(1-MAX(mui, 1e-5))):(-log(MIN(mui, 1-1e-5)))); }
};

// negative binomial regression
class NBinGlm : public PoissonGlm
{
    public:
           NBinGlm(const reg_Method *mm);
	   virtual ~NBinGlm();
	   virtual int regression(gsl_matrix *Y, gsl_matrix *X, gsl_matrix *O)
	      { return nbinfit ( Y, X, O); }
	   int nbinfit(gsl_matrix *, gsl_matrix *, gsl_matrix *);
    private:
	   double weifunc(double mui, double a) const
	        { return (mui/(1+a*mui)); }
	   double varfunc(double mui, double a) const
	        { return (mui*(1+a*mui)); }
	   double llfunc(double yi, double mui, double a) const;
	   double devfunc(double yi, double mui, double a) const
                { return 2*(yi*log(MAX(yi, mintol)/mui)-(yi+1/MAX(a,mintol))*log((1+yi*a)/(1+mui*a))); }

	   int getfAfAdash (double a, size_t id, double *fA, double *fAdash);
};


// base test class
class GlmTest
{
    public: 
            const mv_Method *tm;
            glm *fit;	   
	    gsl_matrix *Xin;
	    gsl_matrix *smryStat, *Psmry;
	    gsl_matrix *anovaStat, *Panova;
	    gsl_matrix *bootID;
            int nSamp;	    
            double *aic;
            int *dfDiff;

            // methods
            GlmTest(const mv_Method *tm, glm *fit);
            virtual ~GlmTest(void);
            void releaseTest(void);	

	    int summary(void);
	    void displaySmry(void);

	    int anova(gsl_matrix *);
            void displayAnova(void);
            
    private:
	    int getBootID(void);

	    int geeCalc(glm *PtrAlt, glm *PtrNull, gsl_matrix *, gsl_matrix *);
	    int geeWald(glm *PtrAlt, gsl_matrix *, gsl_matrix *);
	    int geeScore(glm *PtrAlt, glm *PtrNULL, gsl_matrix *, gsl_matrix *);
	    int geeLR(glm *PtrAlt, glm *PtrNULL, gsl_matrix *);
	    int subGeeWald(glm *, gsl_matrix *, gsl_vector *, gsl_matrix *);
	    int subGeeScore(glm *, glm *, gsl_vector *, gsl_matrix *);
	    int subGeeLR(glm *PtrAlt, glm *PtrNull, gsl_vector *teststat);

            int resampData(glm *, gsl_matrix *, GrpMat *, GrpMat *, size_t i ); // summary
	    int resampNonCase(glm *, gsl_matrix *, size_t i);
	    int resampAnovaCase(glm *, gsl_matrix *, gsl_matrix *, gsl_matrix *, size_t i);
	    int setMonteCarlo(glm *model, gsl_matrix *, gsl_matrix *);

            // intermediate data
	    size_t nRows, nVars, nParam, nModels;

	    // the following used in resampling
            gsl_rng *rnd;
            size_t *permid;   // only useful in permutation test
            double lambda;    // intermediate shrinkage parameter
                
            // the following are used in geeCalc
	    gsl_matrix *L;  // only useful in Wald test
	    gsl_matrix *Rlambda, *bRlambda; 
	    gsl_matrix *Wj;
	    gsl_matrix *XBeta, *Sigma; // used in monte carlo simulation
	    //double *mr, *sr; // mean and variance of model residuals
//	    gsl_vector *mr;

	    GrpMat *GrpHii; // group Hii
	    GrpMat *GrpXs;  // group X0
	    GrpMat *GrpOs;  // group offset

};

// io.cpp - input/output functions
int vector_filesize(FILE *f);
void matrix_filesize(FILE *f, int * row, int * col);
gsl_matrix * load_m(const char * file);
gsl_vector * load_v(const char * file);
void displaymatrix(gsl_matrix * m, const char * name);
void displayvector(gsl_vector * v, const char * name);

// calctest.c - utility functions
double calcDet(gsl_matrix *SS);
int is_sym_matrix(const gsl_matrix *mat);
int subX(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subX1(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subX2(gsl_matrix *X, int id, gsl_matrix *Xi);
int subXrow(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subXrow2(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int subXrow1(gsl_matrix *X, gsl_vector *ref0, gsl_vector *ref1, gsl_matrix *Xi);
int GetR(gsl_matrix *Res, size_t corr, gsl_vector *glmshrink, size_t k, gsl_matrix *R); 
int subtractMean(gsl_matrix *dat);
// calctest.c - manylm related functions
int testStatCalc(mv_mat *H0, mv_mat *H1, mv_Method *mmRef, const int ifcalcH1det, double *stat, gsl_vector *statj);
int calcSS(gsl_matrix *Y, mv_mat *Hat, mv_Method *mmRef, const int ifcalcHat, const int ifcalcCoef, const int ifcalcSS);
int calcAdjustP(const int punit, const int nVars, double *bj, double *sj, double *pj, gsl_permutation *sortid);
int reinforceP(double *p, size_t nVars, gsl_permutation *sortid);
int getHat(gsl_matrix *X, gsl_matrix *W, gsl_matrix *Hat);
int invLSQ(gsl_matrix *A, gsl_vector *b, gsl_vector *x);
int rcalc(gsl_matrix *Res, double, int, gsl_matrix *SS);

// simutility.cpp - functions used in simulation tests 
int GetMean(gsl_matrix *X, gsl_matrix *Y, gsl_matrix *Mu);
int GetPdstbtion(double *p, size_t nVars, size_t *isH0var, size_t *cnt, size_t *cntfwe);
//int GetCov (gsl_matrix *Mu, gsl_matrix *Y, size_t AR1MAT, gsl_matrix *Sigma);
//int GetMeanCov(gsl_matrix *X, gsl_matrix *Y, mv_Method *mm, size_t AR1MAT, gsl_matrix *Mu, gsl_matrix *Sigma);

// rnd.c - functions to generate random numbers from multivariate (normal) distributions
// MVN random number generator
int rmvnorm(const gsl_rng *, const int, const gsl_matrix *, gsl_vector *);
// MVN with positive-semi definite covariance matrix
int semirmvnorm(const gsl_rng *, const int, const gsl_matrix *, gsl_vector *);

#endif
