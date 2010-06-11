// Header file for anova and summary
// Author: Yi Wang (yi dot wang at unsw dot edu dot au)
// Last modified: 20-April-2010

#ifndef _RESAMPTEST_H
#define _RESAMPTEST_H

#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h> 
#include <math.h>

// logic
#define TRUE 1
#define FALSE 0
// shrinkage
#define NOSHRINK 0
#define IDENTITY 1
#define SHRINK 2
// test
#define LOGWILK 0
#define HOTELING 1
// resampling
#define CASEBOOT 0
#define RESIBOOT 1
#define SCOREBOOT 2
#define PERMUTE 3
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
#define TOL 1e-10 
#define MaxIter 1000
#define LAMBDA 0.2236
#define NaN -100000
#define MAX_LINE_LENGTH 65536
#define WRAP 4
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))

typedef struct MethodStruc {
    double tol;
    int nboot;
    int corr;
    double shrink_param;
    int test;
    int resamp;
    int reprand;
    int student;
    int punit;
    int rsquare;
} mv_Method;

typedef struct matStruc {
    gsl_matrix *mat;   // hat(X)
    gsl_matrix *SS;
    gsl_matrix *Coef;
    gsl_matrix *Res;
    gsl_matrix *X;
    gsl_matrix *Y;
    double teststat;
} mv_mat;


// io.c
int vector_filesize(FILE *f);
void matrix_filesize(FILE *f, int * row, int * col);
gsl_matrix * load_m(char * file);
gsl_vector * load_v(char * file);
void displaymatrix(gsl_matrix * m, char * name);
void displayvector(gsl_vector * v, char * name);

// calctest.c
int testStatCalc(mv_mat *H0, mv_mat *H1, mv_Method *mmRef, const int ifcalcH1det, double *stat, gsl_vector *statj);
int calcSS(gsl_matrix *Y, mv_mat *Hat, mv_Method *mmRef, const int ifcalcHat, const int ifcalcCoef, const int ifcalcSS);
double calcDet(gsl_matrix *SS);
int is_sym_matrix(const gsl_matrix *mat);
int subX(gsl_matrix *X, gsl_vector *ref, gsl_matrix *Xi);
int calcAdjustP(const int punit, const int nVars, gsl_vector *bStatj, double *sj, double *pj, gsl_permutation *sortid);
int rcalc(gsl_matrix *Res, double, int, gsl_matrix *SS);

// anova.cpp
class AnovaTest
{
public: mv_Method *mmRef;
	gsl_matrix *Yref;
	gsl_matrix *Xref;
	gsl_matrix *inRef;
	size_t nSamp;

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
//       int SimuTest(gsl_matrix *Y); // for simulation updates

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
	size_t nSamp;

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


#endif
