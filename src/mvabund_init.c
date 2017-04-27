#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mvabund_RtoAnovaCpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mvabund_RtoGlm(SEXP, SEXP, SEXP, SEXP);
extern SEXP mvabund_RtoGlmAnova(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mvabund_RtoGlmSmry(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mvabund_RtoSmryCpp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mvabund_RtoAnovaCpp", (DL_FUNC) &mvabund_RtoAnovaCpp, 5},
    {"mvabund_RtoGlm",      (DL_FUNC) &mvabund_RtoGlm,      4},
    {"mvabund_RtoGlmAnova", (DL_FUNC) &mvabund_RtoGlmAnova, 8},
    {"mvabund_RtoGlmSmry",  (DL_FUNC) &mvabund_RtoGlmSmry,  7},
    {"mvabund_RtoSmryCpp",  (DL_FUNC) &mvabund_RtoSmryCpp,  4},
    {NULL, NULL, 0}
};

void R_init_mvabund(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
//    R_forceSymbols(dll, TRUE);
}
