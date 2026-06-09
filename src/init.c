/* Native routine registration for the RefBasedMI package.
 *
 * Registers the five Fortran entry points of the vendored norm2 engine
 * (called from R/norm2.R via .Fortran) and disables dynamic symbol lookup,
 * as required for CRAN. Argument counts match the .Fortran call sites and the
 * subroutine signatures in src/norm2.f90.
 */

#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(norm_em)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(norm_imp_mean)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(norm_imp_rand)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(norm_logpost)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(norm_mcmc)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"norm_em",       (DL_FUNC) &F77_NAME(norm_em),       38},
    {"norm_imp_mean", (DL_FUNC) &F77_NAME(norm_imp_mean), 13},
    {"norm_imp_rand", (DL_FUNC) &F77_NAME(norm_imp_rand), 14},
    {"norm_logpost",  (DL_FUNC) &F77_NAME(norm_logpost),  16},
    {"norm_mcmc",     (DL_FUNC) &F77_NAME(norm_mcmc),     40},
    {NULL, NULL, 0}
};

void R_init_RefBasedMI(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
