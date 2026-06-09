test_that("all five Fortran routines are registered and loadable", {
  routines <- c("norm_em", "norm_mcmc", "norm_imp_rand", "norm_imp_mean", "norm_logpost")
  for (r in routines) {
    expect_true(is.loaded(r, PACKAGE = "RefBasedMI", type = "Fortran"),
                info = paste("Fortran routine not registered:", r))
  }
})
