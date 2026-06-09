# group_summarize() is the internal Hmisc::summarize() replacement. Its row order
# (grouping variables ascending, first varying slowest) is load-bearing for the
# cumulative case counts that index the imputation loops, so lock it down.

test_that("group_summarize groups, sums, and orders like Hmisc::summarize", {
  gs <- getFromNamespace("group_summarize", "RefBasedMI")
  df <- data.frame(treat = c(2, 1, 2, 1, 1, 2, 1),
                   patt  = c(15, 3, 7, 15, 3, 7, 1),
                   x = 1)
  out <- gs(df$x, by = list(df$treat, df$patt), FUN = sum)

  expect_equal(nrow(out), 5L)
  # ordered by treat then patt, ascending
  expect_equal(out[[1]], c(1, 1, 1, 2, 2))
  expect_equal(out[[2]], c(1, 3, 15, 7, 15))
  # counts per group
  expect_equal(out$value, c(1, 2, 1, 2, 1))
})

test_that("group_summarize rejects a non-list 'by'", {
  gs <- getFromNamespace("group_summarize", "RefBasedMI")
  expect_error(gs(1:3, by = 1:3, FUN = sum), "list")
})
