# Regression tests for smoother()'s per-label threshold lookup.

set.seed(1)
locs <- cbind(runif(300), runif(300))

test_that("smoother handles labels with a gap in the numbering", {
  # Cluster 3 absent, cluster 9 present: 8 distinct labels but a max of 9.
  # Looking the threshold up by position returned NA for label 9 and the
  # comparison in the loop errored.
  labels <- sample(c(1, 2, 4, 5, 6, 7, 8, 9), 300, replace = TRUE)
  expect_no_error(smoother(labels, locs = locs, k = 10, verbose = FALSE))
})

test_that("smoother is unchanged when labels are a contiguous 1..n", {
  labels <- sample(1:9, 300, replace = TRUE)
  out <- smoother(labels, locs = locs, k = 10, verbose = FALSE)
  expect_length(out, length(labels))
  expect_true(all(out %in% labels))
})

test_that("explicit named props are respected", {
  labels <- sample(1:3, 300, replace = TRUE)
  props <- c("1" = 1.5, "2" = 1.5, "3" = 1.5)
  expect_identical(smoother(labels, locs = locs, k = 10, props = props,
                            verbose = FALSE), labels)
})

test_that("unnamed props keep positional lookup", {
  labels <- sample(1:3, 300, replace = TRUE)
  expect_identical(smoother(labels, locs = locs, k = 10,
                            props = rep(1.5, 3), verbose = FALSE), labels)
})

test_that("factor-code labels are unaffected", {
  f <- factor(sample(c("CA1", "GCL", "WM"), 300, replace = TRUE))
  expect_no_error(smoother(as.numeric(f), locs = locs, k = 10, verbose = FALSE))
})
