# Scaling of the projected factors. These use synthetic objects so they run
# without fetching any data.

mk_spe <- function(genes, ncell, seed) {
  set.seed(seed)
  m <- matrix(abs(rnorm(length(genes) * ncell)), length(genes), ncell,
              dimnames = list(genes, paste0("c", seq_len(ncell))))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(logcounts = m),
    spatialCoords = matrix(runif(2 * ncell), ncell, 2,
                           dimnames = list(NULL, c("x", "y"))))
  SummarizedExperiment::rowData(spe)$gene_name <- genes
  spe
}

mk_nmf <- function(ngene, k, nspot = 60, seed = 1) {
  set.seed(seed)
  w <- matrix(abs(rnorm(ngene * k)), ngene, k)
  # h is k x n_spots with rows summing to 1, as singlet returns it. The scale
  # estimated from a small-panel target is expressed per source spot, so
  # ncol(h) is needed and h can no longer be NULL.
  h <- matrix(abs(rnorm(k * nspot)), k, nspot)
  list(w = sweep(w, 2, colSums(w), "/"), d = 10 * 2^seq_len(k),
       h = sweep(h, 1, rowSums(h), "/"))
}

genes <- paste0("g", 1:400)
k     <- 5
src   <- mk_spe(genes, 60, 1)
nmf   <- mk_nmf(length(genes), k)

test_that("the per-factor scale is applied before the transpose", {
  # proj is k x n, so the division must happen while factors are in rows.
  # t(proj)/d recycles down columns and mis-scales all but 1 in k entries.
  tgt  <- mk_spe(genes, 30, 2)
  pr   <- project_raw(src, tgt, "logcounts", nmf)
  fac  <- project_factors(src, tgt, "logcounts", nmf)
  expect_equal(unname(fac), unname(t(pr$proj / nmf$d)))
  expect_false(isTRUE(all.equal(unname(fac), unname(t(pr$proj) / nmf$d))))
})

test_that("a target sharing the source's genes uses the source d", {
  tgt <- mk_spe(genes, 30, 2)
  pr  <- project_raw(src, tgt, "logcounts", nmf)
  expect_equal(target_factor_scale(pr$proj, nmf, pr$n_shared), nmf$d)
})

test_that("a small-panel target uses a scale estimated from itself", {
  tgt <- mk_spe(genes[1:40], 30, 3)          # 10% of the source's genes
  pr  <- project_raw(src, tgt, "logcounts", nmf)
  sc  <- suppressWarnings(target_factor_scale(pr$proj, nmf, pr$n_shared))
  expect_false(isTRUE(all.equal(sc, nmf$d)))
  expect_equal(length(sc), k)
})

test_that("factors with no mass in the target fall back to the source d", {
  proj <- matrix(runif(k * 20), k, 20)
  proj[2, ] <- 0                              # factor 2 absent from the target
  expect_warning(sc <- target_factor_scale(proj, nmf, n_shared = 40))
  expect_equal(sc[2], nmf$d[2])
  expect_equal(sc[-2], (rowMeans(proj) * ncol(nmf$h))[-2])
})

test_that("a supplied pooled scale overrides the per-target estimate", {
  proj   <- matrix(runif(k * 20), k, 20)
  pooled <- rowSums(proj) * 7
  expect_equal(target_factor_scale(proj, nmf, 40, pooled_d_target = pooled), pooled)
})

test_that("the target-estimated scale does not depend on the number of cells", {
  # The reason for using the mean rather than the total: doubling the cells in a
  # section must not change the scale, or the predictions depend on section size.
  proj  <- matrix(runif(k * 20), k, 20)
  twice <- cbind(proj, proj)                  # same section, sampled twice over
  expect_equal(target_factor_scale(proj,  nmf, n_shared = 40),
               target_factor_scale(twice, nmf, n_shared = 40))
})

test_that("a high-overlap target still uses the source d unchanged", {
  # Visium -> Visium needs no size correction; that arm must be untouched.
  proj <- matrix(runif(k * 20), k, 20)
  expect_equal(target_factor_scale(proj, nmf, n_shared = nrow(nmf$w)), nmf$d)
})
