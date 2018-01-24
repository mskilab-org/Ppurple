sessionInfo()
suppressPackageStartupMessages(library(Ppurple))

cov = fread(system.file("extdata", "coverage.csv", package = "Ppurple"))
hets = fread(system.file("extdata", "hets.csv", package = "Ppurple"))
segs = fread(system.file("extdata", "segs.csv", package = "Ppurple"))

test_that("ppurple", {
  pp = ppurple(cov = cov, hets = hets, verbose = TRUE)
  expect_equal(round(pp[1,]$purity,2), 0.49)
  expect_equal(round(pp[1,]$ploidy,2), 3.86)

  pp = ppurple(cov = cov, hets = hets, segs = segs, verbose = TRUE)  
  expect_equal(round(pp[1,]$purity,2), 0.49)
  expect_equal(round(pp[1,]$ploidy,2), 3.86)

  set.seed(42)
  pp = ppurple(cov = sample(dt2gr(cov), 1e4), hets = sample(dt2gr(hets), 1e4), segs = segs, verbose = TRUE)  
  expect_equal(round(pp[1,]$purity,2), 0.49)
  expect_equal(round(pp[1,]$ploidy,2), 3.90)    
})

