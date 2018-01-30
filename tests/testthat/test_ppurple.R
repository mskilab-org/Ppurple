library(Ppurple)
library(data.table)

context('testing Ppurple ops')

## Error in test_files(paths, reporter = reporter, env = env, stop_on_failure = stop_on_failure,  : 
cov = fread(system.file("extdata", "coverage.csv", package = "Ppurple"))
hets = fread(system.file("extdata", "hets.csv", package = "Ppurple"))
segs = fread(system.file("extdata", "segs.csv", package = "Ppurple"))
 

test_that("ppurple() works without seg input", {
  pp = ppurple(cov = cov, hets = hets, verbose = TRUE)
  expect_equal(round(pp[1,]$purity,2), 0.49)
  expect_equal(round(pp[1,]$ploidy,2), 3.86)
})

test_that("ppurple() works with seg input", {
  pp = ppurple(cov = cov, hets = hets, segs = segs, verbose = TRUE)  
  expect_equal(round(pp[1,]$purity,2), 0.49)
  expect_equal(round(pp[1,]$ploidy,2), 3.86)
})

test_that("ppurple() works with gr input on subsample", {
  set.seed(42)
  pp = ppurple(cov = sample(dt2gr(cov), 10000), hets = sample(dt2gr(hets), 10000), segs = segs, verbose = TRUE)
  expect_equal(round(pp[1,]$purity,2), 0.49)
  expect_equal(round(pp[1,]$ploidy,2), 3.86)    
})


