# -----------------------------------------
# Test Script - Error nsplit function
# -----------------------------------------

# Required libraries
library(splitSelect)

# Context of test script
context("Verify computation of number of possible splits.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("Error for the \"nsplit\" function. Unequal number of possible splits.", {
  
  # Number of variables
  p <- 10
  # Number of groups
  G <- 4:5
  expect_equal(sum(nsplit(p,G)), nsplit(p, G[1]) + nsplit(p, G[2]))
})