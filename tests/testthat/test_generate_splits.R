# --------------------------------------------------
# Test Script - Error for generate_subsets Function
# --------------------------------------------------

# Required libraries
library(splitSelect)

# Context of test script
context("Verify generation of all possible splits.")

# There should be an error if we want to compute the IF TS, and no returns are provided
test_that("Error for the \"generate_splits\" function. Unequal number of possible splits.", {
  
  # Number of variables
  p <- 8
  # Number of groups
  G <- 4:5
  expect_equal(sum(nsplit(p,G)), nrow(generate_splits(p, G)))
})