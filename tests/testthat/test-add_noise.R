test_that("basic", {
  
  # load example data
  data("test_sounds_est")
  
  # make it a 'by song' extended selection table
  X <- by_element_est(X = test_sounds_est, pb = FALSE)
  #'
  #' # add noise to the first five rows
  X_noise <-
    add_noise(
      X = X[1:5, ],
      mar = 0.2,
      target.snr = 3,
      cores = 1,
      pb = FALSE,
      max.iterations = 100,
      kind = "pink"
    )
  
  expect_true(is_extended_selection_table(X_noise))
  
  expect_equal(sum(!is.na(X_noise$adjusted.snr)), 4)

  expect_equal(sum(X_noise$adjusted.snr, na.rm = TRUE), 11.905124)
  
  expect_equal(nrow(X_noise), 5)
  
  expect_equal(ncol(X_noise), 10)
  
  # test white noise
  #' # add noise to the first five rows
  X_noise <-
    add_noise(
      X = X[1:5, ],
      mar = 0.2,
      target.snr = 3,
      cores = 1,
      pb = FALSE,
      max.iterations = 100,
      kind = "white"
    )
  
  expect_equal(sum(X_noise$adjusted.snr, na.rm = TRUE), 12.107109)
  
})
