test_that("envelope blur ratio", {
  # load example data
  data("test_sounds_est")
  
  # add reference to X
  X <-
    set_reference_sounds(X = test_sounds_est[test_sounds_est$sound.id == test_sounds_est$sound.id[6],])
  
  unlink(list.files(
    path = tempdir(),
    pattern = "^blur_ratio",
    full.names = TRUE
  ))
  
  # create plots
  plot_blur_ratio(X = X,
                  dest.path = tempdir(),
                  pb = FALSE)
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^blur_ratio",
               full.names = TRUE)
  
  expect_length(fls, 4)
  
  unlink(fls)
})


test_that("only rows with a valid reference are processed (not every row)", {
  data("test_sounds_est")

  X <- set_reference_sounds(X = test_sounds_est, pb = FALSE)

  n_valid_refs <- sum(!is.na(X$reference))

  # sanity check: this data set actually has rows without a reference
  # (ambient noise, or the reference itself), so this test is meaningful
  expect_true(n_valid_refs < nrow(X))

  call_count <- 0
  testthat::local_mocked_bindings(.plot_blur = function(...) {
    call_count <<- call_count + 1
    NULL
  })

  plot_blur_ratio(X = X, dest.path = tempdir(), pb = FALSE)

  expect_equal(call_count, n_valid_refs)

})

test_that("spectrum blur ratio", {
  # load example data
  data("test_sounds_est")
  
  # add reference to X
  X <-
    set_reference_sounds(X = test_sounds_est[test_sounds_est$sound.id == test_sounds_est$sound.id[6],])
  
  unlink(list.files(
    path = tempdir(),
    pattern = "^blur_ratio",
    full.names = TRUE
  ))
  
  # create plots
  plot_blur_ratio(
    X = X,
    dest.path = tempdir(),
    pb = FALSE,
    type = "spectrum"
  )
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^spectrum_blur_ratio",
               full.names = TRUE)
  
  expect_length(fls, 4)
  
  unlink(fls)
})
