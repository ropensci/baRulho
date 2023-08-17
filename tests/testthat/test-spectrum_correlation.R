test_that("using extended table and method 1", {
  data("degradation_est")
  
  X <- degradation_est[degradation_est$sound.files != "master.wav", ]
  
  ec <- spectrum_correlation(X = X, method = 1)
  
  expect_equal(sum(is.na(ec$spectrum.correlation)), 9)
  
  expect_equal(nrow(ec), 25)
  
  expect_equal(ncol(ec), 11)
  
  expect_equal(class(ec)[1], "extended_selection_table")
  
})

test_that("using data frame and method 2", {
  
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[-1])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <- as.data.frame(degradation_est[degradation_est$sound.files != "master.wav", ])
  
  expect_warning(ec <- spectrum_correlation(X = X, method = 2))
  
  expect_equal(sum(is.na(ec$spectrum.correlation)), 9)
  
  expect_equal(nrow(ec), 25)
  
  expect_equal(ncol(ec), 11)
  
  expect_equal(class(ec)[1], "data.frame")
  
})
