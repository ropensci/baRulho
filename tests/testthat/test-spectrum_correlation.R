test_that("using extended table and method 1", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]
  
  X <- set_reference_sounds(X, pb = FALSE)
  
  ec <- spectrum_correlation(X = X, pb = FALSE)
  
  expect_equal(sum(is.na(ec$spectrum.correlation)), 9)
  
  expect_equal(nrow(ec), 25)
  
  expect_equal(ncol(ec), 11)
  
  
  expect_equal(class(ec)[1], "extended_selection_table")
  
})

test_that("using data frame and method 2", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[-1])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <-
    as.data.frame(test_sounds_est[test_sounds_est$sound.files != "master.wav",])
  
  X <- set_reference_sounds(X, method = 2, pb = FALSE)
  
  expect_warning(ec <- spectrum_correlation(X = X,
                                            cores = 1,
                                            pb = FALSE,
                                            cor.method = "pearson",
                                            spec.smooth = 5,
                                            path = td,
                                            n.bins = 100))
  
  expect_equal(sum(is.na(ec$spectrum.correlation)), 13)
  
  expect_equal(nrow(ec), 25)
  
  expect_equal(ncol(ec), 11)
  
  expect_equal(class(ec)[1], "data.frame")
  
})
