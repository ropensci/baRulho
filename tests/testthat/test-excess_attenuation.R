test_that("using extended table", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]
  
  X <- set_reference_sounds(X, pb = FALSE)
  
  ea <- excess_attenuation(X = X)
  
  expect_equal(sum(is.na(ea$excess.attenuation)), 9)
  
  expect_equal(nrow(ea), 25)
  
  expect_equal(ncol(ea), 11)
  
  expect_equal(class(ea)[1], "extended_selection_table")
  
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
  
  X <- set_reference_sounds(X, method = 2)
  
  expect_warning(ea <- excess_attenuation(X = X))
  
  expect_equal(sum(is.na(ea$excess.attenuation)), 13)
  
  expect_equal(nrow(ea), 25)
  
  expect_equal(ncol(ea), 11)
  
  expect_equal(class(ea)[1], "data.frame")
  
})
