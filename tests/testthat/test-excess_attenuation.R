test_that("using extended table", {
  data("degradation_est")
  
  X <- degradation_est[degradation_est$sound.files != "master.wav", ]
  
  ea <- excess_attenuation(X = X, method = 1)
  
  expect_equal(sum(is.na(ea$excess.attenuation)), 9)
  
  expect_equal(nrow(ea), 25)
  
  expect_equal(ncol(ea), 11)
  
  expect_equal(class(ea)[1], "extended_selection_table")
  
})

test_that("using data frame and method 2", {
  
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[-1])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <- as.data.frame(degradation_est[degradation_est$sound.files != "master.wav", ])
  
  expect_warning(ea <- excess_attenuation(X = X, method = 2))
  
  expect_equal(sum(is.na(ea$excess.attenuation)), 13)
  
  expect_equal(nrow(ea), 25)
  
  expect_equal(ncol(ea), 11)
  
  expect_equal(class(ea)[1], "data.frame")
  
})
