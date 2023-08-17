test_that("using extended table method 1", {
  data("degradation_est")
  
  X <- degradation_est[degradation_est$sound.files != "master.wav", ]
  
  xc <- spcc(X = X, method = 1)
  
  expect_equal(sum(is.na(xc$cross.correlation)), 9)
  
  expect_equal(nrow(xc), 25)
  
  expect_equal(ncol(xc), 11)
  
  expect_equal(class(xc)[1], "extended_selection_table")
  
})

test_that("using data frame", {
  
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[-1])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)

  X <- as.data.frame(degradation_est[degradation_est$sound.files != "master.wav", ])
  
  expect_warning(xc <- spcc(X = X, method = 2))
  
  expect_equal(sum(is.na(xc$cross.correlation)), 9)
  
  expect_equal(nrow(xc), 25)
  
  expect_equal(ncol(xc), 11)
  
  expect_equal(class(xc)[1], "data.frame")
  
})
