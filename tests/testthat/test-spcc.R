test_that("using extended table method 1", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]
  
  X <- set_reference_sounds(X, pb = FALSE)
  
  xc <- spcc(X = X, pb = FALSE)
  
  expect_equal(sum(is.na(xc$cross.correlation)), 9)
  
  expect_equal(nrow(xc), 25)
  
  expect_equal(ncol(xc), 11)
  
  expect_equal(class(xc)[1], "extended_selection_table")
  
})

test_that("using data frame", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[-1])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <-
    as.data.frame(test_sounds_est[test_sounds_est$sound.files != "master.wav",])
  
  X <- set_reference_sounds(X, method = 2)
  
  expect_warning(xc <- spcc(X = X, pb = FALSE))
  
  expect_equal(sum(is.na(xc$cross.correlation)), 13)
  
  expect_equal(nrow(xc), 25)
  
  expect_equal(ncol(xc), 11)
  
  expect_equal(class(xc)[1], "data.frame")
  
})
