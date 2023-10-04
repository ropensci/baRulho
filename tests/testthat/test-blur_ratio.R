test_that("using extended table and method 1", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]
  
  X <- set_reference_sounds(X, pb = FALSE)
  
  br <- blur_ratio(X = X, pb = FALSE)
  
  expect_equal(sum(is.na(br$blur.ratio)), 9)
  
  expect_equal(nrow(br), 25)
  
  expect_equal(ncol(br), 11)
  
  expect_equal(class(br)[1], "extended_selection_table")
  
})

test_that("using extended table and method 2 returning envelopes", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]
  
  X <- set_reference_sounds(X, method = 2)
  
  br <- blur_ratio(X = X,
                   envelopes = TRUE,
                   env.smooth = 140, pb = FALSE)
  
  expect_equal(sum(is.na(br$blur.ratio)), 13)
  
  expect_equal(nrow(br), 25)
  
  expect_equal(ncol(br), 11)
  
  expect_equal(class(br)[1], "extended_selection_table")
  
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
  
  expect_warning(br <- blur_ratio(X = X))
  
  expect_equal(sum(is.na(br$blur.ratio)), 13)
  
  expect_equal(nrow(br), 25)
  
  expect_equal(ncol(br), 11)
  
  expect_equal(class(br)[1], "data.frame")
  
})
