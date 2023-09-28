test_that("using extended table and method 1", {
  data("test_sounds_est")
  
  X <- test_sounds_est[test_sounds_est$sound.files != "master.wav", ]
  
  snr <- signal_to_noise_ratio(X = X, mar = 0.1)
  
  expect_equal(sum(is.na(snr$signal.to.noise.ratio)), 5)
  
  expect_equal(nrow(snr), 25)
  
  expect_equal(ncol(snr), 10)
  
  expect_equal(class(snr)[1], "extended_selection_table")
  
})

test_that("using data frame", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(test_sounds_est$sound.files)[-1])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <- as.data.frame(test_sounds_est[test_sounds_est$sound.files != "master.wav", ])
  
  snr <- signal_to_noise_ratio(X = X, mar = 0.1)
  
  expect_equal(sum(is.na(snr$signal.to.noise.ratio)), 5)
  
  expect_equal(nrow(snr), 25)
  
  expect_equal(ncol(snr), 10)
  
  expect_equal(class(snr)[1], "data.frame")
  
})

test_that("using extended table, method 1 and custom noise", {
  data("test_sounds_est")
  
  X <- test_sounds_est[test_sounds_est$sound.files != "master.wav", ]
  
  snr <- signal_to_noise_ratio(X = X, mar = 0.1, noise.ref = "custom")
  
  expect_equal(sum(is.na(snr$signal.to.noise.ratio)), 5)
  
  expect_equal(nrow(snr), 25)
  
  expect_equal(ncol(snr), 10)
  
  expect_equal(class(snr)[1], "extended_selection_table")
  
})