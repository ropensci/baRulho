test_that("using extended table and method 1", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]
  
  snr <- tail_to_signal_ratio(X = X, mar = 0.1)
  
  expect_equal(sum(is.na(snr$tail.to.signal.ratio)), 5)
  
  expect_equal(nrow(snr), 25)
  
  expect_equal(ncol(snr), 10)
  
  expect_equal(class(snr)[1], "extended_selection_table")
  
})

test_that("using data frame", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files))
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <-
    as.data.frame(test_sounds_est[test_sounds_est$sound.files != "master.wav",])
  
  snr <- tail_to_signal_ratio(X = X, mar = 0.1)
  
  expect_equal(sum(is.na(snr$tail.to.signal.ratio)), 5)
  
  expect_equal(nrow(snr), 25)
  
  expect_equal(ncol(snr), 10)
  
  expect_equal(class(snr)[1], "data.frame")

})

test_that("uses the 20 * log10() dB convention, not a different scale", {
  data("test_sounds_est")

  X <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav",]

  snr <- tail_to_signal_ratio(X = X, mar = 0.1, pb = FALSE)

  tsr_vals <- snr$tail.to.signal.ratio[!is.na(snr$tail.to.signal.ratio)]

  # golden value: pins the result to the 20 * log10(tail_RMS / sig_RMS) formula
  # used consistently by signal_to_noise_ratio() and excess_attenuation().
  # A regression to 10 * log() (natural log) would shift this by a factor
  # of ~1.151 and fail this check.
  expect_equal(round(tsr_vals[1:3], 3), c(-15.559, -6.658, -13.342))

})
