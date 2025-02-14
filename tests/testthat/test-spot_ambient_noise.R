test_that("basic", {
  td <- tempdir()  
  data("test_sounds_est")
  
  
  # save example files in working director to recreate a case in which working
  # with sound files instead of extended selection tables.
  # This doesn't have to be done with your own data as you will
  # have them as sound files already.
  for (i in unique(test_sounds_est$sound.files)[1:2]) {
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]],
              file.path(tempdir(), i))
  }
  
  # plot_aligned_sounds(X = alg.tests, path = td, dest.path = td)
  
  test_sounds_df <- as.data.frame(test_sounds_est)
  test_sounds_df <- test_sounds_df[test_sounds_df$sound.id != "ambient", ]
  
  test_sounds_df <- test_sounds_df[test_sounds_df$sound.files %in% unique(test_sounds_est$sound.files)[1:2], ]
  
  
  # closest to mean
  san <- spot_ambient_noise(X = test_sounds_df, path = td, length = 0.12, ovlp = 20)

  expect_true(any(san$sound.id == "ambient"))
  
  expect_equal(nrow(san), 10)
  
})
