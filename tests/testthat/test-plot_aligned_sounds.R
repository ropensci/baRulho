test_that("basic", {
    # load example data
    data("test_sounds_est")

    
    

  # plot (look into temporary working directory `tempdir()`)
  plot_aligned_sounds(X = test_sounds_est[test_sounds_est$sound.files == test_sounds_est$sound.files[1]], dest.path = tempdir(), duration = 1, ovlp = 0)

  fls <- list.files(path = tempdir(), pattern = "^plot_align", full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
  })
