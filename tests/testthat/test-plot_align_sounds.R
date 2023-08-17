test_that("basic", {
    # load example data
    data("degradation_est")

    # create subset of data with only re-recorded files
    rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]

  # plot (look into temporary working directory `tempdir()`)
  plot_align_sounds(X = rerecorded_est[rerecorded_est$sound.files == rerecorded_est$sound.files[1]], dest.path = tempdir(), duration = 1, ovlp = 0)

  fls <- list.files(path = tempdir(), pattern = "^plot_align", full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
  })
