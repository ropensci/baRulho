test_that("basic", {
    # load example data
    data("degradation_est")

    # create subset of data with only re-recorded files
    rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]

    # add reference to X
    X <- set_reference_sounds(X = rerecorded_est[rerecorded_est$sound.id == rerecorded_est$sound.id[2], ])

    # create plots
    plot_blur_ratio(X = X, dest.path = tempdir(), pb = FALSE)
    
    fls <-
      list.files(path = tempdir(),
                 pattern = "^blur_ratio",
                 full.names = TRUE)
    
    expect_length(fls, 4)
    
    unlink(fls)
    })