test_that("basic", {
    # load example data
    data("test_sounds_est")

    
    

    # add reference to X
    X <- set_reference_sounds(X = test_sounds_est[test_sounds_est$sound.id == test_sounds_est$sound.id[2], ])

    # create plots
    plot_blur_ratio(X = X, dest.path = tempdir(), pb = FALSE)
    
    fls <-
      list.files(path = tempdir(),
                 pattern = "^blur_ratio",
                 full.names = TRUE)
    
    expect_length(fls, 4)
    
    unlink(fls)
    })
