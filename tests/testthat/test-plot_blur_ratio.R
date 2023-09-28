test_that("envelope blur ratio", {
    # load example data
    data("test_sounds_est")

    # add reference to X
    X <- set_reference_sounds(X = test_sounds_est[test_sounds_est$sound.id == test_sounds_est$sound.id[2], ])

    unlink(list.files(path = tempdir(),
                         pattern = "^blur_ratio",
                         full.names = TRUE))
    
    # create plots
    plot_blur_ratio(X = X[2:3, ], dest.path = tempdir(), pb = FALSE)
    
    fls <-
      list.files(path = tempdir(),
                 pattern = "^blur_ratio",
                 full.names = TRUE)
    
    expect_length(fls, 1)
    
    unlink(fls)
    })


test_that("spectrum blur ratio", {
  # load example data
  data("test_sounds_est")
  
  # add reference to X
  X <- set_reference_sounds(X = test_sounds_est[test_sounds_est$sound.id == test_sounds_est$sound.id[2], ])
  
  unlink(list.files(path = tempdir(),
                    pattern = "^blur_ratio",
                    full.names = TRUE))
  
  # create plots
  plot_blur_ratio(X = X[2:3, ], dest.path = tempdir(), pb = FALSE, type = "spectrum")
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^spectrum_blur_ratio",
               full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
})
