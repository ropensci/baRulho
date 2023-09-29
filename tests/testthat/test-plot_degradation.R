test_that("basic", {
  # load example data
  data("test_sounds_est")
  
  
  test_sounds_est <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav", ]
  
  # order so spectrograms from same sound id as close in the graph
  test_sounds_est <-
    test_sounds_est[order(test_sounds_est$sound.id), ]
  
  # set directory to save image files
  options(dest.path = tempdir())
  
  X <-
    set_reference_sounds(test_sounds_est[test_sounds_est$sound.files != "master.wav",])
  
  
  # plot (look into temporary working directory `tempdir()`)
  # plot degradation spectrograms
  plot_degradation(
    X = X,
    nrow = 3,
    ovlp = 95,
    colors = viridis::magma(4, alpha = 0.3),
    palette = viridis::magma
  )
  
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^plot_degrad",
               full.names = TRUE)
  
  expect_length(fls, 3)
  
  unlink(fls)
})


test_that("many ros", {
  # load example data
  data("test_sounds_est")
  
  
  test_sounds_est <-
    test_sounds_est[test_sounds_est$sound.files != "master.wav", ]
  
  # order so spectrograms from same sound id as close in the graph
  test_sounds_est <-
    test_sounds_est[order(test_sounds_est$sound.id), ]
  
  # set directory to save image files
  options(dest.path = tempdir())
  
  X <- set_reference_sounds(test_sounds_est)
  
  # plot (look into temporary working directory `tempdir()`)
  # plot degradation spectrograms
  plot_degradation(
    X = X,
    nrow = 3000,
    ovlp = 95,
    colors = viridis::magma(4, alpha = 0.3),
    palette = viridis::magma
  )
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^plot_degrad",
               full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
})
