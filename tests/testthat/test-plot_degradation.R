test_that("basic", {
  # load example data
  data("degradation_est")
  
  # create subset of data with only re-recorded files
  rerecorded_est <-
    degradation_est[degradation_est$sound.files != "master.wav",]
  
  # order so spectrograms from same sound id as close in the graph
  rerecorded_est <-
    rerecorded_est[order(rerecorded_est$sound.id),]
  
  # set directory to save image files
  options(dest.path = tempdir())
  
  X <- set_reference_sounds(rerecorded_est[rerecorded_est$sound.files != "master.wav", ])
  
  
  # plot (look into temporary working directory `tempdir()`)
    # plot degradation spectrograms
    plot_degradation(
      X = X, nrow = 3, ovlp = 95,
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
  data("degradation_est")
  
  # create subset of data with only re-recorded files
  rerecorded_est <-
    degradation_est[degradation_est$sound.files != "master.wav",]
  
  # order so spectrograms from same sound id as close in the graph
  rerecorded_est <-
    rerecorded_est[order(rerecorded_est$sound.id),]
  
  # set directory to save image files
  options(dest.path = tempdir())
  
  X <- set_reference_sounds(rerecorded_est)
  
  # plot (look into temporary working directory `tempdir()`)
  # plot degradation spectrograms
  plot_degradation(
    X = X, nrow = 3000, ovlp = 95,
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