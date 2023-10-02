test_that("basic", {
  # load example data
  data("test_sounds_est")
  
  
  test_sounds_est <-
    test_sounds_est[test_sounds_est$sound.id == "freq:9",]
  
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
    ovlp = 0,
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


test_that("many ros", {
  # load example data
  data("test_sounds_est")
  
  
  test_sounds_est <-
    test_sounds_est[test_sounds_est$sound.id == "freq:9",]
  
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
    ovlp = 0,
    cols = viridis::magma(4, alpha = 0.3),
    palette = viridis::magma,
    env.smooth = 200,
    hop.size = 11.5,
    collevels = seq(-120, 0, 5),
    flim = c("-1", "+1"),
    envelope = TRUE,
    spectrum = TRUE,
    heights = c(4, 1),
    widths = c(5, 1),
    margins = c(2, 1),
    row.height = 2,
    col.width = 2,
    res = 120,
  )
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^plot_degrad",
               full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
})
