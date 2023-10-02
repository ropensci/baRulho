test_that("basic", {
  # load example data
  data("test_sounds_est")
  
  unlink(list.files(tempdir(), pattern = "jpeg$"))
  
  # plot (look into temporary working directory `tempdir()`)
  plot_aligned_sounds(
    X = test_sounds_est[test_sounds_est$sound.files == test_sounds_est$sound.files[1]],
    dest.path = tempdir(),
    duration = 1,
    ovlp = 0,
    hop.size = 11.6,
    wl = 512,
    path = ".",
    cores = 1,
    pb = FALSE,
    collevels = seq(-120, 0, 5),
    palette = viridis::viridis,
    mar = 0.2,
    flim = c(1, 10),
    col = "white",
    width = 7,
    height = 4,
    res = 100,
    label = TRUE,
    fast.spec = FALSE,
    srt = 0,
    cex = 1
  )
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^plot_align",
               full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
})


test_that("srt = 20", {
  # load example data
  data("test_sounds_est")
  
  unlink(list.files(tempdir(), pattern = "jpeg$"))
  
  # plot (look into temporary working directory `tempdir()`)
  plot_aligned_sounds(
    X = test_sounds_est[test_sounds_est$sound.files == test_sounds_est$sound.files[1]],
    dest.path = tempdir(),
    duration = 1,
    ovlp = 0,
    srt = 20
  )
  
  fls <-
    list.files(path = tempdir(),
               pattern = "^plot_align",
               full.names = TRUE)
  
  expect_length(fls, 1)
  
  unlink(fls)
})
