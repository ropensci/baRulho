test_that("basic", {
  # load example data
  data("test_sounds_est")
  
  
  test_sounds_est <-
    test_sounds_est[test_sounds_est$sound.id == "freq9",]
  
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
    cols = viridis::magma(4, alpha = 0.3),
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
    test_sounds_est[test_sounds_est$sound.id == "freq9",]
  
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

test_that("errors on duplicated sound.id within a transect/distance combination", {
  data("test_sounds_est")

  td <- tempdir()
  for (i in unique(test_sounds_est$sound.files))
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  options(sound.files.path = td)

  X <- as.data.frame(test_sounds_est[test_sounds_est$sound.id == "freq9",])
  X <- set_reference_sounds(X)

  # fabricate a collision that is a duplicate only at the transect/distance
  # level (not at the sound.files/distance level): relabel a row from a
  # different file/transect/distance so it collides with an existing
  # transect == "open" & distance == 1 row, while keeping a different
  # 'sound.files'. This specifically targets the transect-aware duplicate
  # check (as opposed to the sound.files-based one used for other functions).
  dup_row <- X[X$transect == "closed" & X$distance == 10,][1, ]
  dup_row$transect <- "open"
  dup_row$distance <- 1
  dup_row$selec <- max(X$selec) + 1
  Xdup <- rbind(X, dup_row)

  # the message is specific to the duplicate check being reached at all;
  # if the check were skipped (e.g. due to a broken 'fun' dispatch) the
  # function instead crashes later with an unrelated, uninformative error.
  expect_error(
    suppressWarnings(plot_degradation(X = Xdup, pb = FALSE, dest.path = tempdir())),
    regexp = "Duplicated 'sound.id' labels are not allowed within a.*transect/distance combination"
  )

})
