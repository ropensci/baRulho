test_that("basic", {
  
   skip_if_not(interactive())  
  
  # load example data
    data(list = c("master_est", "test_sounds_est"))

    # save example files in working director to recreate a case in which working
    # with sound files instead of extended selection tables.
    # This doesn't have to be done with your own data as you will
    # have them as sound files already.
    for (i in unique(test_sounds_est$sound.files)[1:2]) {
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(tempdir(), i))
    }

    # save master file
    writeWave(object = attr(master_est, "wave.objects")[[1]], file.path(tempdir(), "master.wav"))

    # get marker position
    markers <- find_markers(X = master_est, test.files = unique(test_sounds_est$sound.files)[2],
    path = tempdir())

    # align all test sounds
    alg.tests <- align_test_files(X = master_est, Y = markers, path = tempdir())

    # add error to alignment
    lag <- (as.numeric(as.factor(alg.tests$sound.files)) - 2) / 30
    alg.tests$start <- alg.tests$start + lag
    alg.tests$end <- alg.tests$end + lag

    realigned_est <- manual_realign(
      X = alg.tests,
      Y = master_est,
      duration = 2,
      ovlp = 50,
      hop.size = 14,
      collevels = seq(-140, 0, 5),
      palette = viridis::mako,
      ext.window = T,
      marker = "freq:9"
    )
  
  
})
