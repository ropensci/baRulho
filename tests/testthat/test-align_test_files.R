test_that("using extended table and both marker", {
  data("test_sounds_est")
  data("master_est")
  # set temporary directory
  td <- tempdir()
  
  unlink(list.files(
    path = td,
    pattern = ".wav",
    ignore.case = TRUE,
    full.names = TRUE
  ))
  
  for (i in unique(test_sounds_est$sound.files)[1:2])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  options(sound.files.path = td, pb = T)
  
  # save master file
  writeWave(object = attr(master_est, "wave.objects")[[1]],
            file.path(td, "master.wav"))
  
  options(sound.files.path = td, pb = T)
  
  pks <-
    find_markers(X = master_est,
                 test.files = unique(test_sounds_est$sound.files)[2])
  
  alg.tests <- align_test_files(X = master_est, Y = pks, pb = FALSE)
  
  expect_true(is_extended_selection_table(alg.tests))
  
  expect_equal(nrow(alg.tests), 7)
  
  expect_equal(ncol(alg.tests), 8)
  
  expect_equal(class(alg.tests)[1], "extended_selection_table")
  
})

test_that("using data frame and start marker", {
  data("test_sounds_est")
  data("master_est")
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[1:2])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  options(sound.files.path = td, pb = T)
  
  # save master file
  writeWave(object = attr(master_est, "wave.objects")[[1]],
            file.path(td, "master.wav"))
  
  options(sound.files.path = td, pb = FALSE)
  
  master <-
    as.data.frame(test_sounds_est[test_sounds_est$sound.files == "master.wav",])
  
  pks <-
    find_markers(
      X = master_est,
      markers = "start_marker",
      test.files = unique(test_sounds_est$sound.files)[2]
    )
  
  alg.tests <- align_test_files(X = master_est, Y = pks)
  
  expect_true(is_extended_selection_table(alg.tests))
  
  expect_equal(ncol(alg.tests), 8)
  
  expect_equal(nrow(alg.tests), 7)
  
})
