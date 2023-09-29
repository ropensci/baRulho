test_that("using extended table and both marker", {
  data("test_sounds_est")
  data("master_est")
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[1:2])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  options(sound.files.path = td, pb = T)
  
  # save master file
  writeWave(object = attr(test_sounds_est, "wave.objects")[[1]],
            file.path(td, "master.wav"))
  
  options(sound.files.path = td, pb = T)
  
  pks <-
    find_markers(X = master_est,
                 test.files = unique(test_sounds_est$sound.files)[2])
  
  expect_equal(nrow(pks), 2)
  
  expect_equal(ncol(pks), 7)
  
  expect_equal(class(pks)[1], "data.frame")
  
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
  
  options(sound.files.path = td, pb = T)
  
  
  pks <-
    find_markers(
      X = master_est,
      markers = "start_marker",
      test.files = unique(test_sounds_est$sound.files)[2]
    )
  
  expect_equal(nrow(pks), 1)
  
  expect_equal(ncol(pks), 6)
  
  expect_true(pks$scores >  0.6)
  
  
  expect_equal(class(pks)[1], "data.frame")
  
})
