test_that("using extended table and both marker", {
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in 1:2)
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, unique(degradation_est$sound.files)[i]))
  
  options(sound.files.path = td, pb = FALSE)
  
  master <- degradation_est[degradation_est$sound.files == "master.wav", ]
  
  pks <- find_markers(X = master)
  
  alg.tests <- align_test_files(X = master, Y = pks, pb = FALSE)
  
  expect_equal(nrow(alg.tests), 5)
  
  expect_equal(ncol(alg.tests), 9)
  
  expect_equal(class(alg.tests)[1], "extended_selection_table")
  
})

test_that("using data frame and start marker", {
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in 1:2)
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, unique(degradation_est$sound.files)[i]))
  
  options(sound.files.path = td, pb = FALSE)
  
  master <- as.data.frame(degradation_est[degradation_est$sound.files == "master.wav", ])
  
  pks <- find_markers(X = master, markers = "start_marker")
  
  alg.tests <- align_test_files(X = master, Y = pks, remove.markers = FALSE)
  
  expect_equal(ncol(alg.tests), 8)
  
  expect_equal(class(alg.tests)[1], "data.frame")
  
  expect_equal(nrow(alg.tests), 7)
  
})
