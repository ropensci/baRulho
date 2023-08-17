test_that("using extended table and both marker", {

    data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[1:2])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  master <- degradation_est[degradation_est$sound.files == "master.wav", ]
  
  pks <- find_markers(X = master, test.files = unique(degradation_est$sound.files)[2])
  
  expect_equal(nrow(pks), 2)
  
  expect_equal(ncol(pks), 7)
  
  expect_equal(class(pks)[1], "data.frame")
  
})

test_that("using data frame and start marker", {

    data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[1:2])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  master <- as.data.frame(degradation_est[degradation_est$sound.files == "master.wav", ])
  
  pks <- find_markers(X = master, markers = "start_marker", test.files = unique(degradation_est$sound.files)[2])

  expect_equal(nrow(pks), 1)
  
  expect_equal(ncol(pks), 6)
  
  expect_true(pks$scores >  0.6)
  
  
  expect_equal(class(pks)[1], "data.frame")
  
})
