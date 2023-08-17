test_that("using extended table", {
  data("degradation_est")
  
  X <- degradation_est[degradation_est$sound.files != "master.wav", ]
  
  dd <- detection_distance(X = X[X$distance == 1, ], spl.cutoff = 5, mar = 0.05)
  
  expect_equal(sum(is.na(dd$detection.distance)), 1)
  
  expect_equal(nrow(dd), 5)
  
  expect_equal(ncol(dd), 11)
  
  expect_equal(class(dd)[1], "extended_selection_table")
  
})

test_that("using data frame", {
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[-1])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <- as.data.frame(degradation_est[degradation_est$sound.files != "master.wav", ])
  
  dd <- detection_distance(X = X[X$distance == 1, ], spl.cutoff = 5, mar = 0.05)
  
  expect_equal(sum(is.na(dd$detection.distance)), 1)
  
  expect_equal(nrow(dd), 5)
  
  expect_equal(ncol(dd), 11)
  
  expect_equal(class(dd)[1], "data.frame")
  
})
