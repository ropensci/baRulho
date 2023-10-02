test_that("using extended table", {
  data("test_sounds_est")
  
  X <-
    test_sounds_est[test_sounds_est$sound.id == "freq:9",]
  
  
  X <- set_reference_sounds(X, method = 2)
  
  dd <-
    detection_distance(X = X[X$distance %in% c(1, 10),],
                       spl.cutoff = 5,
                       mar = 0.05)
  
  expect_equal(sum(is.na(dd$detection.distance)), 2)
  
  expect_equal(nrow(dd), 3)
  
  expect_equal(ncol(dd), 11)
  
  expect_equal(class(dd)[1], "extended_selection_table")
  
})

test_that("using data frame", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[-1])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  X <-
    as.data.frame(test_sounds_est[test_sounds_est$sound.id == "freq:9",])
  
  X <- set_reference_sounds(X, pb = FALSE)
  
  dd <-
    detection_distance(X = X[X$distance %in% c(1, 10),],
                       spl.cutoff = 5,
                       mar = 0.05)
  
  expect_equal(sum(is.na(dd$detection.distance)), 1)
  
  expect_equal(nrow(dd), 3)
  
  expect_equal(ncol(dd), 11)
  
  expect_equal(class(dd)[1], "data.frame")
  
})
