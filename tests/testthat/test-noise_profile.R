test_that("on extended table", {
  # load example data
  data("test_sounds_est")
  
  # measure on custom noise reference
  np <-
    noise_profile(
      X = test_sounds_est,
      mar = 0.01,
      pb = FALSE,
      noise.ref = "custom"
    )
  
  expect_equal(nrow(np), 50)
  
  expect_equal(ncol(np), 3)
  
  expect_equal(class(np)[1], "data.frame")
})


test_that("on extended table without ambient sound id", {
  # load example data
  data("test_sounds_est")
  
  # remove noise selections so noise is measured right before the signals
  test_sounds_est <-
    test_sounds_est[test_sounds_est$sound.id != "ambient",]
  
  np <-
    noise_profile(
      X = test_sounds_est,
      mar = 0.01,
      pb = FALSE,
      noise.ref = "adjacent"
    )
  
  expect_equal(nrow(np), 50)
  
  expect_equal(ncol(np), 3)
  
  expect_equal(class(np)[1], "data.frame")
  
})


test_that("using data frame", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[2:3])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  np <-
    noise_profile(
      mar = 0.01,
      pb = FALSE,
      noise.ref = "adjacent",
      files = unique(test_sounds_est$sound.files)[2:3]
    )
  
  expect_equal(nrow(np), 20)
  
  expect_equal(ncol(np), 3)
  
  expect_equal(class(np)[1], "data.frame")
  
})


test_that("using entire files", {
  data("test_sounds_est")
  
  # set temporary directory
  td <- tempdir()
  
  for (i in unique(test_sounds_est$sound.files)[2:3])
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  np <-
    noise_profile(mar = 0.01, pb = FALSE)
  
  expect_equal(nrow(np), 60)
  
  expect_equal(ncol(np), 3)
  
  expect_equal(class(np)[1], "data.frame")
  
})