test_that("on extended table", {
  
    # load example data
    data("degradation_est")

    # create subset of data with only re-recorded files
    rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]

    # measure on custom noise reference
    np <- noise_profile(X = rerecorded_est, mar = 0.01, pb = FALSE, noise.ref = "custom")

    expect_equal(nrow(np), 50)
    
    expect_equal(ncol(np), 3)
    
    expect_equal(class(np)[1], "data.frame")
})


test_that("on extended table without ambient sound id", {
  
  # load example data
  data("degradation_est")
  
  # create subset of data with only re-recorded files
  rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]
  
  # remove noise selections so noise is measured right before the signals
  rerecorded_est <- rerecorded_est[rerecorded_est$sound.id != "ambient", ]
  
  np <- noise_profile(X = rerecorded_est, mar = 0.01, pb = FALSE, noise.ref = "adjacent")
  
  expect_equal(nrow(np), 50)
  
  expect_equal(ncol(np), 3)
  
  expect_equal(class(np)[1], "data.frame")
  
})


test_that("using data frame", {
  data("degradation_est")
  
  # set temporary directory
  td <- tempdir()  
  
  for (i in unique(degradation_est$sound.files)[2:3])
    writeWave(object = attr(degradation_est, "wave.objects")[[i]], file.path(td, i))
  
  options(sound.files.path = td, pb = FALSE)
  
  np <-  
    noise_profile(mar = 0.01, pb = FALSE, noise.ref = "adjacent", files = unique(degradation_est$sound.files)[2:3])
  
    expect_equal(nrow(np), 20)

    expect_equal(ncol(np), 3)
  
  expect_equal(class(np)[1], "data.frame")
  
})
