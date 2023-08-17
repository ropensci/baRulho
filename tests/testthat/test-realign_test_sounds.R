test_that("basic", {


    # load example data
    data("degradation_est")

    # create subset of data with only re-recorded files
    rerecorded_est <- degradation_est[degradation_est$sound.files != "master.wav", ]

    # create "unaligned_rerecorded_est" by adding noise to "rerecorded_est" start and end
    unaligned_rerecorded_est <- rerecorded_est
    set.seed(123)
    noise_time <- sample(c(0.005, -0.005, 0.006, -0.006, 0, 0.002, -0.002),
      nrow(unaligned_rerecorded_est),
      replace = TRUE
    )
    attr(unaligned_rerecorded_est, "check.res")$start <-
      unaligned_rerecorded_est$start <- unaligned_rerecorded_est$start + noise_time
    attr(unaligned_rerecorded_est, "check.res")$end <- unaligned_rerecorded_est$end <-
      unaligned_rerecorded_est$end + noise_time

    # re align
    rts <- realign_test_sounds(X = unaligned_rerecorded_est)
  
    
    expect_equal(nrow(rts), 25)
    
    expect_equal(ncol(rts), 9)
    
    expect_true(mean(unaligned_rerecorded_est$start) >  mean(rts$start))
    
    
    expect_equal(class(rts)[1], "extended_selection_table")
    
  })
