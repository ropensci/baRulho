test_that("basic", {


    # load example data
    data("test_sounds_est")

    
    

    # create "unaligned_test_sounds_est" by adding noise to "test_sounds_est" start and end
    unaligned_test_sounds_est <- test_sounds_est
    set.seed(123)
    noise_time <- sample(c(0.005, -0.005, 0.006, -0.006, 0, 0.002, -0.002),
      nrow(unaligned_test_sounds_est),
      replace = TRUE
    )
    attr(unaligned_test_sounds_est, "check.res")$start <-
      unaligned_test_sounds_est$start <- unaligned_test_sounds_est$start + noise_time
    attr(unaligned_test_sounds_est, "check.res")$end <- unaligned_test_sounds_est$end <-
      unaligned_test_sounds_est$end + noise_time

    # re align
    rts <- realign_test_sounds(X = unaligned_test_sounds_est)
  
    
    expect_equal(nrow(rts), 25)
    
    expect_equal(ncol(rts), 9)
    
    expect_true(mean(unaligned_test_sounds_est$start) >  mean(rts$start))
    
    
    expect_equal(class(rts)[1], "extended_selection_table")
    
  })
