test_that("basic", {
 
  # measure attenuation
    att <- attenuation(f = 2000, dist = 50, dist0 = 1)
  
    expect_equal(nrow(att), 1)
    
    expect_equal(ncol(att), 6)

    expect_equal(att$habitat.attenuation, 1.96)
    
    expect_false(anyNA(att))
    
})
