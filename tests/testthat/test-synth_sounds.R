test_that("basic", {

  ss <- synth_sounds(
    mar = 0.01,
    frequencies = c(1, 2, 3, 5),
    durations = 0.1,
    fm = TRUE,
    am = TRUE,
    nharmonics = 4,
    shuffle = TRUE,
    replicates = 3,
  )
  
  expect_equal(nrow(ss), 96)
  
  expect_equal(ncol(ss), 14)
  
  expect_equal(class(ss)[1], "extended_selection_table")
  
})

test_that("shuffle false", {
  
  ss <- synth_sounds(
    mar = 0.01,
    frequencies = c(1, 2, 3, 5),
    durations = 0.1,
    fm = TRUE,
    am = TRUE,
    nharmonics = 4,
    shuffle = FALSE,
    replicates = 3
  )
  
  expect_equal(nrow(ss), 96)
  
  expect_equal(ncol(ss), 14)
  
  expect_equal(class(ss)[1], "extended_selection_table")
  
})
