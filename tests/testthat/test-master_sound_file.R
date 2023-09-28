test_that("save file", {
  
  data("master_est")
  
    # create master sound file
    master.sel.tab <- master_sound_file(
      X = master_est, file.name = "example_master",
      dest.path = tempdir(), gap.duration = 0.3)

    expect_true(file.exists(file.path(tempdir(), "example_master.wav")))
      
  unlink(file.path(tempdir(), "example_master.wav"))
})


test_that("no top.freq, bottom.freq, sound id columns ", {
  
  data("master_est")
  
  # create master sound file
  master.sel.tab <- master_sound_file(
    X = master_est[,1:4], file.name = "example_master",
    dest.path = tempdir(), gap.duration = 0.3)
  
  expect_true(file.exists(file.path(tempdir(), "example_master.wav")))
  
  unlink(file.path(tempdir(), "example_master.wav"))
  
})
