test_that("sav file", {
  
  data("degradation_est")
  

    # create master sound file
    master.sel.tab <- master_sound_file(
      X = degradation_est[degradation_est$sound.files == "master.wav",], file.name = "example_master",
      dest.path = tempdir(), gap.duration = 0.3)

    expect_true(file.exists(file.path(tempdir(), "example_master.wav")))
      
  
})
