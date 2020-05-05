# *baRulho 1.0.2*

* Fix sign error in signal amplitude measurements (`signal_to_noise_ratio()` and `excess_attenuation()`)
* New function `noise_profile()` 
* rename `snr()` to `signal_to_noise_ratio()`
* New function `tail_to_signal_ratio()` to measure reverberations
* Fix bug on `excess_attenuation()` when method = 1
* Added type argument to `excess_attenuation()` to run "Darden" version of excess attenuation 

# *baRulho 1.0.1*

* New function `search_templates()` to find signals in re-recorded sound files
* New function `align_test_files()` to set time of signals in aligned re-recorded files
* Parallel available on internal `prep_X_bRlo_int()` function
* Data frame are also returned by most functions

# *baRulho 1.0.0*

* First release
