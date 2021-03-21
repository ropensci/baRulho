# *baRulho 1.0.4*

* `warbleR::freq_range_detec()` is now used internally to detect frequency range of markers in `master_sound_file()` 

# *baRulho 1.0.3*

* Update requested by CRAN

# *baRulho 1.0.3*

* New argument 'marker' in `align_test_files()` to control if the start or end marker is being used for aliging
* Fix bug when detecting several templates per sound file in `search_templates()`

# *baRulho 1.0.2*

* Windows length is converted to even in all functions
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
