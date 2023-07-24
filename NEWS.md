# *baRulho 1.1.0*

* `align_test_files()` now can take several markers as input and select that with the highest correlation score for aligining test files
* 'markers' argument deprecated in `align_test_files()`
* fix bug in `spcc_align()`
* `find_markers()` compares the time difference between markers to that in the master sound file as a measure of precision
* `find_markers()` can run several markers as templates on the same files
* 'template.rows' argument deprecated in `find_markers()`
* `search_templates()` renamed `find_markers()`
* 'parallel' argument deprecated and replaced by 'cores'
* New functions `attenuation()` and 
* Data frames and selection tables can be used as input data
* Added new methods to `blur_ratio()` and `excess_atenuation()`
* `spectral_correlation()` and `spectral_blur_ratio()` renamed to `spectrum_correlation()`, and `spectrum_blur_ratio()` respectively
* Fix bug in `spcc()` and `excess_attenuation()`

# *baRulho 1.0.6*

* Update requested by CRAN

# *baRulho 1.0.5*

* Update requested by CRAN

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
