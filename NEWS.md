baRulho 2.1.2
=========================

### MINOR IMPROVEMENTS

* Improved tracking of progress messages when more than 1 progress bar is used by a process


baRulho 2.1.1
=========================

* Update requested by CRAN to fix package reference issue. No changes in the code.


baRulho 2.1.0
=========================

## NEW FEATURES

* new function `add_noise()` to modifying signal-to-noise ratio by adding synthetic noise
* new function `manual_realign()` that generates an interactive plot for manual adjustment of time alignments
* new function `plot_blur_ratio()` that generates plots previously created by `blur_ratio()` and `spectrum_blur_ratio()`
* new function `plot_degradation()` to visually compare sounds across distances
* new function `plot_aligned_sounds()` to visually assess precision of test sound alignment
* `align_test_files()` now can take several markers as input and select that with the highest correlation score for aligning test files
* New function `attenuation()`
* Data frames and selection tables can be used as input data
* Added new methods to `blur_ratio()` and `excess_atenuation()`
* `spectral_correlation()` and `spectral_blur_ratio()` renamed to `spectrum_correlation()`, and `spectrum_blur_ratio()` respectively

### MINOR IMPROVEMENTS

* optimize performance in `signal_to_noise_ratio()`
* `atmospheric_attenuation()` is no longer available (use `attenuation()` instead)
* Improved documentation for all functions
* Fix bug in `spcc()` and `excess_attenuation()`
* 'markers' argument deprecated in `align_test_files()`
* fix bug in `auto_realign()`
* `find_markers()` compares the time difference between markers to that in the master sound file as a measure of precision
* `find_markers()` can run several markers as templates on the same files
* 'template.rows' argument deprecated in `find_markers()`
* `spcc_align()` renamed `auto_realign()`
* `search_templates()` renamed `find_markers()`
* 'output' argument deprecated
* 'parallel' argument deprecated and replaced by 'cores'

baRulho 1.0.6
=========================

* Update requested by CRAN

baRulho 1.0.5
=========================

* Update requested by CRAN

baRulho 1.0.4
=========================

### MINOR IMPROVEMENTS

* `warbleR::freq_range_detec()` is now used internally to detect frequency range of markers in `master_sound_file()` 

baRulho 1.0.3
=========================

* Update requested by CRAN

baRulho 1.0.3
=========================

### MINOR IMPROVEMENTS

* New argument 'marker' in `align_test_files()` to control if the start or end marker is being used for aligning
* Fix bug when detecting several templates per sound file in `search_templates()`

baRulho 1.0.2
=========================

### MINOR IMPROVEMENTS

* Windows length is converted to even in all functions
* Fix sign error in signal amplitude measurements (`signal_to_noise_ratio()` and `excess_attenuation()`)
* New function `noise_profile()` 
* rename `snr()` to `signal_to_noise_ratio()`
* New function `tail_to_signal_ratio()` to measure reverberations
* Fix bug on `excess_attenuation()` when method = 1
* Added type argument to `excess_attenuation()` to run "Darden" version of excess attenuation 

baRulho 1.0.1
=========================

### NEW FEATURES

* New function `search_templates()` to find signals in re-recorded sound files
* New function `align_test_files()` to set time of signals in aligned re-recorded files
* Parallel available on internal `prep_X_bRlo_int()` function
* Data frame are also returned by most functions

baRulho 1.0.0
=========================

* First release
