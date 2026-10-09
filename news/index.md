# Changelog

## baRulho 2.2.0

#### MINOR IMPROVEMENTS

- Fix bug in
  [`spcc()`](https://marce10.github.io/baRulho/reference/spcc.md) in
  which the `cor.method` argument was ignored (always computed Pearson
  correlation)
- Fix bug in
  [`tail_to_signal_ratio()`](https://marce10.github.io/baRulho/reference/tail_to_signal_ratio.md)
  that used the wrong dB scale (natural log instead of `20 * log10()`,
  inconsistent with
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
  and
  [`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md))
- Fix bug in
  [`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md)’s
  internal sliding-window search that never actually slid across
  candidate alignments
- Fix bug in
  [`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md)
  and
  [`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md)
  in which an internal dispatch issue caused some of their own argument
  validation checks to be silently skipped
- Fix bug in
  [`synth_sounds()`](https://marce10.github.io/baRulho/reference/synth_sounds.md)
  that crashed when called with a single frequency and a single duration
- Fix bug in
  [`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md)
  in which the `wl` argument was not reaching the internal
  peak-frequency calculation
- Fix misspelled output column `atmopheric.attenuation` in
  [`attenuation()`](https://marce10.github.io/baRulho/reference/attenuation.md)
  (now `atmospheric.attenuation`)
- Fix several smaller documentation inconsistencies (incorrect argument
  defaults, mislocated columns, stale formulas)
- Reformat all exported function documentation using markdown roxygen
  for a more consistent and readable reference site
- Vignettes now pull their citation directly from `inst/CITATION`
  instead of hardcoded text, and no longer describe a stale version of
  the example data set

## baRulho 2.1.7

CRAN release: 2026-07-24

#### MINOR IMPROVEMENTS

- Fix bug in
  [`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md)
  that ask for a ‘scores’ column in Y

- Improve argument validation

## baRulho 2.1.6

CRAN release: 2025-10-24

#### MINOR IMPROVEMENTS

- Improve documentation of several functions

## baRulho 2.1.5

CRAN release: 2025-07-17

### NEW FEATURES

- new argument `freq_range` in
  [`synth_sounds()`](https://marce10.github.io/baRulho/reference/synth_sounds.md)
  to limit the frequency range of modulated sounds
- new function
  [`spot_ambient_noise()`](https://marce10.github.io/baRulho/reference/spot_ambient_noise.md)
  to find segments of ambient noise to be used as reference by other
  functions
- new argument `seed` in
  [`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md)
  to let users set a seed for reproducibility (or not)

#### MINOR IMPROVEMENTS

- Improve error messages for functions working with master annotations
- Remove requirement of “sound.id” column in
  [`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md)
  and
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
- Fix bug in
  [`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md)
  in which ‘…’ were not sent internally to
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
- Improved error message when missing columns in supplied `Y` data
  frames

## baRulho 2.1.4

#### MINOR IMPROVEMENTS

- Fix bug in
  [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
  and
  [`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md)
  when missing one marker that is not used

## baRulho 2.1.3

CRAN release: 2025-01-24

#### MINOR IMPROVEMENTS

- Improved performance issue of
  [`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md)

## baRulho 2.1.2

CRAN release: 2024-08-31

#### MINOR IMPROVEMENTS

- Improved tracking of progress messages when more than 1 progress bar
  is used by a process

## baRulho 2.1.1

- Update requested by CRAN to fix package reference issue. No changes in
  the code.

## baRulho 2.1.0

CRAN release: 2024-04-21

### NEW FEATURES

- new function
  [`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md)
  to modifying signal-to-noise ratio by adding synthetic noise
- new function
  [`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md)
  that generates an interactive plot for manual adjustment of time
  alignments
- new function
  [`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md)
  that generates plots previously created by
  [`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md)
  and
  [`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md)
- new function
  [`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md)
  to visually compare sounds across distances
- new function
  [`plot_aligned_sounds()`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md)
  to visually assess precision of test sound alignment
- [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md)
  now can take several markers as input and select that with the highest
  correlation score for aligning test files
- New function
  [`attenuation()`](https://marce10.github.io/baRulho/reference/attenuation.md)
- Data frames and selection tables can be used as input data
- Added new methods to
  [`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md)
  and `excess_atenuation()`
- `spectral_correlation()` and `spectral_blur_ratio()` renamed to
  [`spectrum_correlation()`](https://marce10.github.io/baRulho/reference/spectrum_correlation.md),
  and
  [`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md)
  respectively

#### MINOR IMPROVEMENTS

- optimize performance in
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
- `atmospheric_attenuation()` is no longer available (use
  [`attenuation()`](https://marce10.github.io/baRulho/reference/attenuation.md)
  instead)
- Improved documentation for all functions
- Fix bug in
  [`spcc()`](https://marce10.github.io/baRulho/reference/spcc.md) and
  [`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md)
- ‘markers’ argument deprecated in
  [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md)
- fix bug in
  [`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md)
- [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
  compares the time difference between markers to that in the master
  sound file as a measure of precision
- [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
  can run several markers as templates on the same files
- ‘template.rows’ argument deprecated in
  [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
- `spcc_align()` renamed
  [`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md)
- `search_templates()` renamed
  [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
- ‘output’ argument deprecated
- ‘parallel’ argument deprecated and replaced by ‘cores’

## baRulho 1.0.6

CRAN release: 2022-03-01

- Update requested by CRAN

## baRulho 1.0.5

CRAN release: 2021-04-21

- Update requested by CRAN

## baRulho 1.0.4

CRAN release: 2021-03-09

#### MINOR IMPROVEMENTS

- [`warbleR::freq_range_detec()`](https://marce10.github.io/warbleR/reference/freq_range_detec.html)
  is now used internally to detect frequency range of markers in
  [`master_sound_file()`](https://marce10.github.io/baRulho/reference/master_sound_file.md)

## baRulho 1.0.3

CRAN release: 2021-02-11

- Update requested by CRAN

## baRulho 1.0.3

CRAN release: 2021-02-11

#### MINOR IMPROVEMENTS

- New argument ‘marker’ in
  [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md)
  to control if the start or end marker is being used for aligning
- Fix bug when detecting several templates per sound file in
  `search_templates()`

## baRulho 1.0.2

CRAN release: 2020-06-07

#### MINOR IMPROVEMENTS

- Windows length is converted to even in all functions
- Fix sign error in signal amplitude measurements
  ([`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
  and
  [`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md))
- New function
  [`noise_profile()`](https://marce10.github.io/baRulho/reference/noise_profile.md)
- rename `snr()` to
  [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
- New function
  [`tail_to_signal_ratio()`](https://marce10.github.io/baRulho/reference/tail_to_signal_ratio.md)
  to measure reverberations
- Fix bug on
  [`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md)
  when method = 1
- Added type argument to
  [`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md)
  to run “Darden” version of excess attenuation

## baRulho 1.0.1

CRAN release: 2020-03-09

#### NEW FEATURES

- New function `search_templates()` to find signals in re-recorded sound
  files
- New function
  [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md)
  to set time of signals in aligned re-recorded files
- Parallel available on internal `prep_X_bRlo_int()` function
- Data frame are also returned by most functions

## baRulho 1.0.0

CRAN release: 2020-02-22

- First release
