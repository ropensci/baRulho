# Package index

## Prepare acoustic data

Functions that format acoustic data for playback experiments

- [`master_sound_file()`](https://marce10.github.io/baRulho/reference/master_sound_file.md)
  : Create a master sound file
- [`synth_sounds()`](https://marce10.github.io/baRulho/reference/synth_sounds.md)
  : Create synthetic sounds

## Aligning test sound files

Aligning test sounds to determine their time position

- [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md)
  : Align test sound files
- [`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md)
  : Fix small misalignments in the time position test sounds
- [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
  : Find acoustic markers on test sound files
- [`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md)
  : Plot spectrograms to check test sound files alignment
- [`plot_aligned_sounds()`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md)
  : Plot spectrograms to check test sound files alignment

## Quantify degradation

Functions for quantifying and visually exploring degradation

- [`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md)
  : Measure blur ratio in the time domain
- [`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md)
  : Measure detection distance of sound
- [`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md)
  : Measure amplitude envelope correlation
- [`excess_attenuation()`](https://marce10.github.io/baRulho/reference/excess_attenuation.md)
  : Measure excess attenuation
- [`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md)
  : Plot blur ratio
- [`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md)
  : Save multipanel plots with reference and test sounds
- [`set_reference_sounds()`](https://marce10.github.io/baRulho/reference/set_reference_sounds.md)
  : Set reference for test sounds
- [`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
  : Measure attenuation as signal-to-noise ratio
- [`spcc()`](https://marce10.github.io/baRulho/reference/spcc.md) :
  Measure spectrographic cross-correlation as a measure of sound
  distortion
- [`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md)
  : Measure blur ratio in the frequency domain
- [`spectrum_correlation()`](https://marce10.github.io/baRulho/reference/spectrum_correlation.md)
  : Measure frequency spectrum correlation
- [`tail_to_signal_ratio()`](https://marce10.github.io/baRulho/reference/tail_to_signal_ratio.md)
  : Measure reverberations as tail-to-signal ratio

## Built in datasets

Datasets included in baRulho

- [`master_est`](https://marce10.github.io/baRulho/reference/master_est.md)
  : Extended selection table of master acoustic data
- [`test_sounds_est`](https://marce10.github.io/baRulho/reference/test_sounds_est.md)
  : Extended selection table with re-recorded playbacks

## Additional functions

- [`attenuation()`](https://marce10.github.io/baRulho/reference/attenuation.md)
  : Estimate attenuation of sound pressure level
- [`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md)
  : Add synthetic noise
- [`noise_profile()`](https://marce10.github.io/baRulho/reference/noise_profile.md)
  : Measure full spectrum sound noise profiles
- [`spot_ambient_noise()`](https://marce10.github.io/baRulho/reference/spot_ambient_noise.md)
  : Find a segment of ambient noise to be used as reference
