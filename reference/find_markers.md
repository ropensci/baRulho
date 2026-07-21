# Find acoustic markers on test sound files

`find_markers` find acoustic markers on test (re-recorded) sound files
using spectrographic cross-correlation.

## Usage

``` r
find_markers(
  X,
  markers = c("start_marker", "end_marker"),
  test.files = NULL,
  path = getOption("sound.files.path", "."),
  pb = getOption("pb", TRUE),
  cores = getOption("mc.cores", 1),
  ...
)
```

## Arguments

- X:

  Object of class 'data.frame', 'selection_table' or
  'extended_selection_table' (the last 2 classes are created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the reference to the sounds in the
  master sound file. Must contain the following columns: 1)
  "sound.files": name of the .wav files, 2) "selec": unique selection
  identifier (within a sound file), 3) "start": start time, 4) "end":
  end time of selections and 5) "sound.id": unique identifier for each
  of the annotated sounds in 'X'. Columns for 'top.freq', 'bottom.freq'
  and 'channel' are optional. The acoustic start and end markers (added
  by
  [`master_sound_file`](https://marce10.github.io/baRulho/reference/master_sound_file.md))
  should be labeled as "start_marker" and "end_marker" respectively.
  Required.

- markers:

  Character vector with the name of the annotations (as in the column
  'sound.id') to be used as templates for cross-correlation. Default is
  `c("start_marker", "end_marker")`. Using more than one marker is
  recommended as the time difference between their position can be used
  to evaluate the precision of the detection (see 'Value' section).

- test.files:

  Character vector of length 1 with the name(s) of the test
  (re-recorded) file(s) in which to search for the marker(s). If not
  supplied all sound files in 'path' are used instead.

- path:

  Character string containing the directory path where test
  (re-recorded) sound files are found.

- pb:

  Logical argument to control if progress bar is shown. Default is
  `TRUE`. Can be set globally for the current R session via the "pb"
  option (see [`options`](https://rdrr.io/r/base/options.html)).

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default is 1
  (i.e. no parallel computing). Can be set globally for the current R
  session via the "mc.cores" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- ...:

  Additional arguments to be passed to
  [`template_correlator`](https://docs.ropensci.org/ohun/reference/template_correlator.html)
  for setting cross-correlation parameters (e.g. 'wl', 'ovlp', etc).

## Value

A data frame with test file names, marker id, maximum cross-correlation
score for each marker and the start and end where it was detected. If
two or more markers are used the function computes an additional column,
'time.mismatch', that compares the time difference between the two
markers in the test-files against that in the master sound file. In a
perfect detection the value must be 0.

## Details

The function takes a master sound file's reference data ('X') and finds
the position of acoustics markers ('markers' argument, included as
selections in 'X') in the re-recorded sound files. This is used to align
signals found in re-recorded sound files according to a master sound
file referenced in 'X'. The position of the markers is determined as the
highest spectrogram cross-correlation value for each marker using the
functions
[`template_correlator`](https://docs.ropensci.org/ohun/reference/template_correlator.html)
and
[`template_detector`](https://docs.ropensci.org/ohun/reference/template_detector.html).
**Make sure the master sound file (that referred to in 'X') is found in
the same folder than the re-recorded sound files**. Take a look at the
package vignette for information on how to incorporate this function
into a sound degradation analysis workflow. In cases in which markers
are not correctly detected editing test sound files to remove audio
segments with no target sounds (before the start marker and after the
end marker) can improve performance. Using a low 'hop.size' or window
length 'wl' (used internally by
[`template_correlator`](https://docs.ropensci.org/ohun/reference/template_correlator.html))
can help to improve precision. Other spectrogram types (argument 'type'
in
[`template_correlator`](https://docs.ropensci.org/ohun/reference/template_correlator.html))
can sometimes show better performance when markers are highly degraded.
If frequency range columns are included ('bottom.freq' and 'top.freq',
in kHz) cross-correlation will be run on those frequency ranges. All
templates must have the same sampling rate and both templates and
'files' (in which to find templates) must also have the same sampling
rate.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`manual_realign`](https://marce10.github.io/baRulho/reference/manual_realign.md);
[`auto_realign`](https://marce10.github.io/baRulho/reference/auto_realign.md);
[`align_test_files`](https://marce10.github.io/baRulho/reference/align_test_files.md);
[`master_sound_file`](https://marce10.github.io/baRulho/reference/master_sound_file.md)

Other test sound alignment:
[`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md),
[`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md),
[`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md),
[`plot_aligned_sounds()`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # set temporary directory
  td <- tempdir()

  # load example data
  data("master_est")

  # save example files in working director to recreate a case in which working
  # with sound files instead of extended selection tables.
  # This doesn't have to be done with your own data as you will
  # have them as sound files already.
  for (i in unique(test_sounds_est$sound.files)[1:2]) {
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(td, i))
  }

  # save master file
  writeWave(object = attr(master_est, "wave.objects")[[1]], file.path(td, "master.wav"))

  # set path and no progress bar in global options
  options(sound.files.path = td, pb = FALSE)

  # get marker position
  markers <- find_markers(X = master_est, test.files = unique(test_sounds_est$sound.files)[2])
}
```
