# Align test sound files

`align_test_files` aligns test (re-recorded) sound files.

## Usage

``` r
align_test_files(
  X,
  Y,
  path = getOption("sound.files.path", "."),
  by.song = TRUE,
  marker = NULL,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  ...
)
```

## Arguments

- X:

  object of class 'data.frame', 'selection_table' or
  'extended_selection_table' (the last 2 classes are created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the master sound file annotations. This
  should be the same data than that was used for finding the position of
  markers in
  [`find_markers`](https://marce10.github.io/baRulho/reference/find_markers.md).
  It should also contain a 'sound.id' column that will be used to label
  re-recorded sounds according to their counterpart in the master sound
  file.

- Y:

  object of class 'data.frame' with the output of
  [`find_markers`](https://marce10.github.io/baRulho/reference/find_markers.md).
  This object contains the position of markers in the re-recorded sound
  files. If more than one marker is supplied for a sound file only the
  one with the highest correlation score ('scores' column in 'X') is
  used.

- path:

  Character string containing the directory path where test
  (re-recorded) sound files are found.

- by.song:

  Logical argument to indicate if the extended selection table should be
  created by song (see 'by.song'
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  argument). Default is `TRUE`.

- marker:

  Character string to define whether a "start" or "end" marker would be
  used for aligning re-recorded sound files. Default is `NULL`.
  DEPRECATED.

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default is 1
  (i.e. no parallel computing). Can be set globally for the current R
  session via the "mc.cores" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- pb:

  Logical argument to control if progress bar is shown. Default is
  `TRUE`. Can be set globally for the current R session via the "pb"
  option (see [`options`](https://rdrr.io/r/base/options.html)).

- ...:

  Additional arguments to be passed to
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  for customizing extended selection table.

## Value

An object of the same class than 'X' with the aligned sounds from test
(re-recorded) sound files.

## Details

The function aligns sounds found in re-recorded sound files (referenced
in 'Y') according to a master sound file (referenced in 'X'). If more
than one marker is supplied for a sound file only the one with the
highest correlation score ('scores' column in 'X') is used. The function
outputs an 'extended selection table' by default.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`manual_realign`](https://marce10.github.io/baRulho/reference/manual_realign.md);
[`find_markers`](https://marce10.github.io/baRulho/reference/find_markers.md);
[`plot_aligned_sounds`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md)

Other test sound alignment:
[`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md),
[`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md),
[`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md),
[`plot_aligned_sounds()`](https://marce10.github.io/baRulho/reference/plot_aligned_sounds.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data(list = c("master_est", "test_sounds_est"))

  # save example files in working director to recreate a case in which working
  # with sound files instead of extended selection tables.
  # This doesn't have to be done with your own data as you will
  # have them as sound files already.
  for (i in unique(test_sounds_est$sound.files)[1:2]) {
    writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], 
              file.path(tempdir(), i))
  }

  # save master file
  writeWave(object = attr(master_est, "wave.objects")[[1]], 
        file.path(tempdir(), "master.wav"))

  # get marker position for the first test file
    markers <- find_markers(X = master_est,
    test.files = unique(test_sounds_est$sound.files)[1],
    path = tempdir())

  # align all test sounds
  alg.tests <- align_test_files(X = master_est, Y = markers, 
  path = tempdir())
}
#> computing correlations (step 1 of 0):
#> all selections are OK 
#> 
```
