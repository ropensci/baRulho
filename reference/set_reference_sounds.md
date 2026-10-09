# Set reference for test sounds

`set_reference_sounds()` sets rows to be used as reference for each test
sound.

## Usage

``` r
set_reference_sounds(
  X,
  method = getOption("method", 1),
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  path = getOption("sound.files.path", ".")
)
```

## Arguments

- X:

  Object of class `data.frame`, `selection_table`, or
  `extended_selection_table` (the last 2 classes are created by
  [`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the **warbleR** package) with the test sound files' annotations.
  Must contain the following columns: 1) `sound.files`: name of the
  `.wav` files, 2) `selec`: unique selection identifier (within a sound
  file), 3) `start`: start time and 4) `end`: end time of selections, 5)
  `bottom.freq`: low frequency for bandpass, 6) `top.freq`: high
  frequency for bandpass, 7) `sound.id`: ID of sounds used to identify
  counterparts across distances (and transects, if more than 1), and 8)
  `distance`: distance (numeric) at which each test sound was
  re-recorded. A `transect` column labeling those sounds recorded in the
  same transect is required if `method = 2`. `X` can only have 1 copy
  for any given sound ID in a distance or a transect-distance
  combination (if the `transect` column is supplied). In addition,
  `selec` column values in `X` cannot be duplicated within a sound file
  (`sound.files` column), as this combination is used to refer to
  specific rows in the output `reference` column.

- method:

  Integer vector of length 1 to indicate the "experimental design" for
  measuring degradation. Two methods are available:

  - **`1`**: compare sounds (by `sound.id`) with their counterpart that
    was recorded at the closest distance to the source (e.g. compare
    sounds recorded at 5m, 10m, and 15m with their counterpart recorded
    at 1m). This is the default method. The function will try to use
    references from the same transect. However, if there is another test
    sound from the same `sound.id` at a shorter distance in other
    transects, it will be used as reference instead. This behavior aims
    to account for the fact that in this type of experiment, reference
    sounds are typically recorded at 1 m and at a single transect.

  - **`2`**: compare all sounds with their counterpart recorded at the
    distance immediately before within a transect (e.g. a sound recorded
    at 10m compared with the same sound recorded at 5m, then the sound
    recorded at 15m compared with the same sound recorded at 10m, and so
    on). The `transect` column in `X` is required.

  Can be set globally for the current R session via the `"method"`
  option (see [`options()`](https://rdrr.io/r/base/options.html)).

- cores:

  Numeric vector of length 1. Controls whether parallel computing is
  applied by specifying the number of cores to be used. Default `1`
  (i.e. no parallel computing). Can be set globally for the current R
  session via the `"mc.cores"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- pb:

  Logical argument to control if progress bar is shown. Default `TRUE`.
  Can be set globally for the current R session via the `"pb"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

## Value

An object similar to `X` with one additional column, `reference`, with
the ID of the sounds to be used as reference by degradation-quantifying
functions in downstream analyses. The ID is created as
`paste(X$sound.files, X$selec, sep = "-")`.

## Details

This function adds a `reference` column defining which sounds will be
used by other functions as reference. Two methods are available (see the
`method` argument description). For method 1, the function will attempt
to use re-recorded sounds from the shortest distance in the same
transect as reference. However, if there is another re-recorded sound
from the same `sound.id` at a shorter distance in other transects, it
will be used as reference instead. This behavior aims to account for the
fact that in this type of experiment, reference sounds are typically
recorded at 1 m and at a single transect. Note that if users want to
define their own reference sound, this can be set manually: `NA`s must
be used to indicate rows to be ignored, and references must be indicated
as the combination of the `sound.files` and `selec` column. For
instance, `"10m.wav-1"` indicates that the row in which the `selec`
column is `1` and the sound file is `"10m.wav"` should be used as
reference. The function also checks that the information in `X` is in
the right format so it won't produce errors in downstream analysis (see
the `X` argument description for details on format). The function will
ignore rows in which the `sound.id` column equals `"ambient"`,
`"start_marker"`, or `"end_marker"`.

## References

Araya-Salas, M., & Smith-Vidaurre, G. (2017). warbleR: An R package to
streamline analysis of animal acoustic signals. Methods in Ecology and
Evolution, 8(2), 184-191.

## See also

[`warbleR::check_sound_files()`](https://marce10.github.io/warbleR/reference/check_sound_files.html)
and
[`warbleR::check_sels()`](https://marce10.github.io/warbleR/reference/check_sels.html),
used internally to validate `X`.

Other quantify degradation:
[`blur_ratio()`](https://marce10.github.io/baRulho/reference/blur_ratio.md),
[`detection_distance()`](https://marce10.github.io/baRulho/reference/detection_distance.md),
[`envelope_correlation()`](https://marce10.github.io/baRulho/reference/envelope_correlation.md),
[`plot_blur_ratio()`](https://marce10.github.io/baRulho/reference/plot_blur_ratio.md),
[`plot_degradation()`](https://marce10.github.io/baRulho/reference/plot_degradation.md),
[`signal_to_noise_ratio()`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md),
[`spcc()`](https://marce10.github.io/baRulho/reference/spcc.md),
[`spectrum_blur_ratio()`](https://marce10.github.io/baRulho/reference/spectrum_blur_ratio.md),
[`spectrum_correlation()`](https://marce10.github.io/baRulho/reference/spectrum_correlation.md),
[`tail_to_signal_ratio()`](https://marce10.github.io/baRulho/reference/tail_to_signal_ratio.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data("test_sounds_est")

# save wav file examples
X <- test_sounds_est[test_sounds_est$sound.files != "master.wav", ]

# method 1
Y <- set_reference_sounds(X = X)

# method 2
Y <- set_reference_sounds(X = X, method = 2)
}
```
