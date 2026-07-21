# Find a segment of ambient noise to be used as reference

`spot_ambient_noise` finds a segment of ambient noise to be used as
reference by other functions.

## Usage

``` r
spot_ambient_noise(
  X,
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  path = getOption("sound.files.path", "."),
  length = NULL,
  ovlp = 0,
  fun = function(x) which.min(abs(x - mean(x)))
)
```

## Arguments

- X:

  Object of class 'data.frame', or 'selection_table' (a class are
  created by the function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the test sound files' annotations
  ('extended_selection_table' are not supported). Must contain the
  following columns: 1) "sound.files": name of the .wav files, 2)
  "selec": unique selection identifier (within a sound file), 3)
  "start": start time and 4) "end": end time of selections, 5)
  "bottom.freq": low frequency for bandpass, 6) "top.freq": high
  frequency for bandpass, 7) "sound.id": ID of sounds used to identify
  counterparts across distances/transects. 'selec' column values in 'X'
  cannot be duplicated within a sound file ('sound.files' column) as
  this combination is used to refer to specific rows.

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

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when 'X' is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the "sound.files.path" option
  (see [`options`](https://rdrr.io/r/base/options.html)).

- length:

  Numeric. Length (in s) of the segments to be used as ambient noise.
  Must be supplied. Default is `NULL`.

- ovlp:

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive segments. Default is 0. Can be set globally
  for the current R session via the "ovlp" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- fun:

  Function to be applied to select the segment to be used as ambient
  noise. It must be a function that takes a numeric vector (peak sound
  pressure level values for each candidate segment) and a single value
  with the index of the value to keep. Default is
  `function(x) which.min(abs(x - mean(x)))`.

## Value

An object similar to 'X' with one additional row for each sound file,
containing the selected 'ambient' reference.

## Details

This function finds a segment of ambient noise to be used as reference
by other functions. The function first finds candidate segments that do
not overlap with annotated sounds in 'X'. Then, it calculates the peak
sound pressure level (SPL) of each candidate segment and applies the
function supplied by the argument 'fun' to select a single segment. By
default 'fun' searches for the segment with the closest value to the
mean peak SPL across all candidate segments. Ambient noise annotations
are added as a new row in 'X'. Ambient noise annotations are used by the
functions
[`signal_to_noise_ratio`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md)
and
[`noise_profile`](https://marce10.github.io/baRulho/reference/noise_profile.md)
to determine background noise levels. Note that this function does not
work with annotations in 'extended_selection_table' format.

## References

\#' Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M.,
Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation
in animal acoustic signals with the R package baRulho. Methods in
Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
Araya-Salas, M., & Smith-Vidaurre, G. (2017). warbleR: An R package to
streamline analysis of animal acoustic signals. Methods in Ecology and
Evolution, 8(2), 184-191.

## See also

[`signal_to_noise_ratio`](https://marce10.github.io/baRulho/reference/signal_to_noise_ratio.md),
[`noise_profile`](https://marce10.github.io/baRulho/reference/noise_profile.md)

Other prepare acoustic data:
[`master_sound_file()`](https://marce10.github.io/baRulho/reference/master_sound_file.md),
[`synth_sounds()`](https://marce10.github.io/baRulho/reference/synth_sounds.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
# set temporary directory
td <- tempdir()  
# load example data
data("test_sounds_est")
########## save acoustic data (This doesn't have to be done 
# with your own data as you will have them as sound files already.)
 # save example files in working director 
for (i in unique(test_sounds_est$sound.files)[1:2]) {
 writeWave(object = attr(test_sounds_est, "wave.objects")[[i]],
           file.path(tempdir(), i))
}
test_sounds_df <- as.data.frame(test_sounds_est)
test_sounds_df <- test_sounds_df[test_sounds_df$sound.id != "ambient", ]
test_sounds_df <- 
 test_sounds_df[test_sounds_df$sound.files %in% 
  unique(test_sounds_est$sound.files)[1:2], ]
####
# closest to mean (default)
spot_ambient_noise(X = test_sounds_df, path = td, length = 0.12, ovlp = 20)

# min peak
spot_ambient_noise(X = test_sounds_df, path = td, length = 0.12, ovlp = 20, fun = which.min)
}
#>       sound.files selec    start      end bottom.freq top.freq sound.id
#> 1  10m_closed.wav     4 1.800045 2.000068       0.422    1.223    freq1
#> 2  10m_closed.wav     3 1.550023 1.750045       3.208    4.069    freq4
#> 3  10m_closed.wav     5 2.050068 2.250091       6.905    7.917    freq7
#> 4  10m_closed.wav     2 1.300000 1.500023       7.875    8.805    freq9
#> 5  10m_closed.wav     6 0.192000 0.312000       0.422    8.805  ambient
#> 6  30m_closed.wav     4 1.800045 2.000068       0.422    1.223    freq1
#> 7  30m_closed.wav     3 1.550023 1.750045       3.208    4.069    freq4
#> 8  30m_closed.wav     5 2.050068 2.250091       6.905    7.917    freq7
#> 9  30m_closed.wav     2 1.300000 1.500023       7.875    8.805    freq9
#> 10 30m_closed.wav     6 3.264000 3.384000       0.422    8.805  ambient
#>    transect distance
#> 1    closed       10
#> 2    closed       10
#> 3    closed       10
#> 4    closed       10
#> 5      <NA>       NA
#> 6    closed       30
#> 7    closed       30
#> 8    closed       30
#> 9    closed       30
#> 10     <NA>       NA
```
