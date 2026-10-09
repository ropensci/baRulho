# Plot spectrograms to check test sound files alignment

`manual_realign()` plots spectrograms to visually inspect, and
interactively adjust, alignment precision on test sound files.

## Usage

``` r
manual_realign(
  X,
  Y,
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 0),
  path = getOption("sound.files.path", "."),
  collevels = seq(-120, 0, 5),
  palette = viridis::viridis,
  duration = 2,
  mar = 0.2,
  step.lengths = c(5, 30),
  flim = NULL,
  label.col = "white",
  ext.window = TRUE,
  width = 10,
  height = 5,
  srt = 0,
  cex = 1,
  fast.spec = TRUE,
  marker = "start_marker",
  grid = 0.2,
  ...
)
```

## Arguments

- X:

  Object of class `data.frame`, `selection_table`, or
  `extended_selection_table` (the last 2 classes are created by
  [`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the **warbleR** package) with the test sound files' annotations
  (typically the output of
  [`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md))
  to be aligned. Must contain the following columns: 1) `sound.files`:
  name of the `.wav` files, 2) `selec`: unique selection identifier
  (within a sound file), 3) `start`: start time and 4) `end`: end time
  of selections, 5) `bottom.freq`: low frequency for bandpass, 6)
  `top.freq`: high frequency for bandpass, and 7) `sound.id`: ID of
  sounds used to identify counterparts across distances. Each sound must
  have a unique ID within a given distance.

- Y:

  Object of class `data.frame`, `selection_table`, or
  `extended_selection_table` (the last 2 classes are created by
  [`warbleR::selection_table()`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the **warbleR** package) with the master sound file annotations.
  This should be the same data used for finding the position of markers
  in
  [`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md).
  It should also contain a `sound.id` column.

- hop.size:

  Numeric vector of length 1 specifying the time window duration (in
  ms). Default `11.6` ms, which is equivalent to 512 `wl` for a 44.1 kHz
  sampling rate. Ignored if `wl` is supplied. Can be set globally for
  the current R session via the `"hop.size"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- wl:

  A vector with a single even integer number specifying the window
  length of the spectrogram. Default `NULL`. If supplied, `hop.size` is
  ignored. Odd integers will be rounded up to the nearest even number.
  Can be set globally for the current R session via the `"wl"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- ovlp:

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive windows, as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `0`. Can be set globally for the current R session via the
  `"ovlp"` option (see
  [`options()`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when `X` is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the `"sound.files.path"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- collevels:

  Numeric vector of length 3. Specifies levels to partition the
  amplitude range of the spectrogram (in dB). The more levels, the
  higher the resolution of the spectrogram. Default `seq(-120, 0, 1)`.

- palette:

  Color palette function for the spectrogram. Default
  [`viridis::viridis()`](https://sjmgarnier.github.io/viridis/reference/reexports.html).
  See
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html)
  for more palettes. Palettes such as
  [`monitoR::gray.2`](https://rdrr.io/pkg/monitoR/man/specCols.html) may
  work better when `fast.spec = TRUE` (an argument that can be passed to
  the internal spectrogram function using `...`).

- duration:

  Numeric vector of length 1. Specifies the overall duration of the clip
  that will be plotted. Notice that only the initial part of the test
  files is plotted, as this is usually enough to tell the precision of
  the alignment.

- mar:

  Numeric vector of length 1. Specifies the minimum margins adjacent
  (before and after) to the start of the marker used for checking
  alignments (see the `marker` argument). Default `0.2`.

- step.lengths:

  Numeric vector of length 2 indicating the time length (in ms) of short
  (`min(step.lengths)`) and long (`max(step.lengths)`) steps for
  manually aligning spectrograms. Default `c(5, 30)`.

- flim:

  Numeric vector of length 2 indicating the highest and lowest frequency
  limits (kHz) of the spectrogram, as in
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html).
  Default `NULL`, which will plot spectrograms in the full frequency
  range (0 - nyquist frequency).

- label.col:

  Character string controlling the color of lines and sound ID labels.

- ext.window:

  Logical. If `TRUE`, an external graphic window is used. Dimensions can
  be set using the `width` and `height` arguments. Default `TRUE`.

- width:

  Numeric vector of length 1. Single value (in inches) indicating the
  width of the output image files. Default `10`.

- height:

  Numeric vector of length 1. Single value (in inches) indicating the
  height of the output image files. Default `5`.

- srt:

  Numeric argument of length 1. The rotation (in degrees) of the sound
  ID labels. Default `0`.

- cex:

  Numeric argument of length 1 controlling the size of sound ID text
  labels. Default `1`.

- fast.spec:

  Logical. If `TRUE`, the `image` function is used internally to create
  spectrograms, which substantially increases performance (much faster),
  although some options become unavailable, such as `collevels`
  (amplitude scale). Default `FALSE`.

- marker:

  Character string with the name of the marker to be used as the main
  reference for checking/adjusting time alignments. Default
  `"start_marker"`. Note that this can take any of the sound IDs in
  `Y$sound.id`.

- grid:

  Numeric vector of length 1 controlling the spacing between vertical
  lines on the spectrogram. Default `0.2` s. Use `0` to remove the grid.

- ...:

  Additional arguments to be passed to the internal spectrogram-creating
  function for customizing graphical output. The function is a modified
  version of
  [`seewave::spectro()`](https://rdrr.io/pkg/seewave/man/spectro.html),
  so it takes the same arguments.

## Value

Creates a multipanel graph with spectrograms of master and test sound
files in which users can interactively adjust their alignment in time.
Returns an object similar to the input object `X`, in which the start
and end of the sounds have been adjusted.

## Details

This function allows the interactive adjustment of the alignment of test
sound files produced by
[`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md).
The function generates a multipanel graph with the spectrogram of the
master sound file on top of that from test sound files, highlighting the
position of corresponding test sounds on both in order to simplify
assessing and adjusting their alignment. Spectrograms include the first
few seconds of the sound files (controlled by `duration`), which is
usually enough to tell the precision of the alignment. The lower
spectrogram shows a series of "buttons" that users can click on to
control if the test sound file spectrogram (lower panel) needs to be
moved to the left (`"<"`) or right (`">"`). Users can also reset the
spectrogram to its original position (`"reset"`), move on to the next
sound file in `X` (test sound file annotations), or stop the process
(stop button).

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md),
for automatic (non-interactive) realignment;
[`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md)
and
[`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md),
used upstream to produce the input for this function.

Other test sound alignment:
[`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md),
[`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md),
[`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md),
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
  writeWave(object = attr(test_sounds_est, "wave.objects")[[i]], file.path(tempdir(), i))
  }

  # save master file
  writeWave(object = attr(master_est, "wave.objects")[[1]], file.path(tempdir(), "master.wav"))

  # get marker position
  markers <- find_markers(X = master_est, test.files = unique(test_sounds_est$sound.files)[2], 
  path = tempdir())

  # align all test sounds
  alg.tests <- align_test_files(X = master_est, Y = markers, path = tempdir())

  # add error to alignment
  lag <- (as.numeric(as.factor(alg.tests$sound.files)) - 2) / 30
  alg.tests$start <- alg.tests$start + lag
  alg.tests$end <- alg.tests$end + lag

  if(interactive()){
  realigned_est <- manual_realign(X = alg.tests, Y = master_est, duration = 2,
  ovlp = 50, hop.size = 14, collevels = seq(-140, 0, 5), palette = viridis::mako,
  ext.window = FALSE)
 }
}
```
