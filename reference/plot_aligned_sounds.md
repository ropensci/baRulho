# Plot spectrograms to check test sound files alignment

`plot_aligned_sounds` plots spectrograms to visually inspect alignment
precision on test sound files.

## Usage

``` r
plot_aligned_sounds(
  X,
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 50),
  path = getOption("sound.files.path", "."),
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  collevels = seq(-120, 0, 5),
  palette = viridis::viridis,
  duration = 2,
  mar = 0.2,
  dest.path = getOption("dest.path", "."),
  flim = NULL,
  col = "white",
  width = 7,
  height = 4,
  res = 100,
  label = TRUE,
  fast.spec = FALSE,
  srt = 0,
  cex = 1,
  ...
)
```

## Arguments

- X:

  Object of class 'data.frame', 'selection_table' or
  'extended_selection_table' (the last 2 classes are created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the test sound files' annotations .
  Must contain the following columns: 1) "sound.files": name of the .wav
  files, 2) "selec": unique selection identifier (within a sound
  file), 3) "start": start time and 4) "end": end time of selections, 5)
  "bottom.freq": low frequency for bandpass, 6) "top.freq": high
  frequency for bandpass and 7) "sound.id": ID of sounds used to
  identify counterparts across distances. Each sound must have a unique
  ID within a distance.

- hop.size:

  A numeric vector of length 1 specifying the time window duration (in
  ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz
  sampling rate. Ignored if 'wl' is supplied. Can be set globally for
  the current R session via the "hop.size" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- wl:

  A numeric vector of length 1 specifying the window length of the
  spectrogram, default is NULL. Ignored if `bp = NULL`. If supplied,
  'hop.size' is ignored.

- ovlp:

  Numeric vector of length 1 specifying the percentage of overlap
  between two consecutive windows, as in
  [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html). Default
  is 0. Can be set globally for the current R session via the "ovlp"
  option (see [`options`](https://rdrr.io/r/base/options.html)).

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when 'X' is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the "sound.files.path" option
  (see [`options`](https://rdrr.io/r/base/options.html)).

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

- collevels:

  A numeric vector of length 3. Specifies levels to partition the
  amplitude range of the spectrogram (in dB). The more levels the higher
  the resolution of the spectrogram. Default is seq(-40, 0, 1).
  seq(-115, 0, 1) will produces spectrograms similar to other acoustic
  analysis software packages.

- palette:

  Color palette function for spectrogram. Default is
  [`viridis`](https://sjmgarnier.github.io/viridisLite/reference/viridis.html).
  See [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html) for more
  palettes. Palettes as
  [`gray.2`](https://rdrr.io/pkg/monitoR/man/specCols.html) may work
  better when `fast.spec = TRUE`.

- duration:

  A numeric vector of length 1. Specifies the overall duration of the
  clip that will be plotted. Notice that only the initial part of the
  test files are plotted as this is enough to tell the precision of the
  alignment.

- mar:

  numeric vector of length 1. Specifies the margins adjacent to the
  start of the first annotation to be included in the plot.

- dest.path:

  Character string containing the directory path where the image files
  will be saved. If not supplied the current working directory will be
  used instead. Can be set globally for the current R session via the
  "dest.path" option (see
  [`options`](https://rdrr.io/r/base/options.html)).

- flim:

  A numeric vector of length 2 indicating the highest and lowest
  frequency limits (kHz) of the spectrogram, as in
  [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html). Default is
  `NULL` which will plot spectrograms in the full frequency range (0 -
  nyquist frequency).

- col:

  Character string controlling the color of lines and sound ID labels.

- width:

  Numeric vector of length 1. Single value (in inches) indicating the
  width of the output image files. Default is 7.

- height:

  Numeric vector of length 1. Single value (in inches) indicating the
  height of the output image files. Default is 4.

- res:

  Numeric argument of length 1. Controls image resolution. Default is
  100 (faster) although 300 - 400 is recommended for
  publication/presentation quality.

- label:

  Logical to control if labels (from 'sound.id' column in 'X') are
  plotted. Default is `TRUE`.

- fast.spec:

  Logical. If `TRUE` then image function is used internally to create
  spectrograms, which substantially increases performance (much faster),
  although some options become unavailable, as collevels (amplitude
  scale). Default is `FALSE`.

- srt:

  Numeric argument of length 1. The rotation (in degrees) of the sound
  id labels. Default is 0.

- cex:

  Numeric argument of length 1controlling the size of sound id text
  labels. Default is 1.

- ...:

  Additional arguments to be passed to the internal spectrogram creating
  function for customizing graphical output. The function is a modified
  version of [`spectro`](https://rdrr.io/pkg/seewave/man/spectro.html),
  so it takes the same arguments.

## Value

Image files in jpeg format with spectrograms in the working directory,
one for each sound file in 'X'. It also returns the file path of the
images invisibly.

## Details

This functions aims to simplify the evaluation of the alignment of test
sound files from
[`align_test_files`](https://marce10.github.io/baRulho/reference/align_test_files.md).
The function creates a single spectrogram for each sound file (saved at
'dest.path'). Spectrograms include the first few seconds of the sound
files (controlled by 'duration') which is usually enough to tell the
precision of the alignment. The plots include vertical lines denoting
the start and end of each sound as well as the sound ID ('sound.id'
column in 'X'). Note that no plot is created in the R graphic device.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`manual_realign`](https://marce10.github.io/baRulho/reference/manual_realign.md);
[`auto_realign`](https://marce10.github.io/baRulho/reference/auto_realign.md);
[`find_markers`](https://marce10.github.io/baRulho/reference/find_markers.md);
[`align_test_files`](https://marce10.github.io/baRulho/reference/align_test_files.md)

Other test sound alignment:
[`align_test_files()`](https://marce10.github.io/baRulho/reference/align_test_files.md),
[`auto_realign()`](https://marce10.github.io/baRulho/reference/auto_realign.md),
[`find_markers()`](https://marce10.github.io/baRulho/reference/find_markers.md),
[`manual_realign()`](https://marce10.github.io/baRulho/reference/manual_realign.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data
  data("test_sounds_est")

  # plot (look into temporary working directory `tempdir()`)
  plot_aligned_sounds(X = test_sounds_est, dest.path = tempdir(), duration = 3, ovlp = 0)
}
#> The image files have been saved in the directory path '/tmp/Rtmp3xa5pl'
```
