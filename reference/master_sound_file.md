# Create a master sound file

`master_sound_file` creates a master sound file to be used in playback
experiments related to sound degradation.

## Usage

``` r
master_sound_file(
  X,
  file.name,
  dest.path = getOption("dest.path", "."),
  overwrite = FALSE,
  delay = 1,
  gap.duration = 1,
  amp.marker = 2,
  flim = c(0, 4),
  cex = 14,
  path = getOption("sound.files.path", ".")
)
```

## Arguments

- X:

  Object of class 'data.frame', 'selection_table' or
  'extended_selection_table' (the last 2 classes are created by the
  function
  [`selection_table`](https://marce10.github.io/warbleR/reference/selection_table.html)
  from the warbleR package) with the reference to the model sounds. Must
  contain the following columns: 1) "sound.files": name of the .wav
  files, 2) "selec": unique selection identifier (within a sound
  file), 3) "start": start time and 4) "end": end time of selections, 5)
  "bottom.freq": low frequency for bandpass and 6) "top.freq": high
  frequency for bandpass. An optional 'sound.id' column can be included
  to use a custom label for each sound in the output. This column must
  contain a unique id for each sound (labels cannot repeated). If not
  supplied the function will make it by combining the sound file and
  selection columns.

- file.name:

  Character string indicating the name of the sound file.

- dest.path:

  Character string containing the directory path where the sound file
  will be saved. If `NULL` (default) then the current working directory
  will be used instead.

- overwrite:

  Logical argument to determine if the function will overwrite any
  existing sound file with the same file name. Default is the current
  working directory.

- delay:

  Numeric vector of length 1 to control the duration (in s) of a silence
  gap at the beginning (and at the end) of the sound file. This can be
  useful to allow some time at the start of the playback experiment.
  Default is 1.

- gap.duration:

  Numeric vector of length 1 to control the duration (in s) of silence
  gaps to be placed in between sounds. Default is 1 s.

- amp.marker:

  Numeric vector of length 1 to use as a constant to amplify markers
  amplitude. This is useful to increase the amplitude of markers in
  relation to those of sounds, so it is picked up at further distances.
  Default is 2.

- flim:

  Numeric vector of length 2 to control the (approximate) frequency
  range in which the markers would be found. If `NULL` markers would
  span across the entire frequency range. Default is c(0, 4).

- cex:

  Numeric vector of length 1 indicating the font size for the start and
  end markers. Default is 14.

- path:

  Character string containing the directory path where the sound files
  are found. Only needed when 'X' is not an extended selection table. If
  not supplied the current working directory is used. Can be set
  globally for the current R session via the "sound.files.path" option
  (see [`options`](https://rdrr.io/r/base/options.html)).

## Value

A .wav file in 'path' as well as a data frame in the R environment with
the annotations (i.e. time position) of sounds in the master sound file
and an additional column 'sound.id' that provides a unique id for each
sound in the sound file. This is useful for identifying/labeling sounds
in test (re-recorded) sound files for downstream analyses.

## Details

The function is intended to simplify the creation of master sound files
for playback experiments in sound degradation studies. The function
clips sounds from sound files (or wave objects from extended selection
tables) and concatenates them in a single sound file. All clips are peak
normalized and rescaled to the maximal possible dynamic range. The
function also adds acoustic markers at the start and end of the playback
that can be used to time-sync test (re-recorded) sounds to facilitate
the streamlining of degradation quantification. There is no predefined
limit to the duration of the output master sound file, although the
creation of long files could be constrained by computer memory. As a
reference, master sound files of up to 10 min have been created in a
16GB RAM laptop computer.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

[`exp_raven`](https://rdrr.io/pkg/Rraven/man/exp_raven.html)

Other prepare acoustic data:
[`spot_ambient_noise()`](https://marce10.github.io/baRulho/reference/spot_ambient_noise.md),
[`synth_sounds()`](https://marce10.github.io/baRulho/reference/synth_sounds.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # load example data from warbleR
  data(list = c(
    "Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4",
    "lbh_selec_table"
  ))

  # save sound files to temporary folder
  writeWave(Phae.long1, file.path(tempdir(), "Phae.long1.wav"))
  writeWave(Phae.long2, file.path(tempdir(), "Phae.long2.wav"))
  writeWave(Phae.long3, file.path(tempdir(), "Phae.long3.wav"))
  writeWave(Phae.long4, file.path(tempdir(), "Phae.long4.wav"))

  # make an extended selection table
  est <- selection_table(
    X = lbh_selec_table, extended = TRUE,
    path = tempdir()
  )

  # create master sound file
  master.sel.tab <- master_sound_file(
    X = est, file.name = "example_master",
    dest.path = tempdir(), gap.duration = 0.3
  )

if (FALSE) { # \dontrun{
  # the following code exports the selection table to Raven
  # using the Rraven package
  Rraven::exp_raven(master.sel.tab, path = tempdir(),
  file.name = "example_master_selection_table")
  } # }
}
#> all selections are OK 
#> 

#> The file example_master.wav has been saved in the directory path '/tmp/Rtmp3xa5pl'
```
