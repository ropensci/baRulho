# Create synthetic sounds

`synth_sounds()` creates synthetic sounds that can be used for playback
experiments to understand the link between signal structure and its
transmission properties.

## Usage

``` r
synth_sounds(
  replicates = 1,
  frequencies,
  durations,
  nharmonics = 1,
  fm = FALSE,
  am = FALSE,
  am.amps = rep(c(1:4, 3:2), length.out = 11),
  mar = 0.05,
  seed = NULL,
  sig2 = 0.3,
  shuffle = FALSE,
  hrm.freqs = c(1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10),
  sampling.rate = 44.1,
  pb = getOption("pb", TRUE),
  freq.range = 2
)
```

## Arguments

- replicates:

  Numeric vector of length 1 indicating the number of replicates for
  each treatment combination. Default `1`. Useful for measuring
  variation in transmission parameters.

- frequencies:

  Numeric vector with the different frequencies (in kHz) to synthesize.
  A Brownian bridge motion stochastic process (`diff.fun == "BB"`) is
  used to simulate frequency modulation (see
  [`warbleR::simulate_songs()`](https://marce10.github.io/warbleR/reference/simulate_songs.html)).

- durations:

  Numeric vector with the different durations (in seconds) to
  synthesize.

- nharmonics:

  Numeric vector of length 1 specifying the number of harmonics to
  simulate. `1` indicates that only the fundamental frequency harmonic
  will be simulated.

- fm:

  Logical to control if both frequency modulated sounds and pure tones
  (i.e. non-modulated sounds) are synthesized. If `FALSE` (default),
  only pure tones are synthesized.

- am:

  Logical to control if both amplitude modulated sounds and
  non-modulated sounds are synthesized. If `FALSE` (default), only
  non-modulated sounds are synthesized.

- am.amps:

  Numeric vector with the relative amplitude for each time step to
  simulate amplitude modulation (only applied to the fundamental
  frequency). The default value (`rep(c(1:4, 3:2), length.out = 11)`)
  has 2 amplitude peaks (although only applied if `am = TRUE`).

- mar:

  Numeric vector with the duration of margins of silence around sounds,
  in seconds. Default `0.05`.

- seed:

  Numeric vector of length 1. This allows users to get the same results
  in different runs (using
  [`set.seed()`](https://rdrr.io/r/base/Random.html) internally).
  Default `NULL`.

- sig2:

  Numeric vector of length 1 defining the sigma value of the Brownian
  motion model (used for simulating frequency modulation). Default
  `0.3`.

- shuffle:

  Logical to control if the position of sounds is randomized. Having all
  sounds from the same treatment in a sequence can be problematic if an
  environmental noise masks them. Hence `shuffle` is useful to avoid
  having sounds from the same treatment next to each other. Default
  `FALSE`.

- hrm.freqs:

  Numeric vector with the frequencies of the harmonics relative to the
  fundamental frequency. The default values are
  `c(1/2, 1/3, 2/3, 1/4, 3/4, 1/5, 1/6, 1/7, 1/8, 1/9, 1/10)`.

- sampling.rate:

  Numeric vector of length 1. Sets the sampling frequency of the wave
  object (in kHz). Default `44.1`.

- pb:

  Logical argument to control if progress bar is shown. Default `TRUE`.
  Can be set globally for the current R session via the `"pb"` option
  (see [`options()`](https://rdrr.io/r/base/options.html)).

- freq.range:

  Numeric vector of length 1 with the frequency range around the
  simulated frequency in which signals will modulate. Default `2`, which
  means that sounds will range +/- 1 kHz around the target frequency.

## Value

An extended selection table, which can be input into
[`master_sound_file()`](https://marce10.github.io/baRulho/reference/master_sound_file.md)
to create the `.wav` file. The table contains columns for each of the
varying features, a `treatment` column (useful to tell the acoustic
features of each sound), and a `replicate` column indicating the
replicates for each `treatment`.

## Details

The function can add variation in signal structure in 5 features:

- **frequency**: continuous, argument `frequencies`.

- **duration**: continuous, argument `durations`.

- **harmonic structure**: binary (harmonics vs no-harmonics), arguments
  `nharmonics` and `hrm.freqs`.

- **frequency modulation**: variation in fundamental frequency across
  time. Binary (modulated vs non-modulated), arguments `fm` and `sig2`.

- **amplitude modulation**: variation in amplitude across time. Binary
  (modulated vs non-modulated), arguments `am` and `am.amps`.

Sounds for all possible combinations of the selected structure
dimensions will be synthesized. The output is an extended selection
table, which can be input into
[`master_sound_file()`](https://marce10.github.io/baRulho/reference/master_sound_file.md)
to create the `.wav` file. The function uses
[`warbleR::simulate_songs()`](https://marce10.github.io/warbleR/reference/simulate_songs.html)
internally for synthesizing individual sounds. A Brownian bridge motion
stochastic process (`diff.fun == "BB"`) is used to simulate frequency
modulation. The output table contains columns for each of the varying
features and a `treatment` column (useful to tell sounds from the same
combination of features apart when using replicates).

## References

Araya-Salas, M., & Smith-Vidaurre, G. (2017). warbleR: An R package to
streamline analysis of animal acoustic signals. Methods in Ecology and
Evolution, 8(2), 184-191.

## See also

[`warbleR::simulate_songs()`](https://marce10.github.io/warbleR/reference/simulate_songs.html),
used internally to synthesize individual sounds.

Other prepare acoustic data:
[`master_sound_file()`](https://marce10.github.io/baRulho/reference/master_sound_file.md),
[`spot_ambient_noise()`](https://marce10.github.io/baRulho/reference/spot_ambient_noise.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
if (FALSE) { # \dontrun{

synthetic_est <- synth_sounds(
  mar = 0.01,
  frequencies = c(1, 2, 3, 5),
  durations = 0.1,
  fm = TRUE,
  am = TRUE,
  nharmonics = 4,
  shuffle = TRUE,
  replicates = 3
)
} # }
```
