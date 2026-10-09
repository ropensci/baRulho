# Estimate attenuation of sound pressure level

`attenuation()` estimates atmospheric attenuation and atmospheric
absorption, calculating the geometric, atmospheric, and habitat
attenuation, as well as the overall expected attenuation (the sum of the
other three), based on temperature, relative humidity, atmospheric
pressure, and sound frequency.

## Usage

``` r
attenuation(
  frequency,
  dist0,
  dist,
  temp = 20,
  rh = 60,
  pa = 101325,
  hab.att.coef = 0.02
)
```

## Arguments

- frequency:

  Numeric vector of length 1 with frequency (in Hertz).

- dist0:

  Numeric vector of length 1 with distance (m) for the reference SPL.

- dist:

  Numeric vector of length 1 with distance (m) over which a sound
  propagates.

- temp:

  Numeric vector of length 1 with temperature (in Celsius). Default
  `20`.

- rh:

  Numeric vector of length 1 with relative humidity (in percentage).
  Default `60`.

- pa:

  Numeric vector of length 1 with atmospheric (barometric) pressure in
  Pa (standard: `101325`, default). Used for atmospheric attenuation.

- hab.att.coef:

  Attenuation coefficient of the habitat (in dB/kHz/m).

## Value

A `data.frame` with the geometric, atmospheric, and habitat attenuation
(in dB), as well as the combined attenuation.

## Details

Attenuation values are given in dB. The function is modified from
<http://www.sengpielaudio.com> and
<https://scikit-maad.github.io/generated/maad.spl.attenuation_dB.html#maad.spl.attenuation_dB>.

## References

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez,
A., & Rico-Guevara, A. (2025). Quantifying degradation in animal
acoustic signals with the R package baRulho. Methods in Ecology and
Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481

## See also

Other miscellaneous:
[`add_noise()`](https://marce10.github.io/baRulho/reference/add_noise.md),
[`noise_profile()`](https://marce10.github.io/baRulho/reference/noise_profile.md)

## Author

Marcelo Araya-Salas (<marcelo.araya@ucr.ac.cr>)

## Examples

``` r
{
  # measure attenuation
  attenuation(frequency = 2000, dist = 50, dist0 = 1)
}
#>   frequency dist geometric.attenuation atmospheric.attenuation
#> 1      2000   50               33.9794               0.4547757
#>   habitat.attenuation combined.attenuation
#> 1                1.96             36.39418
```
