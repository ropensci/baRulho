% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_blur_ratio.R
\name{plot_blur_ratio}
\alias{plot_blur_ratio}
\title{Plot blur ratio}
\usage{
plot_blur_ratio(
  X,
  type = c("envelope", "spectrum"),
  cores = getOption("mc.cores", 1),
  pb = getOption("pb", TRUE),
  env.smooth = getOption("env.smooth", 200),
  spec.smooth = getOption("spec.smooth", 5),
  res = 150,
  flim = c("-1", "+1"),
  hop.size = getOption("hop.size", 11.6),
  wl = getOption("wl", NULL),
  ovlp = getOption("ovlp", 70),
  palette = viridis::viridis,
  collevels = seq(-120, 0, 5),
  dest.path = getOption("dest.path", "."),
  path = getOption("sound.files.path", "."),
  colors = viridis::viridis(3),
  n.samples = if (type == "envelope") 100 else 500
)
}
\arguments{
\item{X}{The output of \code{\link{set_reference_sounds}} which is an object of class 'data.frame', 'selection_table' or 'extended_selection_table' (the last 2 classes are created by the function \code{\link[warbleR]{selection_table}} from the warbleR package) with the test sound files' annotations . Must contain the following columns: 1) "sound.files": name of the .wav files, 2) "selec": unique selection identifier (within a sound file), 3) "start": start time and 4) "end": end time of selections, 5) "bottom.freq": low frequency for bandpass, 6) "top.freq": high frequency for bandpass, 7) "sound.id": ID of sounds used to identify counterparts across distances and 8) "reference": identity of sounds to be used as reference for each test sound (row). See \code{\link{set_reference_sounds}} for more details on the structure of 'X'.}

\item{type}{Character vector of length 1 indicating the type of blur ratio to plot. The two options are 'envelope' (for regular blur ratio as in \code{\link{blur_ratio}}, default) and 'spectrum' (for spectrum blur ratio as in \code{\link{spectrum_blur_ratio}}).}

\item{cores}{Numeric vector of length 1. Controls whether parallel computing is applied by specifying the number of cores to be used. Default is 1 (i.e. no parallel computing). Can be set globally for the current R session via the "mc.cores" option (see \code{\link[base]{options}}).}

\item{pb}{Logical argument to control if progress bar is shown. Default is \code{TRUE}.  Can be set globally for the current R session via the "pb" option (see \code{\link[base]{options}}).}

\item{env.smooth}{Numeric vector of length 1 determining the length of the sliding window (in amplitude samples) used for a sum smooth for amplitude envelope calculation (used internally by \code{\link[seewave]{env}}). Default is 200.}

\item{spec.smooth}{Numeric vector of length 1 determining the length of the sliding window used for a sum smooth for power spectrum calculation (in kHz). Default is 5. Can be set globally for the current R session via the "spec.smooth" option (see \code{\link[base]{options}}).}

\item{res}{Numeric argument of length 1. Controls image resolution. Default is 150 (faster) although 300 - 400 is recommended for publication/presentation quality.}

\item{flim}{A numeric vector of length 2 indicating the highest and lowest frequency limits (kHz) of the spectrograms, as in \code{\link[seewave]{spectro}}. Default is \code{NULL}. Alternatively, a character vector similar to \code{c("-1", "1")} in which the first number is the value to be added to the minimum bottom frequency in 'X' and the second the value to be added to the maximum top frequency in 'X'. This is computed independently for each sound id so the frequency limit better fits the frequency range of the annotated signals. This is useful when test sounds show marked differences in their frequency ranges.}

\item{hop.size}{A numeric vector of length 1 specifying the time window duration (in ms). Default is 11.6 ms, which is equivalent to 512 wl for a 44.1 kHz sampling rate. Ignored if 'wl' is supplied. Can be set globally for the current R session via the "hop.size" option (see \code{\link[base]{options}}).}

\item{wl}{a vector with a single even integer number specifying the window length of the spectrogram, default is \code{NULL}. If supplied, 'hop.size' is ignored. Odd integers will be rounded up to the nearest even number. Can be set globally for the current R session via the "wl" option (see \code{\link[base]{options}}).}

\item{ovlp}{Numeric vector of length 1 specifying the percentage of overlap between two
consecutive windows, as in \code{\link[seewave]{spectro}}. Only used when plotting. Default is 70. Applied to both spectra and spectrograms on image files. Can be set globally for the current R session via the "ovlp" option (see \code{\link[base]{options}}).}

\item{palette}{A color palette function to be used to assign colors in the
plot, as in \code{\link[seewave]{spectro}}. Default is \code{\link[viridis]{viridis}}.}

\item{collevels}{Numeric vector indicating a set of levels which are used to partition the amplitude range of the spectrogram (in dB) as in \code{\link[seewave]{spectro}}. Default is \code{seq(-120, 0, 5)}.}

\item{dest.path}{Character string containing the directory path where the image files will be saved. If not supplied the current working directory will be used instead. Can be set globally for the current R session via the "dest.path" option (see \code{\link[base]{options}}).}

\item{path}{Character string containing the directory path where the sound files are found. Only needed when 'X' is not an extended selection table. If not supplied the current working directory is used. Can be set globally for the current R session via the "sound.files.path" option (see \code{\link[base]{options}}).}

\item{colors}{Character vector of length 4 containing the colors to be used for the color to identify the reference sound  (element 1),  the color to identify the test sound (element 2) and the color of blurred region (element 3).}

\item{n.samples}{Numeric vector of length 1 specifying the number of amplitude samples (or frequency bins if \code{spectrum = TRUE}) to use for representing power distributions. Default is 100 for \code{type = "envelope"} and 500 for \code{type = "spectrum"}. If null the raw power distribution is used (note that this can result in high RAM memory usage for large data sets).}
}
\value{
It returns 1 image file (in 'jpeg' format) for each blur ratio estimation, showing spectrograms of both sounds and the overlaid amplitude envelopes (or power spectra if \code{spectrum = TRUE}) as probability mass functions (PMF). Spectrograms are shown within the frequency range of the reference sound. It also returns the file path of the images invisibly.
}
\description{
\code{plot_blur_ratio} plots time and frequency blur ratio in sounds referenced in an extended selection table.
}
\details{
The function generates image files (in 'jpeg' format) for each possible blur ratio estimation in 'X'. The image files show the spectrograms of both sounds and the overlaid power distribution (either amplitude envelopes or power spectrum, see argument 'type') as probability mass functions (PMF). The output graphs highlight the mismatch between the compared distribution which represent the estimated blur ratio returned by either \code{\link{blur_ratio}} or \code{\link{spectrum_blur_ratio}}. Spectrograms are shown within the frequency range of the reference sound and also show dotted lines with the time (type = "envelope") or frequency range (type = "spectrum") in which energy distributions where computed.
}
\examples{
{
  # load example data
  data("test_sounds_est")

  # add reference to X
  X <- set_reference_sounds(X = test_sounds_est)

  # create plots
  plot_blur_ratio(X = X, dest.path = tempdir())
}

}
\references{
Dabelsteen, T., Larsen, O. N., & Pedersen, S. B. (1993). Habitat-induced degradation of sound signals: Quantifying the effects of communication sounds and bird location on blur ratio, excess attenuation, and signal-to-noise ratio in blackbird song. The Journal of the Acoustical Society of America, 93(4), 2206.

Araya-Salas, M., Grabarczyk, E. E., Quiroz-Oliva, M., Garcia-Rodriguez, A., & Rico-Guevara, A. (2025). Quantifying degradation in animal acoustic signals with the R package baRulho. Methods in Ecology and Evolution, 00, 1-12. https://doi.org/10.1111/2041-210X.14481
}
\seealso{
\code{\link{envelope_correlation}}, \code{\link{spectrum_blur_ratio}}, \code{\link{blur_ratio}}

Other quantify degradation: 
\code{\link{blur_ratio}()},
\code{\link{detection_distance}()},
\code{\link{envelope_correlation}()},
\code{\link{plot_degradation}()},
\code{\link{set_reference_sounds}()},
\code{\link{signal_to_noise_ratio}()},
\code{\link{spcc}()},
\code{\link{spectrum_blur_ratio}()},
\code{\link{spectrum_correlation}()},
\code{\link{tail_to_signal_ratio}()}
}
\author{
Marcelo Araya-Salas (\email{marcelo.araya@ucr.ac.cr})
}
\concept{quantify degradation}
