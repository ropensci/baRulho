pkgname <- "baRulho"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "baRulho-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('baRulho')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("atmospheric_attenuation")
### * atmospheric_attenuation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: atmospheric_attenuation
### Title: Measure atmospheric attenuation and absorption of sound
### Aliases: atmospheric_attenuation

### ** Examples

{
# load example data
data("playback_est")

#' # remove ambient selections
playback_est <- playback_est[playback_est$signal.type != "ambient", ]

# measure atmospheric attenuation formula 1
atmospheric_attenuation(f = 20000, temp = 20, RH = 90, p = 88000, formula = 1)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("atmospheric_attenuation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("blur_ratio")
### * blur_ratio

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blur_ratio
### Title: Measure blur ratio in the time domain
### Aliases: blur_ratio

### ** Examples

{
# load example data
data("playback_est")

# remove ambient selections
playback_est <- playback_est[playback_est$signal.type != "ambient", ]

# using method 1
blur_ratio(X = playback_est)

# using method 2
blur_ratio(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blur_ratio", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("envelope_correlation")
### * envelope_correlation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: envelope_correlation
### Title: Measure amplitude envelope correlation
### Aliases: envelope_correlation

### ** Examples

{
# load example data
data("playback_est")

# remove ambient selections
playback_est <- playback_est[playback_est$signal.type != "ambient", ]

# method 1
envelope_correlation(X = playback_est)

# method 2
envelope_correlation(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("envelope_correlation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("excess_attenuation")
### * excess_attenuation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: excess_attenuation
### Title: Measure excess attenuation
### Aliases: excess_attenuation

### ** Examples

{
# load example data
data("playback_est")

# using method 1
excess_attenuation(X = playback_est)

# using method 2
excess_attenuation(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("excess_attenuation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("master_sound_file")
### * master_sound_file

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: master_sound_file
### Title: Create a master sound file
### Aliases: master_sound_file

### ** Examples

{
# load example data from warbleR
data(list = c("Phae.long1", "Phae.long2", "Phae.long3", "Phae.long4", 
"lbh_selec_table"))

# save sound files to temporary folder
writeWave(Phae.long1, file.path(tempdir(), "Phae.long1.wav"))
writeWave(Phae.long2, file.path(tempdir(), "Phae.long2.wav"))
writeWave(Phae.long3, file.path(tempdir(), "Phae.long3.wav"))
writeWave(Phae.long4, file.path(tempdir(), "Phae.long4.wav"))

# make an extended selection table
est <- selection_table(X = lbh_selec_table, extended = TRUE, confirm.extended = FALSE, 
path = tempdir())

# create master sound file
master.sel.tab <- master_sound_file(X = est, file.name = "example_master", 
dest.path = tempdir(), gap.duration = 0.3)

# the following code exports the selection table to Raven using Rraven package
# Rraven::exp_raven(master.sel.tab, path = tempdir(), file.name = "example_master_selection_table")
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("master_sound_file", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("snr")
### * snr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: snr
### Title: Measure attenuation as signal-to-noise ratio
### Aliases: snr

### ** Examples

{
# load example data
data("playback_est")

# using measure ambient noise reference selections 
snr(X = playback_est, mar = 0.05, noise.ref = 'custom')

# remove ambient selections
playback_est <- playback_est[playback_est$signal.type != "ambient", ]
# using margin for ambient noise of 0.05 and adjacent measure ambient noise reference
snr(X = playback_est, mar = 0.05, noise.ref = 'adjacent')
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("snr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spcc")
### * spcc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spcc
### Title: Measure spectrographic cross-correlation as a measure of signal
###   distortion
### Aliases: spcc

### ** Examples

{
# load example data
data("playback_est")

# method 1
spcc(X = playback_est, method = 1)

# method 2
spcc(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spcc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spcc_align")
### * spcc_align

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spcc_align
### Title: Align start and end of signal using spectrographic
###   cross-correlation
### Aliases: spcc_align

### ** Examples

{
# load example data
data("playback_est_unaligned")

# method 1
spcc_align(X = playback_est_unaligned)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spcc_align", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spectral_blur_ratio")
### * spectral_blur_ratio

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spectral_blur_ratio
### Title: Measure blur ratio in the frequency domain
### Aliases: spectral_blur_ratio

### ** Examples

{
# load example data
data("playback_est")

# remove ambient selections
playback_est <- playback_est[playback_est$signal.type != "ambient", ]

# using method 1
spectral_blur_ratio(X = playback_est)

# using method 2
spectral_blur_ratio(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spectral_blur_ratio", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("spectral_correlation")
### * spectral_correlation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: spectral_correlation
### Title: Measure frequency spectrum correlation
### Aliases: spectral_correlation

### ** Examples

{
# load example data
data("playback_est")

# remove ambient selections
pe <- playback_est[playback_est$signal.type != "ambient", ]

# method 1
spectral_correlation(X = pe)

# method 2
spectral_correlation(X = pe, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("spectral_correlation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
