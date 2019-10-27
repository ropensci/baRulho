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
nameEx("blur_ratio")
### * blur_ratio

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: blur_ratio
### Title: Measure blur ratio
### Aliases: blur_ratio

### ** Examples

{
# First set temporary folder
# setwd(tempdir())

data("playback_est")

# using margin for noise and method 1
blur_ratio(X = playback_est)

# using margin for noise and method 2
blur_ratio(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("blur_ratio", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("env_cor")
### * env_cor

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: env_cor
### Title: Measure amplitude envelope correlation
### Aliases: env_cor

### ** Examples

{
# First set temporary folder
# setwd(tempdir())

data("playback_est")

# method 1
env_cor(X = playback_est)

# method 2
env_cor(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("env_cor", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("excess_att")
### * excess_att

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: excess_att
### Title: Measure excess attenuation
### Aliases: excess_att

### ** Examples

{
# First set temporary folder
# setwd(tempdir())

data("playback_est")

# using margin for noise and method 1
excess_att(X = playback_est, mar = 0.05)

# using margin for noise and method 2
excess_att(X = playback_est, mar = 0.05, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("excess_att", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("xcorr_distortion")
### * xcorr_distortion

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: xcorr_distortion
### Title: Measure spectrographic cross-correlation as a measure of signal
###   distortion
### Aliases: xcorr_distortion

### ** Examples

{
# First set temporary folder
# setwd(tempdir())

data("playback_est")

# method 1
xcorr_distortion(X = playback_est)

# method 2
xcorr_distortion(X = playback_est, method = 2)
}




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("xcorr_distortion", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
