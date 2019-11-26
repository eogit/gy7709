#!/usr/bin/Rscript

#######################################
# Generic package for harmonic analysis in R
#######################################
#
# written by Heiko Balzter, copyright 2014
# contact the author: hb91@le.ac.uk
#
# includes modified functions from source: http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html
#
#######################################

# t is a time index from 0 to n-1 where n is the number of measurements in the time-series
# xt is the vector of time points at which the measurements were taken
# The amplitude of a wave is defined as half the height from the maximum to the minimum point.
# The phase of the wave is defined as the angle by which the sine wave is delayed to its first peak.
# A harmonic term is defined by how many complete waves it has within the defined time series, from start to end,
#    i.e. harmonic term 2 has two full waves (two maxima and two minima) within the time series.
# The fundamental period is the period between the first sample and the last.
# The acquisition frequency is the number of measurements between two successive units of time. 
# The fundamental frequency f_0 is 1/N where N is the number of time steps. 
# The frequencies of the wave components must be integer multiples of the fundamental frequency.
# f_0 is called the first harmonic, the second harmonic is 2*f_0, the third is 3*f_0, etc.

get.trajectory <- function(X.k,xt,acq.freq) {
# Inverse Fourier Transform: 
# Returns the x.n time series for a given time sequence (xt) and a vector with the amount of frequencies k in the signal (X.k)
  n   <- length(xt)
  i   <- complex(real = 0, imaginary = 1)
  x.n <- rep(0,n)
  ks  <- 0:(length(X.k)-1)
  for(j in 0:(n-1)) { # compute each time point x_n based on freqs X.k
    x.n[j+1] <- sum(X.k * exp(i*2*pi*ks*j/n)) / n
  }
  x.n * n
}

plot.fourier <- function(fourier.series, f.0, xt, ...) {
# plot a Fourier series
# ***** This function has been verified. *****
  w <- 2*pi*f.0
  trajectory <- sapply(xt, function(t) fourier.series(t,w))
  plot(xt/length(xt), trajectory, type="l", xlab="time", ylab="f(t)"); abline(h=0,lty=3)
}

convert.fft <- function(x.k, acq.freq=1) {
# convert a FFT to amplitude and phase
# x.k is the vector of complex points to convert
  n <- length(Re(x.k)) # number of points
  x.k <- x.k / n # normalize
  distance.center <- function(c)signif( Mod(c),        4)
  angle           <- function(c)signif( 180*Arg(c)/pi, 3)
  df <- data.frame(cycle    = 0:(n-1),
                   freq     = 0:(n-1) / acq.freq,
                   t        = n / 0:(n-1) / acq.freq, # in time units, not sequential units
                   ampl     = sapply(x.k, distance.center) * 2 * n,
                   phase    = sapply(x.k, angle))
  df
}

plot.frequency.spectrum <- function(X.k, acq.freq=acq.freq, col = 1, lwd = 2, pch = "+", ...) {
# plot a frequency spectrum of a given X_k
  xax <- (0:(length(X.k)-1)) / length(X.k) * acq.freq
  xlimits <- c(0, max(xax)/2)
  plot.data  <- cbind( xax, 2 * Mod(X.k))
  plot(plot.data, t="h", main="Periodogram", xlab="Frequency", ylab="Power spectral density",
       col = col, lwd = lwd, pch = pch, xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
}

plot.harmonic <- function(xk, i, xt, acq.freq, mar=c(1,1,1,1),
  col = 3, lwd = 2, pch = "+", cex.lab = 1, cex.axis = 1, cex.main = 1, cex.sub = 1, ...) {
# plot.harmonic() plots the i-th harmonic on the current plot
# xk: the frequencies computed by the FFt
#  i: which harmonic
# xt: the sampling time points
# acq.freq: the acquisition rate
  xk.h <- rep(0,length(xk))
  xk.h[i+1] <- xk[i+1] # i-th harmonic
  harmonic.trajectory <- 2 * get.trajectory(xk.h, xt, acq.freq=acq.freq)
  points(xt, Re(harmonic.trajectory), type="l", mar=mar,
       col = col, lwd = lwd, pch = pch, cex.lab = cex.lab, cex.axis = cex.axis,
       cex.main = cex.main, cex.sub = cex.sub)
}

get.harmonic <- function(xk, i, xt, acq.freq) {
# Get the values that define the i-th harmonic term.
# xk: the frequencies computed by the FFt
#  i: which harmonic term(s)
# xt: the sampling time points
# acq.freq: the acquisition rate
  xk.h <- rep(0,length(xk))
  xk.h[i+1] <- xk[i+1] # i-th harmonic
  harmonic.trajectory <- 2 * get.trajectory(xk.h, xt, acq.freq=acq.freq)
  Re(harmonic.trajectory)
}

harmonic <- function(xt, x, acq.freq, N, alpha, detrend, which, test, ...) {
# core harmonic analysis function
# xt = a vector of time steps in units of s,min, hr or other time units, does not have to be integers
# x = a vector of time-series observations with the same length as xt
# N = number of the harmonic terms to be included, starting with largest amplitude
# acq.freq = number of measurements between two successive units of time. 
# alpha = type I error probability for statistical significance testing (default 0.05 or 5%)
# detrend = TRUE or FALSE, if TRUE (default) then linear detrending is applied to the data
# which = "strongest": the N strongest harmonic terms are included in the model (default)
#       = "first": the first N terms are included, or
#       = "all": all harmonics are included.
# test = "bonferroni" adjusts type I error by the number of tests N; "holm" adjusts the type I error by N+1-k where k=1:N 
#
# How to understand the harmonic terms: 
# Cycle means the number of waves in the time series, i.e. cycle = 9 the wave fits 9 times into the length of the data
#      which is the annual cycle for a 9-year series. 
# Freq is the position in the frequency domain (periodogram) from 0 to 0.5. 
# t is the position in the time domain from 0 to n-1 where n is the number of measurements. 
# Amplitude is the strength of the wave and 
# phase is the delay of the wave in degrees (0-360).

  dig.aov <- 4 # number of significant digits for ANOVA table
  if (missing(acq.freq)) acq.freq <- 1
  if (missing(N)) N <- 20
  if (missing(alpha)) alpha <- 0.05
  if (missing(detrend)) detrend <- TRUE
  if (missing(which)) which <- "strongest"
  if (!(which %in% c("strongest", "first", "all"))) which <- "strongest"
  if (missing(test)) test <- "bonferroni"
  if (!(test %in% c("bonferroni", "holm"))) test <- "bonferroni"
  t <- 0:(n-1)
  # detrending
  n <- length(x)
  if (detrend) {
    trend <- lm(x~xt) # linear model
    cat("Linear detrending result:\n")
    print(summary(trend))
    detrended.x <- trend$residuals
    } else {
    cat("No detrending.\n\n")
    detrended.x <- x
    trend <- "No detrending"
    }

  detrended.x.k <- fft(detrended.x) / n
  windows()
  plot.frequency.spectrum(detrended.x.k, acq.freq=acq.freq)

  # Calculate the amplitude and phase angle for the N harmonic terms
  # Cycle 9 means that the harmonic wave repeats 9 times over the time series
  # Note that cycle is indexed from 0 and ndx from 1
  x.fft <- convert.fft(detrended.x.k, acq.freq)
  nx <- length(x.fft$cycle)
  # you can get the components of the table from:
  #   x.fft$cycle[ndx]
  #   x.fft$freq[ndx]
  #   x.fft$t[ndx]
  #   x.fft$ampl[ndx]
  #   x.fft$phase[ndx]

  # find the N strongest harmonics
  if (which=="strongest") {
    ndx <- order(x.fft$ampl[1:(nx/2)], decreasing = T)[1:N]
    #print(cbind(ndx, x.fft$ampl[ndx]))
    cat(paste(N, "strongest harmonic terms:\n", sep=" "))
    write.table(round(x.fft[ndx,],4), quote = F, sep = "\t", row.names=F)
    cat("\n")
    }
  if (which=="first") {
    ndx <- 2:(N+1)
    #print(cbind(ndx, x.fft$ampl[ndx]))
    cat(paste(N, "first harmonic terms:\n", sep=" "))
    write.table(round(x.fft[ndx,],4), quote = F, sep = "\t", row.names=F)
    cat("\n")
    }
  if (which=="all") {
    ndx <- 1:n
    cat(paste("All harmonic terms:\n", sep=" "))
    write.table(round(x.fft[ndx,],4), quote = F, sep = "\t", row.names=F)
    cat("\n")
    }

  # test for significance of the individual harmonic terms using the F test
  # Reference: Hartley, H. O. (1949): Tests of Significance in Harmonic Analysis. Biometrika, 36, 194-201.
  # time dimension ts is in arbitrary units, with acq.freq measurements between two successive units 1 and 2, say
  # without loss of generality, for the purpose of the significance testing we treat the 
  #     time dimension as steps of 1, 2, ..., n
  # we do this by adjusting the time index tn <- xt * acq.freq

  if (which=="strongest" || which=="first") {
  bonferroni <- alpha / N   # adjusted alpha probability by the number of comparisons, Bonferroni correction
  holm <- alpha / seq(N, 1, -1)   # adjusted alpha probability, Bonferroni/Holm correction
  gamma = 2*pi/n # in Hartley's paper, but only for time steps of 1
  a0 <- mean(detrended.x)
  # work out the mean squares (MSQ) of each harmonic term for ANOVA table
  ssq <- rep(0, N) # SSQ components 
  df <- rep(2, N) # degrees of freedom
  msq <- rep(0, N) # MSQ components = SSQ / df
  a <- rep(0, N) # ai
  b <- rep(0, N) # bi
  f <- rep(0, N) # F values for each harmonic term
  p <- rep(0, N) # p values for each harmonic term
  sig <- rep("n.s.", N) # significance
  for (i in 1:N) {
    a[i] <- 2/n * sum(detrended.x * cos(x.fft$cycle[ndx[i]] * t * gamma)) # note that we use t here and not xt, see above
    b[i] <- 2/n * sum(detrended.x * sin(x.fft$cycle[ndx[i]] * t * gamma))
  }
  # calculate SSQ terms
  ssq <- n/2*(a^2+b^2)
  # calculate MSQ terms
  msq <- ssq/df
  # Calculate the residual MSQ variance component:
  rssq <- sum((detrended.x-a0)^2) - n/2 * sum(a^2+b^2)
  # The total df is n-1. The residual df is the total n – 2N - 1.
  rdf <- n-2*N-1 # residual df
  rmsq <- rssq/rdf
  # Work out the F values:
  f <- msq / rmsq
  # Each harmonic term has 2 degrees of freedom since it is characterised by 2 parameters ai and bi. 
  # The F ratio is calculated by dividing the MSQ of each harmonic term by the residual MSQ. It should be compared to the F distribution for 2,11 degrees of freedom for the 5%/m point, assuming type I error is controlled at 5%.
  p <- pf(f, df1=2, df2=n-N*2-1, lower=FALSE)
  # rounding for pretty printing:
  ssq <- round(ssq, 2)
  msq <- round(msq, 2)
  f <- round(f, 1)
  p <- round(p, 5)
  bonferroni <- round(bonferroni, 5)
  holm <- round(holm, dig.aov)
  rssq <- round(rssq, dig.aov)
  rmsq <- round(rmsq, dig.aov)
  # now merge all into a data.frame
  if (test=="bonferroni") {
    sig[p<bonferroni] <- "*"
    x.aov <- as.data.frame(cbind(x.fft$cycle[ndx], ssq, df, msq, f, p, bonferroni, sig))
    names(x.aov) <- c("cycle", "SSQ", "df", "MSQ", "F", "p", "pBonf", "Sig")
    }
  if (test=="holm") {
    sig[p<holm] <- "*"
    for (i in 1:(N-1)) if (p[i]>=holm[i]) sig[(i+1):N] <- rep("n.s.", N-i)
    x.aov <- as.data.frame(cbind(x.fft$cycle[ndx], ssq, df, msq, f, p, holm, sig))
    names(x.aov) <- c("cycle", "SSQ", "df", "MSQ", "F", "p", "pHolm", "Sig")
    }
  # print it
  cat("ANOVA table for the selected harmonic terms:\n")
  write.table(format(x.aov, trim = FALSE, justify = "right"), quote = F, sep = "\t", row.names=F)
  write.table(format(cbind("Res.", rssq, rdf, rmsq), trim = FALSE, justify = "right"), quote = F, sep = "\t", row.names=F, col.names=F)
  # cat("Res.", round(rssq, dig.aov), rdf, round(rmsq, dig.aov), "\n")
  cat("\n")
  # print(x.aov)
  x.aov <- merge.data.frame(x.aov, data.frame(c(NA, round(rssq,dig.aov), rdf, round(rmsq,dig.aov), NA,NA,NA,NA), row.names = names(x.aov)))

  if (which=="strongest") { # select only the significant harmonic terms
    ndxs <- ndx[sig=="*"]
    N <- length(ndxs)
    cat("\nRetaining ", N, "significant harmonic terms.\n")
    }
  if (which=="first") { # select the first N harmonic terms
    ndxs <- ndx
    cat("\nRetaining the first", N, " harmonic terms.\n")
    }

  # plot detrended time series and overlay the individual N significant harmonics with the largest amplitudes:
  # only plot up to 40 harmonics
  x.n <- get.trajectory(detrended.x.k, xt, acq.freq)  # create time wave from detrended data (if detrending is switched on)
  windows()
  par(mfrow = c(1,1))
  plot(xt, Re(x.n), type="l", lwd=1)
  abline(h=0,lty=2)
  for (i in 1:min(40, N)) plot.harmonic(detrended.x.k, ndxs[i], xt, acq.freq, col=i+1)

  # Now plot detrended time series and the composite of the first significant N harmonics:
  windows()
  plot(xt, Re(x.n), type="l", lwd=1)
  abline(h=0, lty=2)
  wave <- get.harmonic(detrended.x.k, ndxs, xt, acq.freq)
  lines(xt, wave, col="red")

  # And now plot add the trend back on to the composite of the first N harmonics:
  if (detrend) {
    windows()
    plot(xt, x, type="l",lwd=1)
    abline(trend)
    wave <- wave + trend$coef[1] + trend$coef[2] * xt
    lines(xt, wave, col="red")
    }

  # plot residuals
    windows()
    plot(xt, x - wave, type="p", pch="+")
  }
  else
  { # if which == "all"
  ndxs <- 2:n # in case all harmonics will be included except term 0
  N <- n
  x.aov <- "No ANOVA available. All harmonics are included."
  wave <- x
  }

  # return these components:
  harm <- list(xt)
  names(harm) <- "xt"
  harm$lm <- trend
  harm$detrended <- detrended.x
  harm$Nsig <- N
  harm$ndx <- ndx
  harm$frequency.spectrum <- detrended.x.k
  harm$Amp <- x.fft$ampl[ndx]
  harm$Ph <- x.fft$phase[ndx]
  harm$aov <- x.aov
  harm$fitted.values <- wave
  harm$residuals <- x - wave
  harm$call <- match.call()
  class(harm) <- "harmonic"

  # return results  
  harm
}
##################
# end of functions
##################

###############
# verification:
###############

# set working directory
setwd("C:/Users/localadmin1/Desktop/_todo/Balzter_Lake_Balaton_phenology/harmonic analysis")

# set the handling of warning messages. 
# If warn is negative all warnings are ignored. 
# If warn is zero (the default) warnings are stored until the top–level function returns. 
#     If 10 or fewer warnings were signalled they will be printed otherwise a message saying how many were signalled. An object called last.warning is created and can be printed through the function warnings. If warn is one, warnings are printed as they occur. 
# If warn is two or larger all warnings are turned into errors.
options(warn=0)

# close any open graphics windows
for (i in 1:10) dev.off()

# example from web site:
acq.freq <- 100                    # data acquisition (sample) frequency (Hz)
time     <- 6                      # measuring time interval (seconds)
xt       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
f.0 <- 1/time
n <- length(xt)
dc.component <- 1
component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
component.delay <- c(0,0,0)         # delay of signal components (radians)
component.strength <- c(1.5,.5,.75) # strength of signal components
f   <- function(t,w) { 
  dc.component + 
  sum( component.strength * sin(component.freqs*w*t + component.delay)) 
}
windows()
plot.fourier(f,f.0,xt=xt)
w <- 2*pi*f.0
trajectory <- sapply(xt, function(t) f(t,w))
windows()
plot(xt, trajectory, type="l")
X.k <- fft(trajectory)/n                   # find all harmonics with fft()
windows()
plot.frequency.spectrum(X.k, acq.freq=acq.freq)
x.n <- get.trajectory(X.k,xt,acq.freq)
windows()
plot(xt, Re(x.n), type="l"); abline(h=0,lty=3)
points(xt, Re(trajectory),col="red",type="l") # compare with original
convert.fft(X.k, acq.freq=acq.freq)[1:40,]

####################################
# This is how you call the function:
####################################

# with a cosine wave for checking:
time     <- 50*pi                          # measuring time interval (in time units)
acq.freq <- 20                             # data acquisition frequency (Hz), how many measurements
t <- seq(0,time-1/acq.freq,1/acq.freq)     # vector of sampling time-points (s) 
x <- cos(t)*10
n <- length(x)
plot(t,x, type="l")
x.harm <- harmonic(t, x, acq.freq=acq.freq, N=5, alpha=0.05, detrend=FALSE, test="bonferroni")

# cosine wave with phase shift
time     <- 50*pi                         # measuring time interval (in time units)
acq.freq <- 20                           # data acquisition frequency (Hz), how many measurements
t <- seq(0,time-1/acq.freq,1/acq.freq)   # vector of sampling time-points (s) 
# Shift the phase by -90 degrees
ph2 = -90 * pi / 180  # convert to radians
x <- cos(t+ph2)*10
n <- length(x)
plot(t,x, type="l")
x.harm <- harmonic(t, x, acq.freq=acq.freq, N=5, alpha=0.05, detrend=FALSE, test="holm")

##############################
# let's define an oscillating pattern with two waves of amplitudes 3 and 2.
##############################
time     <- 400                           # measuring time interval (in time units)
acq.freq <- 5                           # data acquisition frequency (Hz), how many measurements
xt <- seq(0,time-1/acq.freq,1/acq.freq)   # vector of sampling time-points (s) 
f.0 <- 1/time # fundamental frequency
x <- 3*sin(xt*0.2)-2*cos(xt*0.5)
n <- length(x)
plot(xt,x, type="l")
x.harm <- harmonic(xt, x, acq.freq=acq.freq, N=5, alpha=0.05, detrend=TRUE, test="holm")

# This perfectly reconstructs the original time series from all harmonics:
xk <- fft(x)/n
windows()
plot.frequency.spectrum(xk, acq.freq=acq.freq)
x.n <- get.trajectory(xk,xt,acq.freq)
windows()
plot(xt, Re(x.n), type="l"); abline(h=0,lty=3)
points(xt, x, col="red",type="l", lty=2, lwd=2) # compare with original

# now add a trend:
time     <- 400                           # measuring time interval (in time units)
acq.freq <- 5                           # data acquisition frequency (Hz), how many measurements
xt <- seq(0,time-1/acq.freq,1/acq.freq)   # vector of sampling time-points (s) 
f.0 <- 1/time # fundamental frequency
x <- 0.2*xt+3*sin(xt*0.2)-2*cos(xt*0.5)
n <- length(x)
plot(xt,x, type="l")
x.harm <- harmonic(xt, x, acq.freq=acq.freq, N=5, alpha=0.05, detrend=TRUE, test="holm")

# now add an offset:
x <- -8+0.2*xt+3*sin(xt*0.2)-2*cos(xt*0.5)
plot(xt, x, type="l")
x.harm <- harmonic(xt, x, acq.freq=acq.freq, N=5, alpha=0.05, detrend=TRUE, test="holm")

# now add random noise:
x <- -8+0.2*xt+3*sin(xt*0.2)-2*cos(xt*0.5)+1.5*rnorm(xt)
plot(xt, x, type="l")
x.harm <- harmonic(xt, x, acq.freq=acq.freq, N=5, alpha=0.05, detrend=TRUE, test="holm")

####################################
# test the ANOVA table with the published Hartley data:
####################################
x <- c(12, 13, 13, 18, 23, 29, 31, 33, 29, 18, 19, 8, 1, 8, 14, 18, 26, 33, 31, 32,
  27, 20, 15, 9)
time     <- 24                          # measuring time interval (in time units)
acq.freq <- 1                           # data acquisition frequency (Hz), how many measurements
xt <- seq(0,time-1/acq.freq,1/acq.freq)  # vector of sampling time-points (s) 
n <- length(x)
plot(xt, x, type="l")
x.harm <- harmonic(xt, x, acq.freq=acq.freq, N=6, alpha=0.05, detrend=TRUE, which="first", test="bonferroni")

####################################
# And now process data from the Sutton Bonington Met Station:
####################################

# file name and number of header lines in the file to skip
sb <- read.table("suttonboningtondata.txt", sep=" ", skip = 7)
# column names of the file
colnames(sb) <- c("DELETE", "yr", "m", "tmax_C", "tmin_C", "af_days", "rain_mm", "sun_hrs")
sb <- sb[,-1] # remove the first data column (blank space in the text file)
tmax <- sb$tmax_C # extract maximum monthly temperature measurements
time     <- length(tmax)                # measuring time interval (in months)
acq.freq <- 1                           # data acquisition frequency (Hz), how many measurements per month
ts <- seq(0,time-1/acq.freq,1/acq.freq)  # vector of sampling time-points in months
n <- length(tmax)
f.0 <- 1/time
plot(tmax, type="l", xlab="t", ylab="tmax") # show the data

# with Bonferroni adjusted alpha error
tmax.harm <- harmonic(ts, tmax, N=100, alpha=0.05, detrend=TRUE, test="bonferroni")

# with Holm adjusted alpha error
tmax.harm <- harmonic(ts, tmax, N=100, alpha=0.05, detrend=TRUE, test="holm")

sum(tmax.harm$Amp)

####################################
# Process Chlorophyll-a data from Lake Balaton by Palmer et al. (2014):
####################################

# close any open graphics windows
for (i in 1:30) dev.off()

file <- "Log_Balaton_Chla_harmonics.txt"
sink(file)
cat("Balaton chlorophyll-a analysis from MERIS\n")
cat("Log file: ", file, "\n")
file <- "Balaton_Chla_HB-harmonics.txt"
cat("Data file: ", file, "\n")
chla <- read.table(file, sep="\t", header=TRUE)
time <- length(chla[,1])/36              # measuring time interval (in years), the original data is in dekads
acq.freq <- 36                           # data acquisition frequency (Hz), how many measurements per year
ts <- seq(0,time-1/acq.freq, 1/acq.freq) # vector of sampling time-points in months
n <- length(chla[,1])
f.0 <- 1/time
chlafilled <- chla
for (site in 1:length(chla[1,])) {
  cat("\nSite: ", names(chla)[site], "\n")
  cat("\nGap-filling the data:\n")
  cat("Number of NAs: ", length(which(is.na(chla[, site]))), "\n")
  windows()
  plot(chla[,site], type="l", xlab="t", ylab="chl-a", lwd=1, main=names(chla)[site]) 
  for (t in 1:length(chla[,1])) { # for each time step
    if (is.na(chla[t,site])) { # is the value missing?
      tlow <- max(which(!is.na(chla[1:(t-1), site]))) # time index of nearest earlier data point
      thi  <- t+min(which(!is.na(chla[(t+1):n, site]))) # time index of nearest later data point
      xlow <- chla[tlow, site] # value of nearest earlier data point
      xhi  <- chla[thi,  site] # value of nearest later data point
      tt <- c(tlow, thi)
      xx <- c(xlow, xhi)
      xlm <- lm(xx~tt)
      chlafilled[t,site] <- xlm$coef[1]+xlm$coef[2]*t # fill in
      lines(c(t-1, t), c(chlafilled[t-1,site],chlafilled[t,site]), type="l",col=2, lwd=2)
      if (!is.na(chlafilled[min(t+1,n),site])) lines(c(t, min(t+1,n)), c(chlafilled[t,site],chlafilled[min(t+1,n),site]), type="l",col=2, lwd=2)
      } # endif
    } # end of for t
  cat("\nCalculating harmonic terms for site .", names(chla)[site],"\n")
  chla.harm <- harmonic(ts, chlafilled[,site], N=40, alpha=0.05, detrend=TRUE, which="strongest", test="holm")
  if (site==1) b1.harm <- chla.harm # store the results
  if (site==2) b2.harm <- chla.harm # note that the site numbers are column numbers
  if (site==3) b3.harm <- chla.harm # they do not necessarily match the data column names
  if (site==4) b4.harm <- chla.harm
} # end for site

#############################
# Now, make publication quality graphs, equal axes etc.
# resolution of jpeg in dpi
res <- 72
# width
width <- 24*res
# height
height <- 14*res
# point size for text
ps <- 30
# colour of second data plots
col <- 2
#############################
for (site in 1:length(chla[1,])) {
  if (site==1) chla.harm <- b1.harm
  if (site==2) chla.harm <- b2.harm
  if (site==3) chla.harm <- b3.harm
  if (site==4) chla.harm <- b4.harm

  cat("Plotting the gap filled data.\n")
  jpeg(filename = paste(names(chla)[site], "_data.jpg", sep=""), width = width, height = height, 
      units = "px", pointsize = ps, quality = 100, bg = "white", res = res, restoreConsole = TRUE)
  plot(chla.harm$xt, chlafilled[,site], type="l", col=col, lwd=3, lty=2, main=names(chla)[site],
    ylim=c(0,50))
  lines(chla.harm$xt, chla[,site], type="l", xlab="t", ylab="chl-a", lwd=3, col=1)
  dev.off()

  cat("Plotting the periodogram.\n")
  jpeg(filename = paste(names(chla)[site], "_periodogram.jpg", sep=""), width = height, height = height, 
      units = "px", pointsize = ps, quality = 100, bg = "white", res = res, restoreConsole = TRUE)
  plot.frequency.spectrum(chla.harm$frequency.spectrum, acq.freq=acq.freq)
  dev.off()

  cat("Plotting the individual harmonic terms.\n")
  jpeg(filename = paste(names(chla)[site], "_harmonics.jpg", sep=""), width = width, height = height, 
      units = "px", pointsize = ps, quality = 100, bg = "white", res = res, restoreConsole = TRUE)
  plot(chla.harm$xt, chlafilled[,site]-chla.harm$lm[[1]][1]-chla.harm$lm[[1]][2]*chla.harm$xt, 
      type="l", lwd=2, ylim=c(-20,40))
  abline(h=0,lty=2)
  N <- chla.harm$Nsig
  for (i in 1:N) {
    plot.harmonic(chla.harm$frequency.spectrum, chla.harm$ndx[i], chla.harm$xt, acq.freq, col=i+1)
    }
  dev.off()

  cat("Plotting the composite harmonic model with trend added back on against the original data.\n")
  jpeg(filename = paste(names(chla)[site], "_predicted.jpg", sep=""), width = width, height = height, 
      units = "px", pointsize = ps, quality = 100, bg = "white", res = res, restoreConsole = TRUE)
  plot(chla.harm$xt, chlafilled[,site], type="l", col=1, lwd=2, ylim=c(0,50))
  lines(chla.harm$xt, Re(chla.harm$fitted.values), type="l", xlab="t", ylab="chl-a", lwd=2, col=col)
  abline(chla.harm$lm, lty = 2)
  dev.off()

  cat("Plotting the residuals.\n")
  jpeg(filename = paste(names(chla)[site], "_residuals.jpg", sep=""), width = width, height = height, 
      units = "px", pointsize = ps, quality = 100, bg = "white", res = res, restoreConsole = TRUE)
  plot(chla.harm$xt, chla.harm$residuals, col=1, lwd=2, pch="+", ylim=c(-20,20))
  dev.off()
}

sink()


####################################
# END
####################################
