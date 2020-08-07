AddMaxTime <- function(wake.df, max.time){
  # add 78 hours as last timepoint
  dat.eeg <- wake.df
  wake.df <- subset(dat.eeg, time.shift >= 0 & time.shift <= max.time)
  wake.df.dup <- wake.df[wake.df$time.shift == max(wake.df$time.shift), ]
  wake.df.dup$time.shift <- max.time
  wake.df <- rbind(wake.df, wake.df.dup)
  return(wake.df)
  
}
GetWakeCollapsed <- function(inf, tstep=4/3600, max.time, filt.time, return.wake.df.only = FALSE){
  # if (missing(inf)){
  #   inf <- "Robjs/eeg_data_merged_72_and_78_hr_mice/wake.df.method.mode.Robj"
  # }
  # load(inf, verbose = T)
  # inf is dummy
  data("wake.df.method.mode")
  wake.df <- AddMaxTime(wake.df, max.time)
  if (return.wake.df.only){
    return(wake.df)
  } else {
    wake.collapsed <- CollapseWake(wake.df, tstep, filter.time = filt.time)
    return(wake.collapsed)
  }
}

GetSubsampledEeg <- function(inf, max.time, downsamp = 0.005){
  # get wake.df
  # inf not useful, force load from /data
  # if (missing(inf)){
  #   inf <- "Robjs/eeg_data_merged_72_and_78_hr_mice/wake.df.method.mode.Robj"
  # }
  # load(inf, verbose = T)
  data("wake.df.method.mode")
  wake.df <- AddMaxTime(wake.df, max.time)
  dat.eeg.plot <- wake.df[seq(1, nrow(wake.df), length.out = downsamp * nrow(wake.df)), ]
  return(dat.eeg.plot)
}




DownSample <- function(dat, downsamp){
  return(dat[seq(1, nrow(dat), length.out = downsamp * nrow(dat)), ])
}

Mode <- function(x) {
  # http://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode 
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

TakeAvgWakeStatus <- function(dat, method = "mean"){
  # take mean and mode of wake status
  # method can be mean or mode
  if (method == "mean"){
    wake <- mean(dat$wake)
  } else if (method == "mode"){
    wake <- Mode(dat$wake)
  } else {
    warning("Method must be mean or mode")
  }
  N <- length(dat$wake)
  return(data.frame(wake = wake, N = N))
}

DownsampleEeg <- function(dat.eeg, downsamp = 0.05){
  # for plotting purposes 
  dat.eeg.plot <- dat.eeg[seq(1, nrow(dat.eeg), length.out = downsamp * nrow(dat.eeg)), ]
  return(dat.eeg.plot)
}
EegToWake <- function(e, rem.is.wake = FALSE){
  # eeg status as per Franken email 2016-04-20 sleep-wake model - meeting
  # 'w' '4' '1' -> wake
  # 'n' '2' '5' -> NREM
  # 'r', '3', '6' -> REM
  # convert signals to either wake or not wake
  # 2016-09-04: set nrem same as wake by default (Franken)
  if (rem.is.wake){
    wake.sig <- c("w", 4, 1, "r", 3, 6)
  } else {
    wake.sig <- c("w", 4, 1)
  }
  if (e %in% wake.sig){
    return(1)
  } else {
    return(0)
  }
}

SmoothWakefulness5Minutes <- function(dat.eeg){
  # get minutes of wakefulness of last 5 minutes
  smooth.window <- 75  # 5 minute intervals (4 sec/interval)
  weights <- rep(5 / smooth.window, smooth.window)  # minutes awake per 5 minute inte
  dat.eeg$w.smooth <- stats::filter(dat.eeg$wake, filter = weights, sides = 1)
  return(dat.eeg)
}

CollapseWake <- function(wake.df, tstep, filter.time, jtol=1e-6){
  # collapse wake.df to include only times when you change wake from 0 to 1 OR 
  # if time == filter.time (expected to be whole integers or rational numbers?) 
  # We do the tstep count in order to be able to quickly subset for times of interest 
  # e.g, times at which samples are taken
  wake.collapsed <- data.frame()
  i.init <- wake.df$wake[[1]]
  time <- 0  # init
  i <- 0
  for (w in wake.df$wake[-1]){
    i <- i + 1
    time.cum <- i * tstep
    time <- time + tstep
    # if ((time.cum %in% filter.time) | (w != i.init)){
    if (any(abs(time.cum - filter.time) <= jtol)){
      # adjust time.cum to be exactly same as filter.time
      t.i <- which(abs(time.cum - filter.time) <= jtol)
      if (length(t.i) != 1){
        stop("Length of t.i must be 1")
      }
      time.cum <- filter.time[t.i]
    }
    if (any(abs(time.cum - filter.time) <= jtol) | w != i.init){
      # record and reupdate
      wake.collapsed <- rbind(wake.collapsed, data.frame(wake = i.init, time.duration = time, time.cum = time.cum))
      i.init <- w
      time <- 0
    }
    # if (abs(time.cum - filter.time) <= jtol | (w != i.init)){
    #   # record and reupdate
    #   wake.collapsed <- rbind(wake.collapsed, data.frame(wake = i.init, time.duration = time, time.cum = time.cum))
    #   i.init <- w
    #   time <- 0
    # }
  }
  return(wake.collapsed)
}

ShiftEeg <- function(dat.eeg){
  dat.eeg$time.shift <- dat.eeg$time - 24
  jxlim <- c(-24, 78)
  dat.eeg <- subset(dat.eeg, time.shift >= jxlim[1] & time.shift <= jxlim[2])
  # add 72 as -24
  row.to.add <- subset(dat.eeg, time.shift == -24); row.to.add$time.shift <- 72
  dat.eeg <- rbind(dat.eeg, row.to.add)
  return(dat.eeg)
}

ShiftCollapseEeg <- function(dat.eeg, filt.time, tstep = 4/3600){
  # collapse eeg into only certain timepoints, defined by filt.time
  # otherwise you get only time duration of when you are wake or sleep
  if (missing(filt.time)){
    filt.time <- c(0, 3, 6, 12, 18, 27, 30, 36, 42, 48, 54, 72)  # sleep deprivation timepoints
  }
  dat.eeg <- ShiftEeg(dat.eeg)
  
  wake.df <- subset(dat.eeg, time.shift >= 0 & time.shift <= 72)
  # add 72 hours as 71.99889
  wake.df.dup <- wake.df[wake.df$time.shift == max(wake.df$time.shift), ]
  wake.df.dup$time.shift <- 72
  wake.df <- rbind(wake.df, wake.df.dup)
  
  wake.collapsed <- CollapseWake(wake.df, tstep, filter.time = filt.time)
  return(wake.collapsed)
}
