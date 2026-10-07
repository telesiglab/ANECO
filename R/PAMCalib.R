#' Sound pressure level calibration
#'
#' It is an internal function of the \code{\link[ANECO]{ANECO}} package used in the function \code{\link[ANECO]{acoustic_indices}}. This function transform the intensity values of sound files into calibrated dB values according to technical parameters of the recorder device
#'

#' @param tempWav Character. Full path to the sound file
#' @param ifile Character. Sound file name
#' @param atype Ver Merchant et al. (2015)
#' @param plottype Ver Merchant et al. (2015)
#' @param envi Ver Merchant et al. (2015)
#' @param calib Ver Merchant et al. (2015)
#' @param ctype Ver Merchant et al. (2015)
#' @param Si Ver Merchant et al. (2015)
#' @param Mh Ver Merchant et al. (2015)
#' @param G Ver Merchant et al. (2015)
#' @param vADC Ver Merchant et al. (2015)
#' @param r Ver Merchant et al. (2015)
#' @param N Ver Merchant et al. (2015)
#' @param winname Ver Merchant et al. (2015)
#' @param lcut Ver Merchant et al. (2015)
#' @param hcut Ver Merchant et al. (2015)
#' @param timestring Ver Merchant et al. (2015)
#' @param outdir Ver Merchant et al. (2015)
#' @param outwrite Ver Merchant et al. (2015)
#' @param disppar Ver Merchant et al. (2015)
#' @param welch Ver Merchant et al. (2015)
#' @param chunksize Ver Merchant et al. (2015)
#' @param linlog Ver Merchant et al. (2015)
#'
#' @return A matrix with the values calibrated in decibels, each row is a temporal value of duration TimeRes and each column is a frequency value
#'
#' @references \itemize{
#' \item{Merchant, N. D., Fristrup, K. M., Johnson, M. P., Tyack, P. L., Witt, M. J., Blondel, P., & Parks, S. E. (2015). Measuring acoustic habitats. Methods in Ecology and Evolution, 6(3), 257-265.}{}}
#'
#' @export
#'
PAMCalib <- function(...,atype='PSD',plottype='Both',envi='Air',calib=0,ctype = 'MF',Si = NULL,Mh = NULL,G = NULL,vADC = NULL,r=50,N=Fs,winname='Hann',lcut=Fs/N,hcut=Fs/2,timestring="",outdir=dirname(fullfile),outwrite=0,disppar=0,welch="",chunksize="",linlog = "Log", tempWav = tempWav, ifile = ifile){

  aid <- 0										#reset metadata code
  if (calib == 0) {aid <- aid + 20				#add calibration element to metadata code
  } else {aid <- aid + 10}
  if (timestring != ""){aid <- aid + 1000} else {aid <- aid + 2000}	#add time stamp element to metadata code


  ## Select input file and get file info-----------------------------------------

  fullfile <- tempWav
  ifile <- ifile
  fIN <- tuneR::readWave(fullfile,header = TRUE)			#read file header
  Fs <- fIN[[1]]									#sampling frequency
  Nbit <- fIN[[3]]								#bit depth
  xl <- fIN[[4]]									#length of file in samples
  xlglo <- xl										#back-up file length


  ## Read time stamp data if provided-----------------------------------------

  if (timestring != "")
  {
    tstamp <- as.POSIXct(strptime(ifile,timestring) ,origin="1970-01-01")
    if(disppar == 1)
    {
      cat('Time stamp start time: ',format(tstamp),'\n')
    }
  }
  if (timestring == "") tstamp <- NULL	#compute time stamp in R POSIXct format



  r<-r/100

  ## Read input file------------------------------------------------------------------

  if (chunksize == ""){nchunks = 1				#if loading whole file, nchunks = 1
  } else if (chunksize != "") {				#if loading file in stages
    nchunks <- ceiling(xl/(Fs*as.numeric(chunksize)))	#number of chunks of length chunksize in file
  }

  for (q in 1:nchunks){							#loop through file chunks
    if (nchunks == 1){							#if loading whole file at once
      t1=proc.time()									#start timer
      xbit <- tuneR::readWave(fullfile)						#read file
      xbit <- xbit@left/(2^(Nbit-1))					#convert to full scale (+/- 1) via bit depth
    } else if (nchunks > 1) {						#if loading file in stages
      if (q == nchunks){							#load qth chunk
        xbit <- tuneR::readWave(fullfile,from=((q-1)*as.numeric(chunksize)*Fs+1),to=xlglo,units="samples")
        xbit <- xbit@left/(2^(Nbit-1))			#convert to full scale (+/- 1) via bit depth
        xl <- length(xbit)						#final chunk length
      } else {
        xbit <- tuneR::readWave(fullfile,from=((q-1)*as.numeric(chunksize)*Fs+1),to=(q*as.numeric(chunksize)*Fs),units="samples")
        xbit <- xbit@left/(2^(Nbit-1))			#convert to full scale (+/- 1) via bit depth
        xl <- length(xbit)						#chunk length
      }
    }

    if (envi == 'Air'){pref<-20; aid <- aid+100}	#set reference pressure depending for in-air or underwater
    if (envi == 'Wat'){pref<-1; aid <- aid+200}


    ## Compute system sensitivity if provided-----------------------------------------------

    if (calib == 1){
      if (ctype == 'EE') {		#'EE' = end-to-end calibration
        S <- Si}
      if (ctype == 'RC') {		#'RC' = recorder calibration with separate transducer sensitivity defined by Mh
        S <- Si + Mh}
      if (ctype == 'TS') {		#'TS' = manufacturer's specifications
        S <- Mh + G + 20*log10(1/vADC);		#EQUATION 4
      }

      if (disppar == 1){cat('System sensitivity correction factor, S = ',sprintf('%.01f',S),' dB\n')}
    } else {S <- 0}				#if not calibrated (calib == 0), S is zero

    ## Compute waveform if selected----------------------------------------------------

    if (atype == 'Waveform') {
      if (calib == 1){
        a <- xbit/(10^(S/20)) 				#EQUATION 21
      } else {a <- xbit/(max(xbit))}
      t <- seq(1/Fs,length(a)/Fs,1/Fs)	#time vector
      tana = proc.time()
    }

    ## Compute DFT-based metrics if selected----------------------------------------

    if (atype != 'Waveform') {
      if (nchunks == 1){
        if (atype == 'PowerSpec'){cat('Computing power spectrum...')}
        if (atype == 'TOL'){cat('Computing TOLs by DFT method...')}
        if (atype == 'Broadband'){cat('Computing broadband level...')}
        tana = proc.time()}

      # Divide signal into data segments (corresponds to EQUATION 5)
      N = round(N)						#ensure N is an integer
      nsam = ceiling((xl)-r*N)/((1-r)*N)	#number of segments of length N with overlap r
      xgrid <- matrix(nrow = N,ncol = nsam)	#initialise grid of segmented data for analysis
      for (i in 1:nsam) {					#make grid of segmented data for analysis
        loind <- (i-1)*(1-r)*N+1
        hiind <- (i-1)*(1-r)*N+N
        xgrid[,i] = xbit[loind:hiind]
      }

      M <- length(xgrid[1,])

      # Apply window function (corresponds to EQUATION 6)
      if (winname == 'Rectangular') {			#rectangular (Dirichlet) window
        w <- matrix(1,1,N)
        alpha <- 1 }					#scaling factor
      if (winname == 'Hann') {			#Hann window
        w <- (0.5 - 0.5*cos(2*pi*(1:N)/N))
        alpha <- 0.5 }					#scaling factor
      if (winname == 'Hamming') {			#Hamming window
        w <- (0.54 - 0.46*cos(2*pi*(1:N)/N))
        alpha <- 0.54 }					#scaling factor
      if (winname == 'Blackman') {		#Blackman window
        w <- (0.42 - 0.5*cos(2*pi*(1:N)/N) + 0.08*cos(4*pi*(1:N)/N))
        alpha <- 0.42 }					#scaling factor

      xgrid <- xgrid*w/alpha 				#apply window function


      #Compute DFT (corresponds to EQUATION 7)
      X <- abs(mvfft(xgrid))

      #Compute power spectrum (EQUATION 8)
      P <- (X/N)^2

      #Compute single-side power spectrum (EQUATION 9)
      # FFT row 1 is DC; subsequent rows are k * Fs / N.
      bins <- 0:floor(N/2)
      Pss <- P[bins + 1, , drop = FALSE]
      sided_scale <- rep(2, length(bins))
      sided_scale[1] <- 1
      if (N %% 2 == 0) sided_scale[length(sided_scale)] <- 1
      Pss <- Pss * sided_scale

      # Keep the same frequency bins in both the spectrum and its labels.
      f <- bins * Fs/N
      selected <- which(f >= lcut & f <= hcut)
      if (!length(selected)) {
        stop("No FFT frequencies fall within lcut and hcut.")
      }
      f <- f[selected]
      Pss <- Pss[selected, , drop = FALSE]
      nf <- length(f)


      #Compute PSD in dB if selected
      if (atype == 'PSD') {
        B <- (1/N)*(sum((w/alpha)^2))		#noise power bandwidth (EQUATION 12)
        delf <- Fs/N;						#frequency bin width
        a <- 10*log10((1/(delf*B))*Pss/(pref^2))-S
      }									#PSD (EQUATION 11)

      #Compute power spectrum in dB if selected
      if (atype == 'PowerSpec') {
        a <- 10*log10(Pss/(pref^2))-S
      }									#EQUATION 10

      #Compute broadband level if selected
      if (atype == 'Broadband') {
        a <- 10*log10(colSums(Pss)/(pref^2))-S
      }									#EQUATION 17

      #Compute 1/3-octave band levels if selected
      if (atype == 'TOL') {
        if (lcut <25){					#limit TOL analysis to > 25 Hz
          lcut <- 25}

        #Generate 1/3-octave freqeuncies
        lobandf <- floor(log10(lcut))	#lowest power of 10 for TOL computation
        hibandf <- ceiling(log10(hcut))	#highest power of 10 for TOL computation
        nband <- 10*(hibandf-lobandf)+1	#number of 1/3-octave bands
        fc <- matrix(0,nband)			#initialise 1/3-octave frequency vector
        fc[1] <- 10^lobandf;

        #Calculate centre frequencies (corresponds to EQUATION 13)
        for (i in 2:nband) {
          fc[i] <- fc[i-1]*10^0.1}
        fc <- fc[which(fc >= lcut)[1]:max(which(fc <= hcut))]

        nfc <- length(fc)				#number of 1/3-octave bands

        #Calculate boundary frequencies of each band (EQUATIONS 14-15)
        fb <- fc*10^-0.05				#lower bounds of 1/3-octave bands
        fb[nfc+1] <- fc[nfc]*10^0.05	#upper bound of highest band
        if (max(fb) > hcut) {			#if upper bound exceeds highest
          nfc <- nfc-1				# frequency in DFT, remove
          fc <- fc[1:nfc]}

        #Calculate TOLs (corresponds to EQUATION 16)
        P13 <- matrix(nrow = M,ncol = nfc)
        for (i in 1:nfc) {
          fli <- which(f >= fb[i])[1]
          fui <- max(which(f < fb[i+1]))
          for (k in 1:M) {
            fcl <- sum(Pss[fli:fui,k])
            P13[k,i] <- fcl
          }
        }
        a <- t(10*log10(P13/(pref^2)))-S
      }

      # Compute time vector
      tint <- (1-r)*N/Fs
      ttot <- M*tint-tint
      t <- seq(0,ttot,tint)
    }

    if (nchunks>1){								#if loading in stages, concatenate analyses in each loop iteration
      if (q == 1){
        newa <- a
        newt <- t
        cat('Chunk size:',chunksize,'s\nAnalysing in',nchunks,'chunks. Analysing chunk 1')
      } else if (q > 1){
        dima <- dim(a)
        newa <- cbind(newa,a)
        if (timestring != ""){
          newt <- c(newt,t)
        } else {
          newt <- c(newt,t+(q-1)*as.numeric(chunksize))
        }
        cat(' ',q)
      }
    }
  }

  if (nchunks>1){a <- newa					#reassign output array
  t <- newt
  cat('\n')}

  # If not calibrated, scale relative dB to zero
  if (calib == 0 & atype != 'Waveform') {a <- a-max(a)}

  if (!is.null(tstamp)){t <- t+tstamp}

  #tdiff <- max(t)-min(t)					#define time format for x-axis of time plot
  #if (tdiff < 10){
  #  tform <- "%H:%M:%S:%OS3"}
  #else if (tdiff > 10 & tdiff < 86400){
  #  tform <- "%H:%M:%S"}
  #else if (tdiff > 86400 & tdiff < 86400*7){
  #  tform <- "%H:%M \n %d %b"}
  #else if (tdiff > 86400*7){tform <- "%d %b %y"}
  #}

  ## Construct output array------------------------------------------------

  if (atype == 'PSD' | atype == 'PowerSpec') {
    A <- cbind(t,t(a))
    A <- rbind(c(0,f),A)
    if (atype == 'PSD'){aid <- aid + 1}
    if (atype == 'PowerSpec'){aid <- aid + 2}
    A[1,1] <- aid
  }

  if (atype == 'TOL') {
    A <- cbind(t,t(a))
    A <- rbind(c(0,fc),A)
    aid <- aid + 3
    A[1,1] <- aid
    f <- fc
  }

  if (atype == 'Broadband') {
    A <- t(rbind(t,a))
    aid <- aid + 4
    A[1,1] <- aid
  }

  if (atype == 'Waveform') {
    A <- t(rbind(t,a))		#define output array
    A <- rbind(c(0,0),A)	#add zero top row for metadata
    aid <- aid + 5			#add index to metadata for Waveform
    A[1,1] <- aid			#encode output array with metadata
  }

  ## Reduce time resolution if selected------------------------------------------

  dimA <- dim(A)				#dimensions of output array

  if (welch != "" && atype != "Waveform"){
    lout <- ceiling(dimA[1]/welch)+1	#length of new, Welch-averaged, output array
    AWelch <- matrix(, nrow = lout, ncol = dimA[2])	#initialise Welch array
    AWelch[1,] <- A[1,]					#assign frequency vector
    tint <- A[3,1] - A[2,1]				#time interval between segments
    cat('Welch factor =',welch,'x\nNew time resolution =',welch,'(Welch factor) x',N/Fs,'s (time segment length) x',r*100,'% (overlap) =',welch*tint,'s\n')
    if (lout == 2){						#if Welch factor too large for number of time data points, abort averaging
      stop(paste('Welch factor is greater than half the number of samples. Set welch <=',dimA[1]/2,', or reduce N.',sep=""),call.=FALSE)
    }	else {
      for (i in 2:lout) {						#loop through Welch segments for averaging
        stt <- A[2,1] + (i-2)*tint*welch	#start time
        ett <- stt + welch*tint				#end time
        stiv <- which(A[2:dimA[1]]>=stt)
        sti <- min(stiv)+1					#start index
        etiv <- which(A[2:dimA[1]]<ett)
        eti <- max(etiv)+1					#end index
        nowA <- 10^(A[sti:eti,]/10)			#take RMS level of segments as per Welch method
        AWelch[i,] <- 10*log10(rowMeans(t(nowA)))	#convert to dB
        AWelch[i,1] <- stt+tint*welch/2		#assign time index
      }
    }
    A <- AWelch										#reassign output array
    dimA <- dim(A)
    t <- A[2:dimA[1],1]
    f <- A[1,2:dimA[2]]
    a <- t(A[2:dimA[1],2:dimA[2]])
  }
  dimA <- dim(A)


  ## Write output array to CSV file if selected----------------------

  A <- data.matrix(A, rownames.force = NA)
  return(A)
}
