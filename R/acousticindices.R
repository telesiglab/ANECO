#' Compute acoustic indices per minute with optional SPL calibration
#'
#' Computes a configurable set of acoustic indices for 1-minute segments from all \code{.wav} audio files found in \code{dir}. Indices can be computed from calibrated spectra (absolute SPL; via \code{\link[ANECO]{PAMCalib}}) and, optionally, a subset of commonly used soundscape indices from \pkg{seewave} and \pkg{soundecology} (computed without the SPL calibration step, using package defaults).
#'
#' The output is a minute-level table including identification fields parsed from file names (using \code{prefix.format}) and the requested indices. The function supports parallelization either across audio files (\code{parallel = "files"}) or across minutes within each file (\code{parallel = "minutes"}), and can optionally write the results to a CSV file in \code{dir}.
#'
#' @param dir Character. Path to a folder containing \code{.wav} audio files.
#'
#' @param calibparam Calibration parameters used to convert amplitude measurements to absolute SPL. Either:
#'   \itemize{
#'     \item a numeric vector \code{c(micsens, gain, vADC)}; or
#'     \item a character scalar selecting a recorder preset: \code{"SM4"},
#'       \code{"SMmini"}, or \code{"SM2+"}.
#'   }
#'   Where \code{micsens} is microphone sensitivity (dB re 1 V/Pa), \code{gain} is recorder gain (dB), and \code{vADC} is the ADC zero-to-peak voltage (V). Sensitivity is converted internally to dB re 1 V/\eqn{\mu}Pa for SPL calibration. For background and recommended calibration practice, see Merchant et al.
#'
#' @param sel.ind Character vector specifying which indices to compute. Use \code{"all"} (default) to compute all indices implemented in this function, or provide a character vector of index names (e.g., \code{c("ACI","ADI","msldB_bio")}). Available index codes can be inspected via \code{aneco_list_indices()} or using \code{ANECO::indices_catalog}. See \strong{Details}.
#'
#' @param exclude Character vector of index names to exclude \emph{only when} \code{sel.ind = "all"}. Use \code{"none"} (default) to exclude nothing.
#'
#' @param noise.ind Logical. If \code{TRUE}, ensures that the subset of indices used by  ANECO’s rain/noise classification workflow is included (e.g., \code{msldB_bio}, \code{msldB_low}, \code{mdBGL}, spectral statistics, \code{dfreq}, and band energies). Default is \code{TRUE}.
#'
#' @param bioband Numeric length-2 vector giving the lower and upper frequency bounds (Hz) used to define the biophonic band for calibrated indices. Default is \code{c(1000, 11000)}.
#'
#' @param noiseband Numeric length-2 vector giving the lower and upper frequency bounds (Hz) used to define the low-frequency band for calibrated indices. Default is \code{c(0, 1000)}.
#'
#' @param channel Character. Channel to analyze when converting to mono.
#'   Must be \code{"left"} or \code{"right"}. Default is \code{"left"}.
#'
#' @param prefix.format Character string describing how site/recorder codes are embedded in the file name prefix. Uses \code{"S"} for site characters and \code{"R"} for recorder characters. For example, for files such as \code{BAJA08_20210523_050000.wav}, where \code{BAJA} is the site and \code{08} is the recorder ID, a suitable format is \code{"SSSSRR"}. Default is \code{"SSRR"}.
#'
#' @param wl Numeric. FFT window length (in samples) used during spectral estimation for calibrated indices. Default is \code{512}.
#'
#' @param ncores Integer. Number of workers used for parallel computation. If \code{NULL} (default), uses physical cores minus one (minimum 1). Availability can be checked with \code{parallel::detectCores(all.tests = TRUE, logical = FALSE)}.
#'
#' @param save.file Logical. If \code{TRUE}, writes the resulting table to a CSV file in \code{dir}. The file name is built from the first \code{Code} value and a timestamp. Default is \code{TRUE}.
#'
#' @param add.prefix Character. If a file name lacks a prefix (i.e., only contains the date/time fields), this value is prepended to allow parsing of \code{Site} and \code{RecorderID}. Default is \code{"XXXX"}.
#'
#' @param parallel Character string controlling the parallelization level:
#'   \itemize{
#'     \item \code{"files"}: process different audio files in parallel (each file processed sequentially by minute).
#'     \item \code{"minutes"}: process minutes within each file in parallel (files are processed sequentially).
#'   }
#'   Default is \code{"files"}.
#'
#' @details
#' \strong{Spectral integration.} Calibrated power spectral density is converted to power per FFT bin using the bin width \code{sample_rate / round(wl)} before calculating sound levels and calibrated indices.
#'
#' \strong{Mixed sample rates.} Results are combined by column name. Energy bands unavailable at a file's sample rate are filled with \code{NA}, rather than zero.
#'
#' \strong{Minute segmentation.} Each \code{.wav} file is read in 60-second blocks. The final segment is included only if its duration is at least 30 seconds.
#'
#' \strong{Index selection.} \code{sel.ind} can be \code{"all"} or a vector of specific index names. Available index codes can be inspected via \code{aneco_list_indices()} or using \code{ANECO::indices_catalog}.  When \code{sel.ind = "all"}, \code{exclude} removes indices from the full set; when \code{exclude = "none"}, nothing is removed. When \code{noise.ind = TRUE}, a predefined subset required for rain/noise classification is added to the selection if missing.
#'
#' \strong{Identification fields.} The output includes:
#' \code{Name} (file name), \code{Code} (prefix), \code{Site}, \code{RecorderID}, \code{Date}, \code{Year}, \code{Month}, \code{Day}, \code{Hour}, \code{Minute}, and \code{Duration} (seconds).
#'
#' @return A data.frame with one row per analyzed minute. Includes identification columns
#'   (\code{Name}, \code{Code}, \code{Site}, \code{RecorderID}, \code{Date}, \code{Year}, \code{Month}, \code{Day}, \code{Hour}, \code{Minute}, \code{Duration}) and the requested index columns. Files/minutes that fail processing are omitted; the function may emit a warning listing file names that produced errors.
#'   If no valid minutes remain, returns an empty data.frame and does not write a CSV file.
#'
#' @references
#' \itemize{
#'   \item Buxton, R. T., et al. (2018). Efficacy of extracting indices from large-scale acoustic recordings to monitor biodiversity. \emph{Conservation Biology}.
#'   \item Merchant, N. D., et al. (2015). Measuring acoustic habitats. \emph{Methods in Ecology and Evolution}, 6(3), 257–265.
#'   \item Sueur, J., et al. (2016). \pkg{seewave}: Sound analysis and synthesis.
#'   \item Villanueva-Rivera, L. J. & Pijanowski, B. C. (2016). \pkg{soundecology}: Soundscape Ecology.
#' }
#'
#' @seealso \code{\link[ANECO]{soundclassifier_fit}}, \code{\link[ANECO]{soundclassifier_predict}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Directory containing .wav files
#' path_folder_audios <- "C:/path/to/your/audio/files"
#'
#' # Calibration parameters (Wildlife acoustics SM4 recorder preset)
#' calibration_parameters <- "SM4"
#'
#' # Enable progress reporting (optional)
#' progressr::handlers(global = TRUE)
#'
#' res <- acoustic_indices(
#'   dir = path_folder_audios,
#'   calibparam = calibration_parameters,
#'   sel.ind = c("ACI", "ADI", "msldB_bio"),
#'   exclude = "none",
#'   noise.ind = TRUE,
#'   parallel = "files"
#' )
#' }

acoustic_indices <- function(dir = NULL, calibparam = "SM4", sel.ind = "all", exclude = "none", noise.ind = TRUE, bioband = c(1000, 11000), noiseband = c(0, 1000), channel = "left", prefix.format = "SSRR", wl = 512, ncores = NULL, save.file = TRUE, add.prefix = "XXXX", parallel = "files"){

#----------------------------------
# Checking requierements
# ---------------------------------

  if(is.null(dir)){
    stop("dir argument not defined. You must specify the path where the audio data is located")
  }

  if(length(calibparam) == 3) {
    micsens = calibparam[1]
    gain = calibparam[2]
    vADC = calibparam[3]
  }
  if(length(calibparam) == 1) {
    if(calibparam == "SM4") {
      micsens = -35
      gain = 48
      vADC = 1
    }
    if(calibparam == "SMmini") {
      micsens = -11
      gain = 18
      vADC = 1.5
    }
    if(calibparam == "SM2+") {
      micsens = -36
      gain = 48
      vADC = 1.41
    }
    if(!calibparam %in% c("SM4", "SMmini", "SM2+")){
      stop("The type of calibration parameters specified is incorrect.")
    }
  }

  if(length(calibparam) != 3 & length(calibparam) != 1 ) {
    stop("Calibration parameters misconfigured, check the function's help page")
  }

if(anyNA(c(micsens, gain, vADC))){
  stop("Some of the calibration parameters have not been defined. Ensure that you have correctly specified the arguments: micsens, gain, and vADC.")
}

  if(!channel %in% c("left", "right")){
    stop("Incorrect channel definition, it must be one of the following options: 'left', 'right' ")
  }

valid_codes <- indices_codes

if ("all" %in% sel.ind) {
  sel.ind <- valid_codes
  if (!identical(exclude, "none")) {
    bad <- exclude[!exclude %in% valid_codes]
    if (length(bad)) {
      stop("Unknown indices in 'exclude': ", paste(bad, collapse = ", "))
    }
    sel.ind <- sel.ind[!sel.ind %in% exclude]
  }
} else {
  bad <- sel.ind[!sel.ind %in% valid_codes]
  if (length(bad)) {
    stop("Unknown indices in 'sel.ind': ", paste(bad, collapse = ", "))
  }
}


options_noise.ind = c("msldB_bio", "msldB_low", "mdBGL", "SpecIQR", "SpecKurt", "SpecSkew", "SpecQ3","SpecCent", "S2N_chi", "dfreq", "efb")

if(noise.ind) {
  sel.ind <- append(x = sel.ind, values = options_noise.ind[!options_noise.ind %in% sel.ind])
}

if(!parallel %in% c("files","minutes")) stop("parallel must be 'files' or 'minutes'")

# ---------------------------------
# Data management
# ---------------------------------


  time1 = proc.time()
  df <- data.frame()

  # Load file names
  files <- list.files(path = dir, pattern = "wav$", ignore.case = T, full.names = T)

  if (!length(files)) {
    warning("No WAV files were found in the specified directory.")
    return(data.frame())
  }
  filenumber <- seq_along(files)


  # Habilita el manejo de progreso

  # Working on multiple cores
  noCores <- parallel::detectCores(all.tests = TRUE, logical = FALSE)
  noCores <- noCores - 1
  if(noCores < 1) {
    noCores = 1
  }

  if(!is.null(ncores)){noCores = ncores}


  if(parallel == "files") {
    #if handlers not defined, use default aneco_progress_files
    if(is.null(names(progressr::handlers()))){progressr::handlers(aneco_progress_files)}
    progress_files <- progressr::progressor(along = filenumber)
    # Selecciona automáticamente la mejor opción según el sistema operativo
    if (.Platform$OS.type == "windows") {
      future::plan(multisession, workers = noCores)  # Para Windows
    } else {
      future::plan(multicore, workers = noCores)  # Para Mac y Linux
    }}else{
      future::plan(sequential)
    }

  ap <- future.apply::future_lapply(X = filenumber, FUN = function(xi){

    if(parallel == "files") {
      progress_files()
    }

    file <- files[xi]

#    message("Processing file ", xi, " of ", length(filenumber))

    name <- basename(as.character(file))
    nameFull <- file

    check_audio <- function(filepath) {
      tryCatch({
        audio <- tuneR::readWave(file, header = TRUE)
        return(TRUE)
      }, error = function(e) {
        return(FALSE)
      })
    }

    file_status <- check_audio(file)

    if(file_status){

      wavFile <- tuneR::readWave(file, header = TRUE)

      durS <- wavFile$samples/wavFile$sample.rate
      dur <- floor(durS/60)

      minutos <- seq(from = 0, to = dur, by = 1)
      mmin <- max(minutos)

#      progress_minutos <- progressr::progressor(along = minutos)

      if(parallel == "minutes"){
        #if handlers not defined, use default aneco_progress_minutes
        if(is.null(names(progressr::handlers()))){progressr::handlers(aneco_progress_minutes)}
        progress_minutes <- progressr::progressor(along = minutos)
        # Selecciona automáticamente la mejor opción según el sistema operativo
        if (.Platform$OS.type == "windows") {
          future::plan(multisession, workers = noCores)  # Para Windows
        } else {
          future::plan(multicore, workers = noCores)  # Para Mac y Linux
        }
      }else{
        future::plan(sequential)
      }

        mincl <- future.apply::future_lapply(X = minutos, FUN = function(i){

          if(parallel == "minutes") {
            progress_minutes(message = paste0("Processing file ", xi, " of ", length(filenumber)))
          }

#        progress_minutos(message = paste0("Processing file ", xi, " of ", length(filenumber), "  "))

        if(i == mmin){
          absdur <- durS - (i*60)
          if(absdur < 30){return(NULL)}
          wav <- tuneR::readWave(nameFull, from = (i*60), to = durS, units = "seconds")

        }else{
          wav <- tuneR::readWave(nameFull, from = (i*60), to = (i*60) + 60, units = "seconds")
          absdur <- 60
        }

        wav <- tuneR::mono(wav, which = channel)
        sampling_rate <- wav@samp.rate

        tempWav <- tempfile(fileext = ".wav")
        tuneR::writeWave(wav, filename = tempWav)



        # Calibration process
        # PAMCalib uses a pressure reference in microPa, while recorder
        # sensitivities supplied through calibparam are expressed in V/Pa.
        # 1 Pa = 1e6 microPa, so the equivalent sensitivity is 120 dB lower.
        datCalib <- ANECO::PAMCalib(atype = 'PSD', N = wl, timestring = "", r = 0, outwrite = 1, plottype = "None", calib = 1, envi = "Air", ctype = "TS", Mh = micsens - 120, G = gain, vADC = vADC, tempWav = tempWav, ifile = "", welch = "")


        unlink(tempWav)

        # Integrate PSD (power per Hz) into power per FFT bin before
        # computing band levels, A-weighted levels, and calibrated indices.
        bin_width <- sampling_rate / round(wl)
        datCalib[-1, -1] <- datCalib[-1, -1, drop = FALSE] +
          10 * log10(bin_width)

        #Compute dBA
        f_dBA <- datCalib[1, 2:ncol(datCalib)]
        # Decibel reference
        aweight <- vector()
        for ( w in 1:length(f_dBA)) {
          aweight[w] <- seewave::dBweight(f_dBA[w])$A
        }
        a <- datCalib[2:nrow(datCalib),2:ncol(datCalib)]
        aA = t(t(a) + aweight)


        # convert to pressure
        press <- rowSums(10^(aA/10))
        dBA = 10*log10(press) #hist(dBA)

        #Data manipulation
        FullMatrixA = as.data.frame(aA)
        rm(aA)
        FullMatrix <- as.data.frame(datCalib)[-1, -1]
        colnames(FullMatrix) <- round(datCalib[1, 2:ncol(datCalib)])
        colnames(FullMatrixA) <- round(datCalib[1, 2:ncol(datCalib)])
        rm(datCalib)
        vFreq <- as.numeric(colnames(FullMatrix))

        # Select "Biological" frequency bands
        PosBio <- which(vFreq >= bioband[1] & vFreq <= bioband[2])
        PosBio_standard <- which(vFreq >= 2000 & vFreq <= 8000)
        BioMatrix <- FullMatrix[, PosBio]
        BioMatrixA <- FullMatrixA[, PosBio]
        BioMatrix_standard <- FullMatrix[, PosBio_standard]


        # Select anthropic frequency bands
        PosNoise <- which(vFreq >= noiseband[1] & vFreq < noiseband[2])
        PosNoise_standard <- which(vFreq >= 0 & vFreq < 2000)
        AntMatrix <- FullMatrix[, PosNoise]
        AntMatrixA <- FullMatrixA[, PosNoise]
        AntMatrix_standard <- FullMatrix[, PosNoise_standard]


        # Number of frequency bands
        nFQ = as.numeric( dim(BioMatrix)[2])
        nFQbk =  as.numeric( dim(AntMatrix)[2])
        # Number of time frames
        fileDur = as.numeric( dim(BioMatrix)[1])

        # Matrix of sound pressure values
        BioMatrixPress <- 10^(BioMatrix/10)
        BioMatrixPress_standard <- 10^(BioMatrix_standard/10)
        AntMatrixPress <- 10^(AntMatrix/10)
        AntMatrixPress_standard <- 10^(AntMatrix_standard/10)

        normdBA <- function (x) { (x-(-10)) / (80 - (-10)) }


        # Function to calculate background noise
        # Modified here as a level for each frequency band

        # recordar que se puede quitar lo de los minimos
        BGL <- function(datos) {
          output <- NULL
          for (ff in 1:dim(datos)[2]) { #loop through each frequency band
            #(1) compute a histogram
            pretmp <- sort(datos[,ff]) #ordenar
            tmp <- pretmp[pretmp <= stats::quantile(pretmp, 0.99)] # excludes the top 1% of the values to avoid outliers
            if(min(tmp) == -Inf){
              output <- 0
              break()
            }
            brks = seq(from = min(tmp), to = max(tmp), length.out = length(tmp)/8)
            histo <- graphics::hist(tmp, breaks = brks, plot = FALSE)
            tab <- data.frame(counts = histo[["counts"]], mids = histo[["mids"]])
            tmpval <- tab[,1]
            tmpnam <- tab[,2]
            colmax <- as.numeric(which.max(tmpval))

            #(4) accumulating counts in histogram bins below (3) until 68% of counts below (3) is reached
            cutoff = sum( tmpval[ 1:colmax-1] ) * .68 #value to stop the summation at
            # got an error when the first bin had the max # of values..... so just  use the min value
            if (cutoff == 0) {stploc = 0}
            cntbk <- 0
            for (h in 1:length(tmpval[1:colmax])-1 )
            {
              if (cntbk > cutoff) {
                break }
              loc = colmax - h #find the location to start the summation
              cntbk = cntbk + tmpval[loc]
              stploc = tmpnam[loc-1]
            }
            # (5) final calc: (3) + N(4)
            if (cutoff == 0) {output[ff] = tmpnam[colmax] } else { output[ff] = tmpnam[colmax] + stploc*0.1 }
          }
          return(output)
        }

        # Background noise for each frequency band
        BK_Towsey <- BGL(datos = BioMatrix)
        BK_Towsey_anth <- BGL(datos = AntMatrix)

        if (sum(BK_Towsey) == 0 | sum(BK_Towsey_anth) == 0){
          return(NULL)
        }

        # --------------------------------------------
        # Median background level
        # --------------------------------------------
        median_BGL <- NA
        if("mdBGL" %in% sel.ind) {
          median_BGL <- stats::median(BK_Towsey)
        }

        # --------------------------------------------
        # Acoustic complexity index (ACI) calibrado
        # --------------------------------------------
        ACIout = NA

        if("ACI" %in% sel.ind){
          t <- array(0, dim(BioMatrix) - c(1, 0))
          for (frq in 1:(dim(t)[2])){
            t[, frq] <- abs(diff(BioMatrixPress[, frq])) / sum(BioMatrixPress[, frq])
          }
          ACIout = sum(t)
          rm(t)
        }

        # --------------------------------------------
        # Mean sound level
        #----------------------------------------------
        msldBA_low = NA
        msldBA_bio = NA
        msldB_low = NA
        msldB_bio = NA

        if("msldBA_low" %in% sel.ind){
          # dBA
          msldBA_low = 10*log10( sum( 10^(AntMatrixA/10))/ nrow(AntMatrixA))
        }
        if("msldBA_bio" %in% sel.ind){
          msldBA_bio = 10*log10( sum( 10^(BioMatrixA/10))/ nrow(BioMatrixA) )
        }

        if("msldB_low" %in% sel.ind){
          # dB
          msldB_low  = 10*log10( sum( AntMatrixPress)/ nrow(AntMatrix))
        }
        if("msldB_bio" %in% sel.ind){
          msldB_bio = 10*log10( sum(BioMatrixPress)/ nrow(BioMatrix) )
        }

        # --------------------------------------------
        #  Average signal amplitude
        # --------------------------------------------
        avgAMP = NA

        #Towsey et al., 2013
        # modified as the mean dBA value for the timestep, normalized by the min/max (set at: -10 to 80 dB)
        # modified to fit data

        if("avgAMP" %in% sel.ind){
          avgAMP = normdBA(mean(dBA))
        }

        # --------------------------------------------
        # L10 exceedance level
        # --------------------------------------------
        L10AMP = NA
        if("L10AMP" %in% sel.ind){
          L10AMP = normdBA(quantile (dBA, .1) ) #10 porcentil
        }

        # --------------------------------------------
        # Entropy indices, code developed by Buxton et al. 2018
        # --------------------------------------------

        Hf = NA
        Ht = NA
        EI = NA
        Pres = colMeans(BioMatrixPress)  #Leq for each Fq band over time period (mean of the pressures)
        Pres2 = Pres/sum(Pres)
        if("Hf" %in% sel.ind){
          Hf = -sum( (Pres2 * log(Pres2))) /log(nFQ)
        }
        if("Ht" %in% sel.ind){
          Leqt = rowMeans( BioMatrixPress) #Leq for each second over entire band (could also used dBA)
          Leqt2 = Leqt/sum(Leqt)
          Ht = -sum( (Leqt2 * log(Leqt2)) / log(fileDur) )
        }
        if("EI" %in% sel.ind){
          if(is.na(Hf)){
            Hf = -sum( (Pres2 * log(Pres2))) /log(nFQ)
          }
          if(is.na(Ht)){
            Leqt = rowMeans( BioMatrixPress)
            Leqt2 = Leqt/sum(Leqt)
            Ht = -sum( (Leqt2 * log(Leqt2)) / log(fileDur) )
            }
          EI = Ht * Hf
        }


        # --------------------------------------------
        # Relationship between biophonies and anthropophonies
        # --------------------------------------------
        NDSI = NA
        BioPh = NA
        AntPh = NA
        Bio_anth = NA

        if("BioPh" %in% sel.ind){
          BioPh = sum(colMeans(BioMatrixPress_standard))
        }
        if("AntPh" %in% sel.ind){
          AntPh = sum(colMeans(AntMatrixPress_standard))
        }
        if("NDSI" %in% sel.ind){
          if(is.na(BioPh)){
            BioPh = sum(colMeans(BioMatrixPress_standard))
          }
          if(is.na(AntPh)){
            AntPh = sum(colMeans(AntMatrixPress_standard))
          }
          NDSI = (BioPh -AntPh ) / (BioPh + AntPh )
        }
        if("Bio_anth" %in% sel.ind){
          if(is.na(BioPh)){
            BioPh = sum(colMeans(BioMatrixPress_standard))
          }
          if(is.na(AntPh)){
            AntPh = sum(colMeans(AntMatrixPress_standard))
          }
          Bio_anth = (BioPh /AntPh)
        }

        # --------------------------------------------
        # Difference with background noise
        # --------------------------------------------

        AAanth = NA
        AA = NA
        AAcanth = NA
        AAc = NA
        AAduranth = NA
        AAdur = NA

        if("AA" %in% sel.ind){
          for (f in 1:length(BioMatrix) ) {
            AA[f] =  length( BioMatrix[, f][ BioMatrix[, f] > BK_Towsey[f] + 3 ] )/ length (BioMatrix[, f]) #proporcion
          }
          AA = sum(AA)
        }
        if("AAc" %in% sel.ind){
          AAc = NULL
          for (f in 1:length(BioMatrix) ) {
            AAc[f] = (length( BioMatrix[, f][ BioMatrix[, f] > BK_Towsey[f] + 3 ] ))
          }
          AAc = sum(AAc)
        }

        if("AAdur" %in% sel.ind){
          AAdur = NULL
          for (f in 1:length(BioMatrix) ) {
            #logical matrix...
            temp <- BioMatrix[,f] > BK_Towsey[f] + 3

            temp2 <- temp*1 # convertir a 0/1
            #encontrar numeros consecutivos
            rl <- rle(temp2)
            len = rl$lengths
            v =  rl$value
            if (length(v) == 1 )
            {
              if   (v == 0) {AAdur[f] = 0 }
              else {AAdur[f] = len }
              next
            }
            cumsum = NULL
            cntAA = 0
            for ( qq in  seq(from = 2, to = length(v), by =2) )
            {
              cntAA = cntAA + 1
              cumsum[cntAA] = len[qq]
            }
            AAdur[f] = mean(cumsum)
          }
          AAdur = median(AAdur)
        }


          if("AAanth" %in% sel.ind){
            AAanth = NULL
            for (f in 1:length(AntMatrix) ) {
              AAanth[f] =  length( AntMatrix[, f][ AntMatrix[, f] > BK_Towsey[f] + 3 ] )/ length (AntMatrix[, f]) #proporcion
            }
            AAanth = sum(AAanth)
          }
          if("AAcanth" %in% sel.ind){
            AAcanth = NULL
            for (f in 1:length(AntMatrix) ) {
              AAcanth[f] = (length( AntMatrix[, f][ AntMatrix[, f] > BK_Towsey[f] + 3 ]))
            }
            AAcanth = sum(AAcanth)
          }

          if("AAduranth" %in% sel.ind){
            AAduranth = NULL
            for (f in 1:length(AntMatrix) ) {
              #logical matrix...
              temp <- AntMatrix[,f] > BK_Towsey[f] + 3

              temp2 <- temp*1 # convertir a 0/1
              #encontrar numeros consecutivos
              rl <- rle(temp2)
              len = rl$lengths
              v =  rl$value
              if (length(v) == 1 )
              {
                if   (v == 0) {AAduranth[f] = 0 }
                else {AAduranth[f] = len }
                next
              }
              cumsum = NULL
              cntAA = 0
              for ( qq in  seq(from = 2, to = length(v), by =2) )
              {
                cntAA = cntAA + 1
                cumsum[cntAA] = len[qq]
              }
              AAduranth[f] = mean(cumsum)
            }
            AAduranth = median(AAduranth)
          }

        # --------------------------------------------
        # Roughness
        # - ind.Rough = TRUE
        # --------------------------------------------

        Rough = NA

        if("Rough" %in% sel.ind){
          Rough = NULL
          for (f in 1:length(BioMatrixPress) ){
            x = BioMatrixPress[, f]
            x <- x/max(x)
            deriv2 <- diff(x, 1, 2)
            Rough[f] <- sum(deriv2^2, na.rm = TRUE)
          }
          Rough = median(Rough)
        }

        # --------------------------------------------
        # Acoustic diversity, acoustic eveness
        # - sel.ind = c("adi", "aei")
        # --------------------------------------------

        ADI = NA
        AEI = NA

        f_prop = NULL

        if("ADI" %in% sel.ind){
          rd <- floor(as.numeric(names(BioMatrix))/1000)

          fs <- unique(rd)
          f_perc <- NA
          c1 = 1

          for (frqb in fs) {
            scorefi = NA
            c2 = 1
            fb <- which(rd == frqb)
            for (fi in fb) {
              scorefi[c2] <- length(BioMatrix[, fi][BioMatrix[, fi] > (max(BioMatrix) -50)])
              c2 <- c2 + 1
            }
            f_perc[c1] = sum(scorefi)/(length(fb)*nrow(BioMatrix))

            c1 <- c1 + 1
          }

          f_prop <- f_perc/sum(f_perc)

          ADI = vegan::diversity(f_prop, index = "shannon")
        }
        if("AEI" %in% sel.ind){
          if(is.null(f_prop)){
            rd <- floor(as.numeric(names(BioMatrix))/1000)

            fs <- unique(rd)
            f_perc <- NA
            c1 = 1

            for (frqb in fs) {
              scorefi = NA
              c2 = 1
              fb <- which(rd == frqb)
              for (fi in fb) {
                scorefi[c2] <- length(BioMatrix[, fi][BioMatrix[, fi] > (max(BioMatrix) -50)])
                c2 <- c2 + 1
              }
              f_perc[c1] = sum(scorefi)/(length(fb)*nrow(BioMatrix))

              c1 <- c1 + 1
            }
          }
          f_prop <- f_perc/sum(f_perc)
          AEI = ineq::Gini(f_prop)
        }

        # --------------------------------------------
        # Energy distribution in the spectrum
        #  - ind.freqdist = TRUE
        # --------------------------------------------

        pk = NA
        pkd = NA
        pks = NA
        Hm = NA
        HvPres = NA
        HvSPL = NA


        if(any(c("pk", "pkd", "pks", "Hm") %in% sel.ind)){
          peakf = NULL
          for (j in 1:dim(BioMatrix)[1] )
          {   peakf[j] = (which.max( BioMatrix[j, ] ) )  }
          pk2 = matrix(0, 1, dim(BioMatrix)[2])
          for (uu in 1:nFQ)  { pk2[uu] = sum(peakf == uu)  }
          colnames(pk2) = colnames(BioMatrix)
          pk2nor = pk2/fileDur
        }

         if("pk" %in% sel.ind){
          pk = as.numeric(gsub("H", "", colnames(BioMatrix[which.max(pk2)])))
         }

        if("pkd" %in% sel.ind){
          # kurtosis
          pkd = moments::kurtosis(as.vector(pk2nor))
        }

        if("pks" %in% sel.ind){
          # skewness
          pks = moments::skewness(as.vector(pk2nor))
        }
        if("Hm" %in% sel.ind){
        # Entropy of Spectral Maxima
          pk2[pk2 == 0] <- 1e-07
          pk_prob = pk2/(sum(pk2)) # normalize
          Hm = -sum((pk_prob * log2(pk_prob))) / log2(nFQ)
        }

        if("HvPress" %in% sel.ind){
          # Entropy of spectral variance- pressure and dB
          Press = (BioMatrixPress)
          Pv = NULL
          for (v in 1:dim(Press)[2] ) {Pv[v] = stats::var(Press[ , v]) }
          Pv2 = Pv/sum(Pv)
          HvPres = -sum( (Pv2 * log2(Pv2)) ) / log2(nFQ)
        }

        if("HvSPL" %in% sel.ind){
          Pv = NULL
          for (v in 1:dim(BioMatrix)[2] ) {Pv[v] = var(BioMatrix[,v])  }
          Pv2 = Pv/sum(Pv)
          HvSPL = -sum( (Pv2 * log2(Pv2)) ) / log2(nFQ)
        }


        # --------------------------------------------
        # Normalize exceedence levels for dBA
        # --------------------------------------------
        dif_L10L90 = NA

        if("dif_L10L90" %in% sel.ind){
          Exceed_norm = normdBA(quantile(dBA, c(0.05, 0.1, 0.5, 0.9, 0.95),na.rm = TRUE))
          L10n = (Exceed_norm[4]) #L10 = Porcentil 90
          L90n = (Exceed_norm[2]) #L90 = Porcentil 10
          dif_L10L90 = L10n - L90n
        }

        # --------------------------------------------
        # Median sound level
        # --------------------------------------------
        Mamp = NA

        if("Mamp" %in% sel.ind){
          Mamp = quantile( 10*log10(rowMeans(BioMatrixPress)), 0.5)
        }

        # --------------------------------------------
        # Indices for rain and cicada detection
        # --------------------------------------------

        ind_noise = NA
        if("S2N_chi" %in% sel.ind){
          # Prediccion de ruidos
          fr_chi <- which(vFreq >= 3500 & vFreq <= 5500)
          BioMatrix_chi <- FullMatrix[, fr_chi]
          q75_chi <- quantile(as.matrix(10^(BioMatrix_chi/10)), 0.75)
          BG_amp <- 10^(median(BK_Towsey)/10)
          ind_noise <- (q75_chi - BG_amp)/(q75_chi + BG_amp)
        }

        # Frequence energy
        freq.energy <- NA
        if(any(c("efb") %in% sel.ind)){
          freq.energy <- list()
          fqmax <- max(vFreq)
          lband <- 1000
          countpos <- 1
          while (lband <= fqmax) {
            nband <- which(vFreq >= (lband - 1000) & vFreq < lband)
            fband <- FullMatrix[, nband]
            energydB = 10 * log10(sum(10^(fband/10))/nrow(fband))
            freq.energy[[countpos]] <- energydB
            names(freq.energy)[countpos] <- paste0("E", countpos - 1, countpos)
            countpos = countpos + 1
            lband = lband + 1000
          }
        }

        #///////////////////////////////////////////////////////////
        # The following are indices that do not go through the calibration process
        #//////////////////////////////////////////////////////////

        # Extract the corresponding audio minute

        if(any(c("ADI_soundecology", "AEI_soundecology", "ACI_soundecology", "ACI_seewave", "BIO_soundecology", "NDSI_soundecology", "MAE_seewave", "Hf_seewave", "Ht_seewave", "TE_seewave", "rough_seewave", "NP", "SpecCent", "SpecIQR", "SpecQ3",  "SpecSkew", "SpecKurt", "dfreq", "efb", "dfb", "S2N_chi") %in% sel.ind)){

        if(i == max(minutos)){
          absdur <- durS - (i*60)
          if(absdur < 30){return(NULL)}
          minwav <- tuneR::readWave(nameFull, from = (i*60), to = durS, units = "seconds")
        }else{
          minwav <- tuneR::readWave(nameFull, from = (i*60), to = (i*60) + 60, units = "seconds")
        }


        minwav <- tuneR::channel(minwav, which = channel)

        }



        # --------------------------------------------
        # Dominant frequencies
        # --------------------------------------------

        Prom.Dom.frec = NA
        length.fd.01 <- NA
        length.fd.12 <- NA
        length.fd.23 <- NA
        length.fd.23 <- NA
        length.fd.34 <- NA
        length.fd.45 <- NA
        length.fd.56 <- NA
        length.fd.67 <- NA
        length.fd.78 <- NA
        length.fd.89 <- NA
        length.fd.910 <- NA
        length.fd.910 <- NA
        length.fd.1011 <- NA
        length.fd.11 <- NA

        fmin_domfreq <- 100
        fmax_domfreq <- ifelse(sampling_rate/2 >= 11000, 11000, sampling_rate/2)

        if(any(c("dfb", "dfreq")  %in% sel.ind)){
        dom.freq <- as.data.frame(seewave::dfreq(minwav, channel = 1, f = minwav@samp.rate, wl = 512, threshold = 0.015, plot = FALSE, bandpass = c(fmin_domfreq, fmax_domfreq)))

          Prom.Dom.frec <- mean(dom.freq$y, na.rm = TRUE)

          if("dfb" %in% sel.ind){
          dm.01 <- dom.freq[dom.freq$y <= 1, ]
          dm.12 <- dom.freq[dom.freq$y > 1 & dom.freq$y <= 2, ]
          dm.23 <- dom.freq[dom.freq$y > 2 & dom.freq$y <= 3, ]
          dm.34 <- dom.freq[dom.freq$y > 3 & dom.freq$y <= 4, ]
          dm.45 <- dom.freq[dom.freq$y > 4 & dom.freq$y <= 5, ]
          dm.56 <- dom.freq[dom.freq$y > 5 & dom.freq$y <= 6, ]
          dm.67 <- dom.freq[dom.freq$y > 6 & dom.freq$y <= 7, ]
          dm.78 <- dom.freq[dom.freq$y > 7 & dom.freq$y <= 8, ]
          dm.89 <- dom.freq[dom.freq$y > 8 & dom.freq$y <= 9, ]
          dm.910 <- dom.freq[dom.freq$y > 9 & dom.freq$y <= 10, ]
          dm.1011 <- dom.freq[dom.freq$y > 10 & dom.freq$y <= 11, ]
          dm.11 <- dom.freq[dom.freq$y > 11, ]

          length.fd.01 <- length(dm.01$y)
          length.fd.12 <- length(dm.12$y)
          length.fd.23 <- length(dm.23$y)
          length.fd.23 <- length(dm.23$y)
          length.fd.34 <- length(dm.34$y)
          length.fd.45 <- length(dm.45$y)
          length.fd.56 <- length(dm.56$y)
          length.fd.67 <- length(dm.67$y)
          length.fd.78 <- length(dm.78$y)
          length.fd.89 <- length(dm.89$y)
          length.fd.910 <- length(dm.910$y)
          length.fd.1011 <- length(dm.1011$y)
          length.fd.11 <- length(dm.11$y)
          }
        }

        # --------------------------------------------
        # Traditional acoustic indices
        # --------------------------------------------

        AEIprom = NA
        ADIprom = NA
        ACIprom = NA
        ACIseewave = NA
        BIOprom = NA
        NDSIprom = NA
        TEraw = NA
        Ht_raw = NA
        Hf_raw = NA
        Rough_see = NA
        MAE = NA
        NP = NA


        if("AEI_soundecology" %in% sel.ind){
          capture.output(AEIraw <- soundecology::acoustic_evenness(minwav), file = nullfile())
          AEIprom <- AEIraw$aei_left
        }
        if("ADI_soundecology" %in% sel.ind){
          capture.output(ADIraw <- soundecology::acoustic_diversity(minwav), file = nullfile())
          ADIprom <- ADIraw$adi_left
        }
        if("ACI_soundecology" %in% sel.ind){
          capture.output(ACIraw <- soundecology::acoustic_complexity(minwav), file = nullfile())
          ACIprom <- ACIraw$AciTotAll_left
        }
        if("ACI_seewave" %in% sel.ind){
          ACIseewave <- seewave::ACI(minwav)
        }

        if("BIO_soundecology" %in% sel.ind){
          capture.output(BIOraw <- soundecology::bioacoustic_index(minwav), file = nullfile())
          BIOprom <- BIOraw$left_area
        }

        if("NDSI_soundecology" %in% sel.ind){
          capture.output(NDSIraw <- soundecology::ndsi(minwav), file = nullfile())
          NDSIprom <- NDSIraw$ndsi_left
        }

        if("TE_seewave" %in% sel.ind){
          TEraw <- seewave::H(minwav)
        }

        if("Ht_seewave" %in% sel.ind){
          envorni <- seewave::env(minwav, f = minwav@samp.rate, plot = FALSE)
          Ht_raw <- seewave::th(envorni)
        }

        if(any(c("Hf_seewave", "rough_seewave") %in% sel.ind)){
          suppressWarnings(speca <- seewave::spec(minwav, f = minwav@samp.rate, plot = FALSE))
        }
        if("Hf_seewave" %in% sel.ind){
          Hf_raw <- seewave::sh(speca)
        }

        if("rough_seewave" %in% sel.ind){
          Rough_see <- seewave::roughness(speca) #diferente por defecto, usa todo el espectro
        }

        if("MAE_seewave" %in% sel.ind){
          MAE <- seewave::M(minwav)
        }


        if(any(c("NP", "SpecCent", "SpecIQR", "SpecQ3",  "SpecSkew", "SpecKurt") %in% sel.ind)){
          specm <- seewave::meanspec(minwav, plot = FALSE)
        }

        if("NP" %in% sel.ind){
          peaks <- seewave::fpeaks(specm, plot = FALSE)
          NP <- length(peaks)/2
        }

        # --------------------------------------------
        # Spectral indices
        # --------------------------------------------

        cent = NA
        iqr = NA
        tq = NA
        skew = NA
        kurt = NA

        if(any(c("SpecCent", "SpecIQR", "SpecQ3",  "SpecSkew", "SpecKurt") %in% sel.ind)){
          specStats <- seewave::specprop(specm, f = minwav@samp.rate, flim = c(fmin_domfreq/1000, fmax_domfreq/1000), plot = FALSE)
        }

        if("SpecCent" %in% sel.ind){
          cent <- specStats$cent
        }
        if("SpecIQR" %in% sel.ind){
          iqr <- specStats$IQR
        }
        if("SpecQ3" %in% sel.ind){
          tq <- specStats$Q75
        }
        if("SpecSkew" %in% sel.ind){
          skew <- specStats$skewness
        }
        if("SpecKurt" %in% sel.ind){
          kurt <- specStats$kurtosis
        }

        # -------------------------------------
        # Identification data
        #--------------------------------------
        namesplit <- unlist(strsplit(name, "_"))
        if(length(namesplit == 3)){
          conwav = namesplit[3]
          hora <- as.character(unlist(strsplit(conwav, ".wav"))[1])
          codigo <- namesplit[1]
          sitioindex <- unlist(gregexpr('S', prefix.format))
          grabadoraindex <- unlist(gregexpr('R', prefix.format))
          sitio <- substr(codigo, start = min(sitioindex), stop = max(sitioindex))
          grabadora <- substr(codigo, start = min(grabadoraindex), stop = max(grabadoraindex))
          fecha <- namesplit[2]
          ano <- substr(fecha, start = 1, stop = 4)
          mes <- substr(fecha, start = 5, stop = 6)
          dia <- substr(fecha, start = 7, stop = 8)
        }

        if(length(namesplit) == 2) {
          name2 <- paste(add.prefix, name, sep = "_")
          namesplit <- unlist(strsplit(name2, "_"))
          conwav = namesplit[3]
          hora <- as.character(unlist(strsplit(conwav, ".wav"))[1])
          codigo <- namesplit[1]
          sitioindex <- unlist(gregexpr('S', prefix.format))
          grabadoraindex <- unlist(gregexpr('R', prefix.format))
          sitio <- substr(codigo, start = min(sitioindex), stop = max(sitioindex))
          grabadora <- substr(codigo, start = min(grabadoraindex), stop = max(grabadoraindex))
          fecha <- namesplit[2]
          ano <- substr(fecha, start = 1, stop = 4)
          mes <- substr(fecha, start = 5, stop = 6)
          dia <- substr(fecha, start = 7, stop = 8)
        }

        z <- list(Name = name, Code = codigo, Site = sitio, RecorderID = grabadora, Date = fecha, Year = ano, Month = mes, Day = dia, Hour = hora, Minute = i, Duration = absdur, TE_seewave = TEraw,  Ht_seewave = Ht_raw, Hf_seewave = Hf_raw, MAE_seewave = MAE, NP = NP, ACI_soundecology = ACIprom, ACI_seewave = ACIseewave, rough_seewave = Rough_see, ADI_soundecology = ADIprom, AEI_soundecology = AEIprom, BIO_soundecology = BIOprom, NDSI_soundecology = NDSIprom, msldBA_low = msldBA_low, msldBA_bio = msldBA_bio, msldB_low = msldB_low, msldB_bio = msldB_bio, BioPh = BioPh, AntPh = AntPh, avgAMP = avgAMP, L10AMP = L10AMP, AAanth = AAanth, AA = AA, AAcanth = AAcanth, AAc = AAc, AAduranth = AAduranth, AAdur = AAdur, pk = pk, pkd = pkd, pks = pks, Hf = Hf,Ht = Ht, NDSI = NDSI, Bio_anth = Bio_anth,  rough = Rough, EI = EI, Hm = Hm, HvPres = HvPres, HvSPL = HvSPL, ACI = ACIout, ADI = ADI, AEI = AEI, Mamp = Mamp, dif_L10L90 = dif_L10L90,  mdBGL = median_BGL, SpecCent = cent, SpecIQR = iqr, SpecQ3 = tq, SpecSkew = skew, SpecKurt = kurt, S2N_chi = ind_noise, dfreq = Prom.Dom.frec, freq.energy, F01 = length.fd.01, F12 = length.fd.12, F23 = length.fd.23, F34 = length.fd.34, F45 = length.fd.45, F56 = length.fd.56, F67 = length.fd.67, F78 = length.fd.78, F89 = length.fd.89, F910 = length.fd.910, F1011 = length.fd.1011, F11 = length.fd.11)

        # Creating the result dataframe

        return(data.frame(z))
      }, future.seed = TRUE, future.scheduling = 1)
      future::plan(sequential)
      df <- as.data.frame(dplyr::bind_rows(mincl))
    }else{
      df <- data.frame(Name = name, Code = NA, Site = NA, RecorderID = NA, Date = NA, Year = NA, Month = NA, Day = NA, Hour = NA, Minute = NA, Duration = NA, TE_seewave = NA,  Ht_seewave = NA, Hf_seewave = NA, MAE_seewave = NA, NP = NA, ACI_soundecology = NA, ACI_seewave = NA, rough_seewave = NA, ADI_soundecology = NA, AEI_soundecology = NA, BIO_soundecology = NA, NDSI_soundecology = NA, msldBA_low = NA, msldBA_bio = NA, msldB_low = NA, msldB_bio = NA, BioPh = NA, AntPh = NA, avgAMP = NA, L10AMP = NA, AAanth = NA, AA = NA, AAcanth = NA, AAc = NA, AAduranth = NA, AAdur = NA, pk = NA, pkd = NA, pks = NA, Hf = NA,Ht = NA, NDSI = NA, Bio_anth = NA,  rough = NA, EI = NA, Hm = NA, HvPres = NA, HvSPL = NA, ACI = NA, ADI = NA, AEI = NA, Mamp = NA, dif_L10L90 = NA,  mdBGL = NA, SpecCent = NA, SpecIQR = NA, SpecQ3 = NA, SpecSkew = NA, SpecKurt = NA, S2N_chi = NA, dfreq = NA, freq.energy = NA, F01 = NA, F12 = NA, F23 = NA, F34 = NA, F45 = NA, F56 = NA, F67 = NA, F78 = NA, F89 = NA, F910 = NA, F1011 = NA, F11 = NA)
    }
    return(df)
  }, future.seed = TRUE, future.scheduling = 1)
  future::plan(sequential)
  nd <- as.data.frame(dplyr::bind_rows(ap))
  cat("!Analysis completed! \n")

  ptimeS <- (proc.time()-time1)[3]
  phours <- floor(ptimeS/60/60)
  pmin <- floor(ptimeS/60) - phours*60
  cat('Processing time: ', phours, "h ", pmin, 'min \n')

  rownames(nd) <- NULL
  if (!nrow(nd)) {
    warning("No valid minutes were processed.")
    return(nd)
  }

  # Exclude Name when checking failed rows, but inspect every value in columns.
  NArow <- apply(nd[, setdiff(names(nd), "Name"), drop = FALSE], 1,
                 function(x) all(is.na(x)))
  errores <- nd$Name[NArow]
  if (length(errores)) {
    warning("The following files produced errors: \n", paste(errores, collapse = ", "))
  }
  filtro2 <- nd[!NArow, , drop = FALSE]
  if (!nrow(filtro2)) return(filtro2)

  NAcol <- vapply(filtro2, function(x) all(is.na(x)), logical(1))
  filtro2 <- filtro2[, !NAcol, drop = FALSE]

  if(save.file) {
    codigoarchivo <- filtro2$Code[1]
    file_path <- file.path(dir, paste0(codigoarchivo, "_", format(Sys.time(), "%Y%m%d_%H%M") ,".csv"))
    utils::write.csv(filtro2, file_path,row.names = FALSE)
  }
  return(filtro2)
}
