pam_test_spectrum <- function(signal, sample_rate = 8192, wl = 512,
                              lcut = sample_rate / wl,
                              hcut = sample_rate / 2, atype = "PSD") {
  file <- tempfile(fileext = ".wav")
  on.exit(unlink(file), add = TRUE)
  wave <- tuneR::Wave(left = as.integer(round(signal * 2^15)),
                     samp.rate = sample_rate, bit = 16)
  tuneR::writeWave(wave, file)
  PAMCalib(tempWav = file, ifile = "", atype = atype,
           calib = 1, ctype = "EE", Si = 0, envi = "Wat",
           N = wl, r = 0, lcut = lcut, hcut = hcut,
           plottype = "None")
}

test_that("FFT tone frequencies match their spectrum columns", {
  time <- (0:2047) / 8192
  spectrum <- pam_test_spectrum(0.25 * sin(2 * pi * 1024 * time))
  frequencies <- unname(spectrum[1, -1])
  expect_equal(frequencies, seq(16, 4096, by = 16))
  peaks <- max.col(spectrum[-1, -1, drop = FALSE], ties.method = "first")
  expect_equal(unname(frequencies[peaks]), rep(1024, 4))
})

test_that("DC and Nyquist are retained without doubling their power", {
  signals <- list(DC = rep(0.25, 2048),
                  Nyquist = 0.25 * (-1)^(0:2047))
  for (name in names(signals)) {
    spectrum <- pam_test_spectrum(signals[[name]], lcut = 0)
    frequencies <- unname(spectrum[1, -1])
    expect_equal(frequencies[c(1, length(frequencies))], c(0, 4096))
    peaks <- max.col(spectrum[-1, -1, drop = FALSE], ties.method = "first")
    expected_peak <- if (name == "DC") 0 else 4096
    expect_equal(unname(frequencies[peaks]), rep(expected_peak, 4))
    power <- rowSums(10^(spectrum[-1, -1, drop = FALSE] / 10)) * 16
    expect_equal(unname(power), rep(0.25^2, 4), tolerance = 1e-10)
  }
})

test_that("frequency cuts keep the matching FFT bins and matrix shape", {
  time <- (0:511) / 8192
  signal <- 0.25 * sin(2 * pi * 1024 * time)
  full <- pam_test_spectrum(signal)
  single <- pam_test_spectrum(signal, lcut = 1024, hcut = 1024)
  expect_equal(dim(single), c(2L, 2L))
  expect_equal(unname(single[1, 2]), 1024)
  expect_equal(single[2, 2], full[2, which(full[1, ] == 1024)])
  expect_error(pam_test_spectrum(signal, lcut = 1025, hcut = 1030),
               "No FFT frequencies")
})

test_that("odd FFT lengths retain and double the final positive bin", {
  wl <- 511
  time <- (0:(4 * wl - 1)) / 8192
  frequency <- 255 * 8192 / wl
  spectrum <- pam_test_spectrum(0.25 * sin(2 * pi * frequency * time),
                                wl = wl, lcut = 0)
  frequencies <- unname(spectrum[1, -1])
  expect_equal(length(frequencies), 256L)
  expect_equal(tail(frequencies, 1), frequency)
  peaks <- max.col(spectrum[-1, -1, drop = FALSE], ties.method = "first")
  expect_equal(unname(frequencies[peaks]), rep(frequency, 4))
  # At the final odd-length bin, Hann leakage overlaps its negative partner.
  # Compare against windowed time-domain power via Parseval's identity.
  digital <- round(0.25 * sin(2 * pi * frequency * time) * 2^15) / 2^15
  window <- 1 - cos(2 * pi * (1:wl) / wl)
  expected <- colMeans((matrix(digital, nrow = wl) * window)^2) /
    mean(window^2)
  power <- rowSums(10^(spectrum[-1, -1, drop = FALSE] / 10)) * 8192 / wl
  expect_equal(unname(power), unname(expected), tolerance = 1e-10)
})

test_that("odd sample rates preserve exact FFT frequencies", {
  sample_rate <- 8193
  time <- (0:2047) / sample_rate
  spectrum <- pam_test_spectrum(0.25 * sin(2 * pi * 64 * time * sample_rate / 512),
                                sample_rate = sample_rate)
  expect_equal(unname(spectrum[1, -1]), (1:256) * sample_rate / 512)
})

test_that("third-octave bands use the filtered spectrum frequencies", {
  time <- (0:2047) / 8192
  spectrum <- pam_test_spectrum(0.25 * sin(2 * pi * 1024 * time),
                                lcut = 500, hcut = 2000, atype = "TOL")
  frequencies <- unname(spectrum[1, -1])
  peaks <- max.col(spectrum[-1, -1, drop = FALSE], ties.method = "first")
  expect_equal(unname(frequencies[peaks]), rep(1000, 4))
  # PowerSpec/TOL use amplitude-corrected Hann power (bandwidth factor 1.5).
  band <- which.min(abs(frequencies - 1000)) + 1L
  power <- 10^(spectrum[-1, band] / 10)
  expect_length(power, 4L)
  expect_true(all(abs(power - 1.5 * 0.25^2 / 2) < 1e-5))
})
