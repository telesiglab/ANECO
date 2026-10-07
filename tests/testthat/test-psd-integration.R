test_that("band SPL and energy are independent of FFT bin width", {
  set.seed(456)
  audio_dirs <- character()
  on.exit(unlink(audio_dirs, recursive = TRUE), add = TRUE)
  dirs <- list()
  peak_pressure <- 0.001

  for (sample_rate in c(8192, 16384)) {
    time <- (seq_len(120 * sample_rate) - 1) / sample_rate
    envelope <- 1 + 0.1 * sin(2 * pi * time)
    pressure <- peak_pressure * envelope *
      (sin(2 * pi * 256 * time) + sin(2 * pi * 1536 * time)) +
      rnorm(length(time), sd = peak_pressure * 0.02)
    adc_scale <- 10^((-35 + 48) / 20)
    wave <- tuneR::Wave(
      left = as.integer(round(pressure * adc_scale * 2^15)),
      samp.rate = sample_rate, bit = 16
    )
    audio_dir <- tempfile("aneco-psd-integration-")
    dir.create(audio_dir)
    audio_dirs <- c(audio_dirs, audio_dir)
    dirs[[as.character(sample_rate)]] <- audio_dir
    tuneR::writeWave(wave, file.path(audio_dir, "AR01_20260101_000000.wav"))
  }

  configurations <- data.frame(sample_rate = c(8192, 8192, 8192, 16384),
                               wl = c(256, 512, 1024, 512))
  expected_rms <- peak_pressure * sqrt((1 + 0.1^2 / 2) / 2)
  expected_spl <- 20 * log10(expected_rms / 20e-6)
  expected_low_A <- expected_spl + seewave::dBweight(256)$A
  expected_bio_A <- expected_spl + seewave::dBweight(1536)$A
  means <- numeric(nrow(configurations))

  for (i in seq_len(nrow(configurations))) {
    configuration <- configurations[i, ]
    result <- acoustic_indices(
      dir = dirs[[as.character(configuration$sample_rate)]],
      calibparam = "SM4", wl = configuration$wl,
      sel.ind = c("msldB_bio", "msldB_low", "msldBA_bio", "msldBA_low", "efb"),
      noise.ind = FALSE, bioband = c(1000, 2000), noiseband = c(0, 1000),
      ncores = 1, save.file = FALSE
    )
    context <- paste("Fs =", configuration$sample_rate, "wl =", configuration$wl)
    expect_equal(nrow(result), 2L, info = context)
    columns <- c("msldB_bio", "msldB_low", "E01", "E12")
    levels <- as.matrix(result[, columns, drop = FALSE])
    expect_true(all(is.finite(levels)), info = context)
    expect_true(all(abs(levels - expected_spl) < 0.15), info = context)
    expect_true(all(abs(result$msldBA_low - expected_low_A) < 0.15), info = context)
    expect_true(all(abs(result$msldBA_bio - expected_bio_A) < 0.15), info = context)
    means[i] <- mean(result$msldB_bio)
  }
  expect_lt(diff(range(means)), 0.1)
})
