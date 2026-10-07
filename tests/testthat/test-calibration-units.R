test_that("SPL calibration preserves preset and numeric inputs in V/Pa", {
  skip_if_not_installed("tuneR")
  skip_if_not_installed("future.apply")
  skip_if_not_installed("progressr")

  # A known pressure signal, converted to ADC samples using the SM4 settings.
  # Amplitude modulation and a small noise floor avoid constant/silent spectra.
  set.seed(123)
  sample_rate <- 8192
  time <- (seq_len(120 * sample_rate) - 1) / sample_rate
  peak_pressure <- 0.001 # Pa
  pressure <- peak_pressure * (1 + 0.1 * sin(2 * pi * time)) *
    sin(2 * pi * 1024 * time) +
    rnorm(length(time), sd = peak_pressure * 0.02)
  adc_scale <- 10^((-35 + 48) / 20)
  samples <- as.integer(round(pressure * adc_scale * 2^15))
  wave <- tuneR::Wave(left = samples, samp.rate = sample_rate, bit = 16)

  audio_dir <- tempfile("aneco-calibration-")
  dir.create(audio_dir)
  on.exit(unlink(audio_dir, recursive = TRUE), add = TRUE)
  tuneR::writeWave(wave, file.path(audio_dir, "AR01_20260101_000000.wav"))

  compute <- function(calibration) {
    acoustic_indices(
      dir = audio_dir, calibparam = calibration,
      sel.ind = c("msldB_bio", "efb"), noise.ind = FALSE,
      bioband = c(500, 1500), wl = 512,
      ncores = 1, save.file = FALSE
    )
  }
  preset <- compute("SM4")
  custom <- compute(c(-35, 48, 1))
  expect_equal(nrow(preset), 2L)
  expect_equal(preset$msldB_bio, custom$msldB_bio, tolerance = 1e-8)

  # Returned band levels must already include the PSD integration factor.
  expected_rms <- peak_pressure * sqrt((1 + 0.1^2 / 2) / 2)
  expected_spl <- 20 * log10(expected_rms / 20e-6)
  expect_true(all(is.finite(preset$msldB_bio)))
  expect_true(all(abs(preset$msldB_bio - expected_spl) < 0.15))
})
