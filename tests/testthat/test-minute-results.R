test_that("a single valid minute retains its columns and can be saved", {
  set.seed(789)
  audio_dir <- withr::local_tempdir()
  rate <- 8192
  # The second file contributes no rows; values exist only in the first row.
  for (duration in c(45, 20)) {
    samples <- as.integer(round(rnorm(duration * rate, sd = 200)))
    name <- sprintf("AR01_20260101_%06d.wav", if (duration == 45) 0 else 100)
    tuneR::writeWave(tuneR::Wave(left = samples, samp.rate = rate, bit = 16),
                     file.path(audio_dir, name))
  }
  result <- acoustic_indices(
    dir = audio_dir, sel.ind = c("msldB_bio", "efb"), noise.ind = FALSE,
    bioband = c(1000, 3000), ncores = 1, save.file = TRUE
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_equal(result$Name, "AR01_20260101_000000.wav")
  expect_equal(result$Minute, 0)
  expect_equal(result$Duration, 45)
  expect_true(all(c("Name", "Code", "Site", "RecorderID", "Date", "Year",
                    "Month", "Day", "Hour", "Minute", "Duration",
                    "msldB_bio", "E01") %in% names(result)))
  expect_true(is.finite(result$msldB_bio))
  csv <- list.files(audio_dir, pattern = "[.]csv$", full.names = TRUE)
  expect_length(csv, 1L)
  saved <- read.csv(csv)
  expect_equal(nrow(saved), 1L)
  expect_equal(names(saved), names(result))
  expect_equal(saved$msldB_bio, result$msldB_bio, tolerance = 1e-8)
})

test_that("an unreadable WAV returns zero rows without writing a CSV", {
  audio_dir <- withr::local_tempdir()
  writeLines("not a WAV", file.path(audio_dir, "AR01_20260101_000000.wav"))
  expect_warning(result <- acoustic_indices(dir = audio_dir, ncores = 1),
                 "following files produced errors")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
  expect_true(all(c("Name", "Code", "Minute") %in% names(result)))
  expect_length(list.files(audio_dir, pattern = "[.]csv$"), 0L)
})

test_that("a WAV shorter than 30 seconds returns an empty data.frame", {
  audio_dir <- withr::local_tempdir()
  wave <- tuneR::Wave(left = as.integer(rep(100, 20 * 8192)),
                     samp.rate = 8192, bit = 16)
  tuneR::writeWave(wave, file.path(audio_dir, "AR01_20260101_000000.wav"))
  expect_warning(result <- acoustic_indices(dir = audio_dir, ncores = 1),
                 "No valid minutes")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
  expect_length(list.files(audio_dir, pattern = "[.]csv$"), 0L)
})

test_that("a directory without WAVs returns an empty data.frame", {
  audio_dir <- withr::local_tempdir()
  expect_warning(result <- acoustic_indices(dir = audio_dir, ncores = 1),
                 "No WAV files")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0L)
  expect_length(list.files(audio_dir, pattern = "[.]csv$"), 0L)
})
