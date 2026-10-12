test_that("efb combines mixed sample rates by column name", {
  set.seed(654)
  mixed_dir <- withr::local_tempdir()
  separate_dirs <- c(withr::local_tempdir(), withr::local_tempdir())
  rates <- c(16000, 22050)
  filenames <- c("AR01_20260101_000000.wav", "AR02_20260101_000000.wav")
  for (i in seq_along(rates)) {
    wave <- tuneR::Wave(
      left = as.integer(round(rnorm(90 * rates[i], sd = 200))),
      samp.rate = rates[i], bit = 16
    )
    file <- file.path(separate_dirs[i], filenames[i])
    tuneR::writeWave(wave, file)
    file.copy(file, file.path(mixed_dir, filenames[i]))
  }
  writeLines("not a WAV", file.path(mixed_dir, "ZZ01_20260101_000000.wav"))
  compute <- function(dir, mode = "files", save = FALSE) {
    acoustic_indices(dir = dir, sel.ind = "efb", noise.ind = FALSE,
                     ncores = 1, parallel = mode, save.file = save)
  }
  separate <- lapply(separate_dirs, compute)
  upper_bands <- c("E89", "E910", "E1011")
  common_bands <- paste0("E", 0:7, 1:8)
  for (mode in c("files", "minutes")) {
    expect_warning(result <- compute(mixed_dir, mode, save = mode == "files"),
                   "ZZ01_20260101_000000.wav")
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 4L)
    expect_true(all(c(common_bands, upper_bands) %in% names(result)))
    low <- result[result$Name == filenames[1], , drop = FALSE]
    high <- result[result$Name == filenames[2], , drop = FALSE]
    expect_equal(nrow(low), 2L)
    expect_equal(nrow(high), 2L)
    expect_true(all(is.na(as.matrix(low[, upper_bands, drop = FALSE]))))
    expect_true(all(is.finite(as.matrix(high[, upper_bands, drop = FALSE]))))
    expect_true(all(is.finite(as.matrix(low[, common_bands, drop = FALSE]))))
    expect_false("freq.energy" %in% names(result))
    for (i in seq_along(filenames)) {
      actual <- result[result$Name == filenames[i], names(separate[[i]]), drop = FALSE]
      rownames(actual) <- NULL
      expect_equal(actual, separate[[i]], tolerance = 1e-10)
    }
  }
  csv <- list.files(mixed_dir, pattern = "[.]csv$", full.names = TRUE)
  expect_length(csv, 1L)
  saved <- read.csv(csv)
  expect_equal(nrow(saved), 4L)
  expect_true(all(is.na(as.matrix(saved[saved$Name == filenames[1], upper_bands]))))
})
