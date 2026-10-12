# ANECO

**ANalysis of ECOacoustics: acoustic indices, sound pressure calibration, and heavy-rain screening in R.**

ANECO connects acoustic-index computation, recorder calibration, and supervised sound classification in one workflow for passive acoustic monitoring. It produces minute-level tables for soundscape comparisons and downstream ecological analyses.

Developed by the TELESIG laboratory at the Instituto Internacional en Conservación y Manejo de Vida Silvestre (ICOMVIS), Universidad Nacional, Costa Rica, this repository accompanies **“ANECO: an R package for easy calculation of acoustic indices and identification of heavy rain”** by Jimy Barrantes-Madrigal and Mónica Retamosa-Izaguirre.

## Features

- Acoustic indices covering complexity, diversity, entropy, sound levels, acoustic activity, and spectral structure.
- Integrated sound pressure level (SPL) calibration using recorder presets or microphone sensitivity, gain, and ADC voltage.
- Energy and dominant-frequency summaries in 1 kHz bands.
- Parallel processing across recordings or across minutes within a recording.
- A retrainable Random Forest workflow for identifying minutes dominated by heavy rain.
- Bundled recordings, index tables, labels, and article-related analysis scripts.

ANECO includes its own calibrated calculations and implementations from `seewave` and `soundecology`. Codes such as `ACI`, `ACI_seewave`, and `ACI_soundecology` distinguish these alternatives. The index catalog describes each implementation and its calibration status.

## Installation

```r
install.packages("remotes")
remotes::install_github("telesiglab/ANECO")
library(ANECO)
```

The package declares R >= 3.5.0. Dependencies are listed in [DESCRIPTION](DESCRIPTION) and installed by `remotes` as needed. Compatibility also depends on the R requirements of those dependencies.

## Main functions

| Function | Purpose |
| --- | --- |
| `acoustic_indices()` | Segment WAV recordings and calculate indices and metadata. |
| `aneco_list_indices()` | Inspect index codes, descriptions, references, and implementations. |
| `soundclassifier_fit()` | Train a classifier using labeled intervals and matching indices. |
| `soundclassifier_predict()` | Return class probabilities and predicted labels. |
| `PAMCalib()` | Perform the lower-level calibration used by `acoustic_indices()`. |

## Calculate acoustic indices

Place WAV recordings in a folder and select the indices to calculate:

```r
indices <- acoustic_indices(
  dir = "path/to/wav_folder",
  calibparam = "SM4",
  sel.ind = c("ACI", "ADI", "msldB_bio"),
  noise.ind = FALSE,
  prefix.format = "SSRR",
  channel = "left",
  bioband = c(1000, 11000),
  noiseband = c(0, 1000),
  wl = 512,
  parallel = "files",
  ncores = 1,
  save.file = FALSE
)
head(indices)
```

Inspect the full catalog:

```r
catalog <- aneco_list_indices()
catalog[, c("code", "label", "family", "calibrated", "source")]
```

Use `sel.ind = "all"` to select all available codes; `exclude` can then remove specific codes. Set `noise.ind = TRUE` to add the predictors required by rain classification, even when absent from `sel.ind`.

### File names and output

Use `PREFIX_YYYYMMDD_HHMMSS.wav`. For example, `AR01_20260101_060000.wav` with `prefix.format = "SSRR"` identifies site `AR` and recorder `01`. For `BAJA08_20260101_060000.wav`, use `"SSSSRR"`. Files without a prefix can use `add.prefix` with the corresponding format.

WAV files are read directly inside the specified folder, without recursively searching subfolders. Files are divided into 60-second segments; a final segment is included if it lasts at least 30 seconds.

The returned `data.frame` contains requested indices and identification fields: `Name`, `Code`, `Site`, `RecorderID`, `Date`, `Year`, `Month`, `Day`, `Hour`, `Minute`, and `Duration`. `Minute` starts at zero within each file; `Duration` is in seconds. `Hour` retains the recording's start-time string, so use it with `Minute` when reconstructing segment times.

With `save.file = TRUE`, results are written to a CSV in the audio folder. If no valid minutes remain, an empty table is returned without creating a CSV.

### Calibration, frequency bands, and parallel processing

Recorder presets are `"SM4"`, `"SMmini"`, and `"SM2+"`. Alternatively, supply `c(micsens, gain, vADC)`: microphone sensitivity in dB re 1 V/Pa, gain in dB, and ADC zero-to-peak voltage in volts. For example, `calibparam = c(-35, 48, 1)` supplies the SM4 preset values. Match parameters to the microphone, recorder, and actual recording settings. Calibration does not remove all differences in microphone response, condition, or device noise.

`bioband` and `noiseband` are in Hz. The highest available frequency is half the sample rate. With `efb`, results from different sample rates are combined by column name; unavailable bands become `NA`, not zero. The classifier requires complete energy predictors through `E1011` (10–11 kHz), so choose compatible recordings before training or prediction.

Use `parallel = "files"` across recordings or `parallel = "minutes"` within recordings. Start with `ncores = 1` and increase according to available resources. With `ncores = NULL`, ANECO selects physical cores minus one, with a minimum of one worker.

Optional progress reporting:

```r
progressr::handlers("txtprogressbar")
progressr::handlers(global = TRUE)
```

### Try the bundled recordings

This example requires no external audio download:

```r
data("example_wavs", package = "ANECO")
example_dir <- tempfile("aneco_example_")
dir.create(example_dir)
tuneR::writeWave(example_wavs$morning,
                 file.path(example_dir, "AR01_20260101_060000.wav"))
tuneR::writeWave(example_wavs$noon,
                 file.path(example_dir, "AR01_20260101_120000.wav"))

example_indices <- acoustic_indices(
  dir = example_dir,
  calibparam = "SM4",
  sel.ind = c("ACI", "ADI", "msldB_bio"),
  noise.ind = FALSE,
  ncores = 1,
  save.file = FALSE
)
head(example_indices)
```

## Train and apply a rain classifier

The classifier targets minutes dominated by rain that masks other sounds. `"Clean"` denotes the comparison class, not necessarily the absence of rainfall. Additional classes can be supplied, but the predictors and intended use are optimized for heavy-rain screening.

This self-contained example uses bundled index tables and labels. The sample cap keeps the demonstration small; `nmax = "even"` balances classes to the smallest matched class, while `NULL` retains all matched samples.

```r
data("indices_data", package = "ANECO")
data("labeled_data", package = "ANECO")
data("newdata", package = "ANECO")
set.seed(456)
fit <- soundclassifier_fit(
  ref = labeled_data,
  ind = indices_data,
  classes = c("Clean", "Rain"),
  p = 0.8,
  nmax = 100
)
fit$info

# Evaluate on held-out observations before applying the model.
test_predictors <- fit$testdata[, setdiff(names(fit$testdata),
                                         c("Name", "Category")), drop = FALSE]
test_predictions <- predict(fit$model, newdata = test_predictors, type = "raw")
caret::confusionMatrix(test_predictions, fit$testdata$Category, positive = "Rain")

predictions <- soundclassifier_predict(data = newdata, mod = fit$model)
head(predictions)

# Prediction preserves input row order; retain the original indices when filtering.
without_rain <- newdata[predictions$predicted_class != "Rain", , drop = FALSE]
```

Training returns `model`, per-class training counts in `info`, and held-out observations in `testdata`. Prediction returns identification columns, one `prob_*` column per learned class, and `predicted_class`.

### Train with your own annotations

Calculate indices with `noise.ind = TRUE` and supply a label table with these fields:

| Column | Meaning |
| --- | --- |
| `File` | Recording filename matching `Name` in the index table. |
| `Category` | Class label, such as `Clean` or `Rain`. |
| `Start` | Start minute, counted from zero within the recording. |
| `End` | Exclusive end minute. |

For example, `Start = 2` and `End = 5` label minutes 2, 3, and 4. Use integer minute boundaries and resolve conflicting or overlapping annotations before training. Replace the bundled labels and index table with your matched local data.

Use the same index implementations, calibration settings, frequency bands, and package revision for training and prediction. When calibration or spectral calculations change, regenerate training predictors and retrain rather than mixing processing versions.

Validate with representative local recordings. Inspect false `Rain` predictions, which remove usable recordings, and false `Clean` predictions, which retain rain-dominated recordings. Retrain when recorder types, habitats, or noise conditions change.

## Article data and analysis scripts

The article demonstrates the workflow using 180 one-minute observations from mature and secondary tropical forests in Costa Rica. The [Amistosa_soundscape dataset is available on Zenodo](https://zenodo.org/records/18778160).

| Resource | Contents |
| --- | --- |
| [Amistosa_example.R](Example_usage/Amistosa_example.R) | Data download, indices, habitat comparisons, rain classification, evaluation, and filtering. |
| [Amistosa_labeled_example.csv](Example_usage/Amistosa_labeled_example.csv) | Manual labels for the Amistosa example. |
| [Benchmark ANECO vs soundecology.R](Example_usage/Benchmark%20ANECO%20vs%20soundecology.R) | Runtime comparisons across CPU counts. |
| [scikit-maad_bench.ipynb](Example_usage/scikit-maad_bench.ipynb) | Python benchmark workflow. |
| [Calibracion.R](Example_usage/Calibracion.R) | Cross-recorder calibration comparisons and bootstrap summaries. |
| [CalibrationTest/](Example_usage/CalibrationTest/) | Calibration-analysis tables included in the repository. |

These scripts have additional analysis dependencies and may require adjustments to local paths. The raw recordings for the calibration assessment are not distributed here; availability is described in the calibration script.

For article reproduction, retain the code revision, input tables, random seed, and parameters used for that analysis. The article used ANECO 1.0 and R 4.4.2. Results from later revisions can differ because index implementations and calibration calculations have been revised.

Acoustic indices describe soundscape structure and require ecological validation. The article's classifier performance and benchmarks refer to its evaluated datasets and settings; they do not establish universal model accuracy or runtime.

## Repository layout and development

| Directory | Contents |
| --- | --- |
| [R/](R/) | Package functions. |
| [man/](man/) | Function and dataset documentation. |
| [data/](data/) | Bundled data and index catalog. |
| [tests/testthat/](tests/testthat/) | Automated tests. |
| [Example_usage/](Example_usage/) | Article analyses and benchmarks. |

From a local checkout, install the test tools and run the tests:

```r
install.packages(c("testthat", "pkgload"))
testthat::test_local(".")
```

Use `?acoustic_indices`, `?soundclassifier_fit`, and `?soundclassifier_predict` for function documentation. Report bugs or suggestions through [GitHub Issues](https://github.com/telesiglab/ANECO/issues), including the package version, parameters, and a reproducible example.

## Citation

When using ANECO in research, cite the associated article and record the package version or commit used:

> Barrantes-Madrigal, J., & Retamosa-Izaguirre, M. (2027). ANECO: an R package for easy calculation of acoustic indices and identification of heavy rain. *Ecological Frontiers*. https://doi.org/10.3724/j.ecofro-20260345

## License and contact

The software is distributed under the GNU General Public License, version 3; see [LICENSE](LICENSE). Publication and external dataset reuse are governed by their respective terms.

Maintainer: Jimmy Barrantes Madrigal, [jimmybarrantesm@gmail.com](mailto:jimmybarrantesm@gmail.com).
