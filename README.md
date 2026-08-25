# ANECO

ANECO (Analysis of ECOacoustics) is an R package for:

- Computation of ~30 acoustic indices
- Integrated SPL calibration
- Random Forest classification of rain-dominated segments
- Parallel processing support

## Installation

```r
# install.packages("remotes")
remotes::install_github("telesiglab/ANECO")

library(ANECO)

indices <- acoustic_indices(
  dir = "path/to/wavs",
  calibparam = "SM4",
  noise.ind = TRUE
)
