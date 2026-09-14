# ============================================================
# ANECO - Reproducible script
# ============================================================

# --------------------------
# 0) Setup
# --------------------------
# install.packages(c("remotes", "ANECO", "ggplot2", "dplyr", "ggeffects",
#                    "progressr", "caret", "rio", "ggsci", "future"))

# Optional: install ANECO from GitHub (only if needed)
remotes::install_github("telesiglab/ANECO")

library(ANECO)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(progressr)
library(caret)
library(future)

set.seed(456)

progressr::handlers("txtprogressbar")
progressr::handlers(global = TRUE)

# --------------------------
# 1) Download & unpack data
# --------------------------
tar_url <- "https://zenodo.org/records/18778160/files/Amistosa_soundscape.tar"

tmp_dir <- tempfile("amistosa_audio_")
dir.create(tmp_dir, recursive = TRUE)
tar_file <- file.path(tmp_dir, "Amistosa_soundscape.tar")

options(timeout = 60 * 20)  # 20 min

download.file(tar_url, tar_file, mode = "wb", quiet = FALSE)
untar(tar_file, exdir = tmp_dir)

# quick sanity check
wav_files <- list.files(tmp_dir, pattern = "\\.wav$", full.names = FALSE)
print(head(wav_files, 5))

# --------------------------
# 2) Run acoustic_indices()
# --------------------------
sel_indices <- c("ACI", "ADI", "AAc", "msldB_bio")

ind_df <- acoustic_indices(
  dir = tmp_dir,
  calibparam = "SMmini",
  sel.ind = sel_indices,       # keep same indices used in the manuscript example
  prefix.format = "SSRRR",
  bioband = c(1000, 10000),
  noiseband = c(0, 1000),
  noise.ind = TRUE,            # keep TRUE (needed later for classifier)
  channel = "left",
  parallel = "files",
  ncores = NULL,               # NULL = all cores - 1 
  wl = 512,
  save.file = FALSE            # avoids writing CSVs in temp folder
)


print(head(ind_df))

# --------------------------
# 3) Forest-type comparison (linear models + Figure 1)
# --------------------------
fit_site_models <- function(df) {
  list(
    ACI       = lm(ACI ~ Site, data = df),
    ADI       = lm(ADI ~ Site, data = df),
    AAc       = lm(AAc ~ Site, data = df),
    msldB_bio = lm(msldB_bio ~ Site, data = df)
  )
}

effects_long <- function(model_list, rain_label = NULL) {
  out <- lapply(names(model_list), function(nm) {
    m <- model_list[[nm]]
    x <- as.data.frame(ggeffects::ggeffect(m))
    x$indice <- nm
    if (!is.null(rain_label)) x$rain <- rain_label
    x
  })
  bind_rows(out)
}

models_all <- fit_site_models(ind_df)
df_indices <- effects_long(models_all)

p_fig1 <- ggplot(df_indices, aes(x = Site.x, y = Site.predicted)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Site.conf.low, ymax = Site.conf.high),
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ indice, scales = "free_y") +
  xlab("") + ylab("Index value") +
  theme_bw() +
  theme(legend.position = "none")

print(p_fig1)

# --------------------------
# 4) Rain classifier training (package example data)
# --------------------------
data("indices_data", package = "ANECO")
data("labeled_data", package = "ANECO")

# (Optional) quick class balance inspection
labeled_data %>%
  group_by(Category) %>%
  summarise(sum(dif))

model_list <- soundclassifier_fit(
  ref = labeled_data,
  ind = indices_data,
  p = 0.80,
  nmax = "even",
  classes = c("Rain", "Clean")
)

print(model_list$info)
print(model_list$model)

# --------------------------
# 5) Internal test performance (hold-out from labeled_data)
# --------------------------
test <- model_list$testdata
pred <- predict(model_list$model, newdata = test[, -c(1, 2)])

cm_internal <- caret::confusionMatrix(data = pred, reference = test$Category, positive = "Rain")
print(cm_internal)

metrics_from_cm <- function(cm) {
  data.frame(
    Accuracy  = unname(cm$overall["Accuracy"]),
    Precision = unname(cm$byClass["Pos Pred Value"]),
    Recall    = unname(cm$byClass["Sensitivity"]),
    F1        = unname(cm$byClass["F1"])
  )
}

perf_internal <- metrics_from_cm(cm_internal)
print(perf_internal)

# --------------------------
# 6) Predict rain over Amistosa example + external evaluation
# --------------------------

df2 <- soundclassifier_predict(data = ind_df, mod = model_list$model)

# GitHub raw CSV URL
manual_labels_url <- "https://raw.githubusercontent.com/telesiglab/ANECO/main/Example_usage/Amistosa_labeled_example.csv"

# download ground-truth data
ground_truth <- read.csv(manual_labels_url, stringsAsFactors = TRUE, sep = ";")

ground_truth$Category <- as.factor(ground_truth$Category)

# Confusion matrix
cm_external <- caret::confusionMatrix(
  data = df2$predicted_class,
  reference = ground_truth$Category,
  positive = "Rain"
)

print(cm_external)

# performance metrics
perf_external <- data.frame(
  Accuracy  = unname(cm_external$overall["Accuracy"]),
  Precision = unname(cm_external$byClass["Pos Pred Value"]),
  Recall    = unname(cm_external$byClass["Sensitivity"]),
  F1        = unname(cm_external$byClass["F1"])
)

print(perf_external)


# --------------------------
# 6.1) Classifier performance figure
# --------------------------

# Order used throughout the figure
class_order  <- c("Clean", "Rain")
metric_order <- c("Accuracy", "Precision", "Recall", "F1")

# ---- Panel A: internal hold-out confusion matrix ----

cm_internal_df <- as.data.frame(cm_internal$table) %>%
  rename(
    Predicted = Prediction,
    Observed  = Reference,
    Count     = Freq
  ) %>%
  mutate(
    Observed = factor(Observed, levels = class_order),
    # Reverse order so the correct classifications appear
    # from the upper-left to the lower-right
    Predicted = factor(Predicted, levels = rev(class_order))
  )

p_confusion <- ggplot(cm_internal_df, aes(x = Observed, y = Predicted, fill = Count)) +
  geom_tile(color = "white",linewidth = 1.2) +
  geom_text(aes(label = Count, color = Count > max(Count) / 2), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#F1F5F9", high = "#3C5488") +
  scale_color_manual( values = c("FALSE" = "#263238", "TRUE" = "white"), guide = "none") +
  labs(x = "Observed class", y = "Predicted class", fill = "Count", tag = "a",) +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")


# ---- Panel B: comparison of performance metrics ----

make_metric_long <- function(metric_df, dataset_name) {
  data.frame(
    Dataset = dataset_name,
    Metric = metric_order,
    Value = as.numeric(metric_df[1, metric_order]),
    stringsAsFactors = FALSE
  )
}

performance_df <- bind_rows(
  make_metric_long(perf_internal, "Internal hold-out"),
  make_metric_long( perf_external, "Amistosa example")
  ) %>%
  mutate(Dataset = factor(Dataset, levels = c("Internal hold-out", "Amistosa example")), Metric = factor(Metric, levels = metric_order))

library(ggplot2)

p_metrics <- ggplot(performance_df, aes(x = Metric, y = Value, color = Dataset, group = Dataset)) +
  geom_linerange( aes(ymin = 0.5, ymax = Value), position = position_dodge(width = 0.78), linewidth = 0.8 ) +
  geom_point(position = position_dodge(width = 0.78), size = 4) +
  geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(width = 0.78), vjust = -1.2, size = 3.5, show.legend = FALSE ) +
    scale_color_manual( values = c(
      "Internal hold-out" = "#3C5488",
      "Amistosa example" = "#00A087")) +
  scale_y_continuous(
    limits = c(0.5, 1.08),
    breaks = seq(0.5, 1, by = 0.2),
    labels = function(x) paste0(round(x * 100), "%"),
    expand = expansion(mult = c(0, 0.01))) +
  labs(x = NULL, y = "Performance metric", color = NULL, tag = "b") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

# Mostrar el gráfico
print(p_metrics)


# ---- Combine both panels ----

p_classifier_performance <- p_confusion + p_metrics 

print(p_classifier_performance)


# ---- Export publication-quality files ----

ggsave(
  filename = "performance_plot.png",
  plot = p_metrics,
  width = 10,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# --------------------------
# 7) Filter out predicted rain minutes and refit models (Figure 2)
# --------------------------
merged_df <- left_join(
  df2,
  ind_df,
  by = c("Name", "Minute", "Site", "RecorderID", "Date", "Hour")
)

filtered_df <- merged_df %>%
  filter(predicted_class != "Rain")

models_no_rain <- fit_site_models(filtered_df)

df_clean <- effects_long(models_no_rain, rain_label = "Without Rain")
df_rain  <- df_indices %>% mutate(rain = "With Rain")

dd <- bind_rows(df_clean, df_rain)

library(ggsci)

p_fig2 <- ggplot(dd, aes(x = rain, y = Site.predicted, colour = Site.x)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Site.conf.low, ymax = Site.conf.high),
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_wrap(~ indice, scales = "free_y") +
  xlab("") + ylab("Index value") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_color_npg()

print(p_fig2)
