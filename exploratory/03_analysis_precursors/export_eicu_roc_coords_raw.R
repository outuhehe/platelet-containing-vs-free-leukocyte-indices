# Load shared paths/helpers; pROC is added for the ROC calculations.
.this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
.ra_root <- normalizePath(file.path(dirname(.this_file), "..", ".."), winslash = "/")
Sys.setenv(RA_ROOT = .ra_root)
suppressMessages(source(file.path(.ra_root, "R", "00_common.R")))
# Exploration can generate row-level records and observed ROC thresholds.
# Keep every exploratory artifact inside the restricted, git-ignored tree.
P_OUT <- file.path(P_INTERIM, "exploratory")
dir.create(P_OUT, recursive = TRUE, showWarnings = FALSE)
load_analysis_pkgs("pROC")

in_path <- file.path(P_TABLES, "eicu_sicu_external_validation_dataset.csv")
df <- fread(in_path)
if (nrow(df) != 2855L || sum(df$death, na.rm = TRUE) != 276L) {
  stop("Unexpected eICU input size or event count")
}

metric_cols <- c(NLR = "NLR", PLR = "PLR", SII = "SII", SIRI = "SIRI",
                PIV = "PIV", ANC = "neutrophil", ALC = "lymphocyte")
missing <- setdiff(unname(metric_cols), names(df))
if (length(missing)) stop("Missing eICU fields: ", paste(missing, collapse = ", "))

out <- rbindlist(lapply(names(metric_cols), function(indicator) {
  x <- df[[metric_cols[[indicator]]]]
  ok <- is.finite(x) & !is.na(df$death)
  if (any(x[ok] < 0)) stop("Negative values in ", indicator)
  direction <- if (indicator == "ALC") ">" else "<"
  roc_obj <- pROC::roc(
    response = factor(df$death[ok], levels = c(0, 1)), predictor = x[ok],
    levels = c("0", "1"), direction = direction, quiet = TRUE
  )
  co <- as.data.table(pROC::coords(
    roc_obj, "all", ret = c("specificity", "sensitivity", "threshold"),
    transpose = FALSE
  ))
  data.table(
    database = "eICU-CRD", indicator = indicator,
    fpr = 1 - co$specificity, tpr = co$sensitivity,
    threshold = co$threshold, auc = as.numeric(pROC::auc(roc_obj)),
    n = sum(ok), events = sum(df$death[ok] == 1L)
  )
}))

fwrite(out, file.path(P_OUT, "roc_coords_eICU.csv"))
mimic <- fread(file.path(P_OUT, "roc_coords_MIMIC.csv"))
fwrite(rbindlist(list(mimic, out), use.names = TRUE, fill = TRUE),
       file.path(P_OUT, "roc_coords_both.csv"))

cat("Raw-value eICU coordinates; high-risk markers direction='<', ALC direction='>'.\n")
print(unique(out[, .(indicator, auc, n, events)]))
