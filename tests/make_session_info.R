# Generate session_info.txt: R version and exact versions of every package the
# pipeline uses. Mirrors the reporting practice in Ding (R 4.4.0 + named pkgs)
# and Jin (R 4.3.3 + per-package versions). No database connection required.
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
suppressMessages(source(file.path(ra_root, "R", "00_common.R")))

# Load every analysis package so sessionInfo reflects the full run environment.
load_analysis_pkgs("pROC", "survRM2")

# Map absolute library paths to shareable logical labels so the report does not
# expose the local user directory.
logical_lib <- function(p) {
  p <- normalizePath(p, winslash = "/", mustWork = FALSE)
  root_n <- normalizePath(ra_root, winslash = "/")
  if (startsWith(p, root_n)) {
    rel <- sub("^/", "", substr(p, nchar(root_n) + 1L, nchar(p)))
    return(paste0("project-local library: ", rel))
  }
  if (grepl("win-library", p) || grepl("site-library", p))
    return("user library (R default user library)")
  if (grepl("Program Files", p) || grepl("/library$", p))
    return("R system library (R installation)")
  "other library"
}

# Resolve which entry of .libPaths() a package is loaded from, robust to both
# versioned user libraries (win-library/4.5/<pkg>) and the project Rlib.
lib_for_package <- function(pkg) {
  pkg_path <- normalizePath(system.file(package = pkg), winslash = "/")
  libs <- normalizePath(.libPaths(), winslash = "/")
  hit <- libs[startsWith(pkg_path, paste0(libs, "/"))]
  if (!length(hit)) return(NA_character_)
  hit[which.max(nchar(hit))]
}

pkgs <- c("DBI", "RPostgres", "data.table", "pROC", "survRM2")
out <- c(
  "Session information for the reproducible analysis",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## R and platform",
  capture.output(R.version.string),
  paste("Platform:", R.version$platform),
  paste("Architecture:", R.version$arch),
  paste("R OS string:", R.version$os),
  "(Exact operating system is reported in the sessionInfo() section below.)",
  "",
  "## Library paths (logical labels; absolute local paths are not recorded)",
  vapply(.libPaths(), logical_lib, character(1), USE.NAMES = FALSE),
  "",
  "## Exact package versions used by the pipeline"
)
for (p in pkgs) {
  loc <- logical_lib(lib_for_package(p))
  out <- c(out, sprintf("%-12s %s   [%s]", p, as.character(packageVersion(p)), loc))
}
out <- c(out, "", "## Full sessionInfo()", capture.output(sessionInfo()))

writeLines(out, file.path(ra_root, "session_info.txt"))
cat("session_info.txt written\n")
cat(paste(out[seq_len(length(out) - length(capture.output(sessionInfo())) - 1)],
          collapse = "\n"), "\n")
