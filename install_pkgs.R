# =============================================================================
# install_pkgs.R -- Install the R packages the pipeline needs, pinned to the
# exact versions recorded in session_info.txt. Idempotent: a package already
# installed at the pinned version is left untouched. A local Rlib/ is used first
# if present (offline); otherwise every package, including survRM2, is installed
# at the pinned version from CRAN.
# Run once on a fresh machine:
#   Rscript install_pkgs.R
# =============================================================================
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file)), winslash = "/")

# Project-local library first so the bundled survRM2 is discoverable.
rlib <- file.path(ra_root, "Rlib")
if (dir.exists(rlib)) .libPaths(c(rlib, .libPaths()))

# A reachable CRAN mirror; override with options(repos=...) if needed.
if (is.null(getOption("repos")) || identical(getOption("repos"), c(CRAN = "@CRAN@"))) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Pinned CRAN versions (keep in sync with session_info.txt). Versions are
# compared as package_version objects, so equivalent spellings such as
# "3.8-3" and "3.8.3" are treated as the same version rather than re-installed.
cran_pinned <- c(
  DBI        = "1.3.0",
  RPostgres  = "1.4.10",
  data.table = "1.18.4",
  pROC       = "1.19.0.1",
  survival   = "3.8-3",
  survRM2    = "1.0-4"
)
pinned_pv <- lapply(cran_pinned, package_version)

# Nothing to install means no network access at all: return before touching CRAN.
need <- character(0)
for (p in names(cran_pinned)) {
  if (!ver_satisfies(installed_pv(p), pinned_pv[[p]])) need <- c(need, p)
}
if (length(need) == 0) {
  cat("All required R packages are present at the pinned versions; no network access needed.\n")
  status <- vapply(names(cran_pinned), function(p) {
    sprintf("  %-12s %-9s %s", p, as.character(packageVersion(p)),
            dirname(dirname(system.file(package = p))))
  }, character(1))
  cat("Package status:\n"); cat(paste(status, collapse = "\n"), "\n")
  quit(save = "no", status = 0L)
}

# remotes::install_version drives the pinned install; ensure remotes exists.
if (!requireNamespace("remotes", quietly = TRUE)) {
  cat("Installing remotes (required for version pinning)...\n")
  install.packages("remotes", quiet = TRUE)
}

installed_pv <- function(p) {
  if (requireNamespace(p, quietly = TRUE)) packageVersion(p) else NULL
}
ver_satisfies <- function(have, want) !is.null(have) && isTRUE(have == want)

for (p in names(cran_pinned)) {
  want <- pinned_pv[[p]]
  have <- installed_pv(p)
  if (ver_satisfies(have, want)) {
    cat(sprintf("  %-11s %s already at pinned version\n", p, cran_pinned[[p]]))
  } else {
    cat(sprintf("  %-11s installing pinned %s (current: %s)\n", p, cran_pinned[[p]],
                if (is.null(have)) "absent" else as.character(have)))
    remotes::install_version(p, version = cran_pinned[[p]],
                             upgrade = "never", quiet = TRUE)
  }
}

# Verify every required package at its pinned version.
required <- names(cran_pinned)
status <- vapply(required, function(p) {
  ok <- requireNamespace(p, quietly = TRUE)
  sprintf("  %-12s %-9s %s", p, as.character(if (ok) packageVersion(p) else "MISSING"),
          if (ok) dirname(dirname(system.file(package = p))) else "")
}, character(1))
cat("Package status:\n"); cat(paste(status, collapse = "\n"), "\n")

mismatch <- character(0)
for (p in names(cran_pinned)) {
  if (!ver_satisfies(installed_pv(p), pinned_pv[[p]]))
    mismatch <- c(mismatch, paste0(p, "!=", cran_pinned[[p]]))
}
if (length(mismatch) > 0) stop("Version pinning not satisfied: ", paste(mismatch, collapse = ", "))
cat("All required R packages are present at the pinned versions.\n")
