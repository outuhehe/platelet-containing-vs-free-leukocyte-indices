# Parse every R script to confirm syntax validity after the audit fixes.
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")

dirs <- c(file.path(ra_root, "R"), file.path(ra_root, "tests"), file.path(ra_root, "figures"))
files <- c(Sys.glob(file.path(ra_root, "install_pkgs.R")),
           unlist(lapply(dirs, function(d) Sys.glob(file.path(d, "*.R")))),
           list.files(file.path(ra_root, "exploratory"), pattern = "\\.R$",
                      recursive = TRUE, full.names = TRUE))
files <- unique(files)
bad <- character(0)
for (f in files) {
  res <- tryCatch(parse(f), error = function(e) e)
  if (inherits(res, "error")) {
    bad <- c(bad, paste(basename(f), "->", conditionMessage(res)))
    cat(sprintf("parse FAILED: %s\n", basename(f)))
  } else cat(sprintf("parsed OK: %s\n", basename(f)))
}
if (length(bad) > 0) { cat("SYNTAX ERRORS:\n"); cat(bad, sep = "\n"); stop("parse failures") }
cat("All", length(files), "R scripts parse cleanly.\n")
