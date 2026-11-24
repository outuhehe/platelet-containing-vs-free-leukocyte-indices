# =============================================================================
# 00_common.R -- Single source of truth for packages, config, paths and DB access.
# Every numbered script starts with:  source("../R/00_common.R")  (or via run_all)
# Design follows MIT-LCP/mimic-code (config separated from code, read-only access).
# =============================================================================

# ---- locate the package root -------------------------------------------------
# Prefer RA_ROOT set by run_all.ps1; otherwise derive from this file's location.
.this_file <- function() {
  args <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(args)) normalizePath(args[1], winslash = "/", mustWork = FALSE) else NA_character_
}
RA_ROOT <- Sys.getenv("RA_ROOT", unset = NA_character_)
if (is.na(RA_ROOT) || RA_ROOT == "") {
  f <- .this_file()
  RA_ROOT <- if (!is.na(f)) normalizePath(file.path(dirname(f), ".."), winslash = "/") else getwd()
}
RA_ROOT <- normalizePath(RA_ROOT, winslash = "/")

# Project-local library first (self-contained deps such as survRM2), then global.
P_RLIB <- file.path(RA_ROOT, "Rlib")
if (dir.exists(P_RLIB)) .libPaths(c(P_RLIB, .libPaths()))

suppressPackageStartupMessages({
  library(DBI)
  library(RPostgres)
  library(data.table)
})

# Load analysis-only packages on demand (not every script needs pROC/survRM2).
load_analysis_pkgs <- function(...) {
  need <- c(...)
  suppressPackageStartupMessages(for (p in need) library(p, character.only = TRUE))
}

# Standard paths (replaces scattered "outputs/..." relative strings)
P_CFG     <- file.path(RA_ROOT, "config")
P_SQL     <- file.path(RA_ROOT, "sql")
P_R       <- file.path(RA_ROOT, "R")
P_TESTS   <- file.path(RA_ROOT, "tests")
P_OUT     <- file.path(RA_ROOT, "outputs")
P_INTERIM <- file.path(P_OUT, "interim")
P_TABLES  <- file.path(P_OUT, "tables")
P_FIGURES <- file.path(P_OUT, "figures")
P_LOGS    <- file.path(P_OUT, "logs")
for (p in c(P_INTERIM, P_TABLES, P_FIGURES, P_LOGS)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

# ---- load config/db.env ------------------------------------------------------
read_env_file <- function(path) {
  if (!file.exists(path)) return(list())
  ln <- readLines(path, warn = FALSE, encoding = "UTF-8")
  ln <- trimws(ln[!grepl("^\\s*#", ln) & grepl("=", ln)])
  out <- list()
  for (z in ln) {
    k <- trimws(sub("=.*$", "", z)); v <- trimws(sub("^[^=]*=", "", z))
    out[[k]] <- v
  }
  out
}
cfg <- read_env_file(file.path(P_CFG, "db.env"))
getcfg <- function(key, default = "") {
  v <- cfg[[key]]
  if (is.null(v) || is.na(v) || v == "") Sys.getenv(key, unset = default) else v
}

PGHOST <- getcfg("PGHOST", "127.0.0.1")
PGPORT <- as.integer(getcfg("PGPORT", "5431"))
MIMIC_DB <- getcfg("MIMIC_DB", "mimic")
EICU_DB  <- getcfg("EICU_DB", "eicu")
PG_READONLY_OPTIONS <- getcfg("PG_READONLY_OPTIONS",
                              "-c default_transaction_read_only=on -c statement_timeout=0")
PG_WORK_MEM <- getcfg("PG_WORK_MEM", "1GB")
PGUSER <- getcfg("PGUSER", "")
PGPASSWORD <- getcfg("PGPASSWORD", "")

# Fall back to pgpass.conf when credentials are not supplied in db.env.
# Parse every record and match host/port/database/user with '*' wildcards,
# honoring backslash escapes and a ':' that may occur inside the password.
parse_pgpass_line <- function(line) {
  cur <- ""; fields <- character(0); sep_seen <- 0L; i <- 1L; n <- nchar(line)
  while (i <= n) {
    ch <- substr(line, i, i)
    if (ch == "\\" && i < n) {
      cur <- paste0(cur, substr(line, i + 1L, i + 1L)); i <- i + 2L; next
    }
    if (ch == ":" && sep_seen < 4L) {
      fields <- c(fields, cur); cur <- ""; sep_seen <- sep_seen + 1L; i <- i + 1L; next
    }
    cur <- paste0(cur, ch); i <- i + 1L
  }
  fields <- c(fields, cur)   # 5th field = password (may itself contain ':')
  if (length(fields) < 5L) NULL else fields
}
pg_field_match <- function(pattern, value) pattern == "*" || pattern == value
lookup_pgpass <- function(path = NULL) {
  # 'path' is injectable for unit tests; otherwise it defaults to the real file.
  if (is.null(path)) path <- file.path(Sys.getenv("APPDATA"), "postgresql", "pgpass.conf")
  if (!file.exists(path)) return(NULL)
  pgpass <- path
  lines <- readLines(pgpass, warn = FALSE)
  lines <- trimws(lines[nzchar(trimws(lines)) & !grepl("^#", trimws(lines))])
  for (ln in lines) {
    rec <- parse_pgpass_line(ln)
    if (is.null(rec)) next
    ok_host <- pg_field_match(rec[[1]], PGHOST)
    ok_port <- pg_field_match(rec[[2]], as.character(PGPORT))
    ok_db   <- rec[[3]] == "*" || rec[[3]] %in% c(MIMIC_DB, EICU_DB)
    ok_user <- rec[[4]] == "*" || (PGUSER != "" && rec[[4]] == PGUSER) || PGUSER == ""
    if (ok_host && ok_port && ok_db && ok_user)
      return(list(user = rec[[4]], password = rec[[5]]))
  }
  NULL
}
if (PGUSER == "" || PGPASSWORD == "") {
  pgp <- lookup_pgpass()
  if (!is.null(pgp)) {
    if (PGUSER == "" && pgp$user != "*") PGUSER <- pgp$user
    if (PGPASSWORD == "") PGPASSWORD <- pgp$password
  }
}

# ---- DB connection factory (read-only, verified) -----------------------------
.connect <- function(dbname) {
  con <- dbConnect(
    Postgres(),
    host = PGHOST, port = PGPORT, dbname = dbname,
    user = PGUSER, password = PGPASSWORD,
    options = PG_READONLY_OPTIONS,
    bigint = "numeric"   # counts/int8 returned as plain numeric, avoids integer64 print traps
  )
  ro <- dbGetQuery(con, "SHOW default_transaction_read_only")[[1]][[1]]
  if (!isTRUE(ro == "on")) {
    dbDisconnect(con); stop("Connection to ", dbname, " did not enter read-only mode")
  }
  con
}
connect_mimic <- function() .connect(MIMIC_DB)
connect_eicu  <- function() {
  con <- .connect(EICU_DB)
  dbExecute(con, paste("SET work_mem =", shQuote(PG_WORK_MEM)))
  con
}

# ---- IO helpers (data.table), routed to standard folders ---------------------
read_interim <- function(name, ...) {
  fread(file.path(P_INTERIM, if (grepl("\\.", name)) name else paste0(name, ".csv")), ...)
}
write_interim <- function(x, name, ...) {
  fwrite(x, file.path(P_INTERIM, paste0(name, ".csv")), ...)
}
read_table <- function(name, ...) {
  fread(file.path(P_TABLES, if (grepl("\\.", name)) name else paste0(name, ".csv")), ...)
}
write_table <- function(x, name, ...) {
  fwrite(x, file.path(P_TABLES, paste0(name, ".csv")), ...)
}
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(...))
  cat(msg, "\n")
  msg
}
