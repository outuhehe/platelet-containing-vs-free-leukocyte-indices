# Unit tests for pgpass parsing/matching without touching the real file.
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")
ra_root <- normalizePath(file.path(dirname(this_file), ".."), winslash = "/")
Sys.setenv(RA_ROOT = ra_root)
suppressMessages(source(file.path(ra_root, "R", "00_common.R")))

fails <- 0L
check <- function(name, cond) {
  if (isTRUE(cond)) cat("PASS:", name, "\n") else { cat("FAIL:", name, "\n"); fails <<- fails + 1L }
}

# 1. Standard five-field line
r <- parse_pgpass_line("127.0.0.1:5431:mimic:alice:secret")
check("standard line has 5 fields", length(r) == 5L)
check("standard user/password", r[[4]] == "alice" && r[[5]] == "secret")

# 2. Password containing a colon is kept intact
r2 <- parse_pgpass_line("h:5432:db:bob:p:a:s:s")
check("colon inside password preserved", r2[[5]] == "p:a:s:s" && length(r2) == 5L)

# 3. Backslash escape (escaped colon and backslash)
r3 <- parse_pgpass_line("h:5432:db:bob:pa\\:ss\\\\word")
check("backslash escapes decoded", r3[[5]] == "pa:ss\\word")

# 4. Wildcard field matching
check("wildcard matches any value", pg_field_match("*", "anything"))
check("exact match", pg_field_match("mimic", "mimic"))
check("non-match returns false", !pg_field_match("mimic", "eicu"))

# 5. Multi-record selection against the production lookup_pgpass().
# Fix deterministic connection variables for the test and restore them after.
saved <- list(PGHOST = PGHOST, PGPORT = PGPORT, MIMIC_DB = MIMIC_DB,
              EICU_DB = EICU_DB, PGUSER = PGUSER)
PGHOST <- "127.0.0.1"; PGPORT <- 5431
MIMIC_DB <- "mimic"; EICU_DB <- "eicu"; PGUSER <- ""

# 5a. First non-matching host is skipped; wildcard host then matches.
tmp <- tempfile()
writeLines(c(
  "10.0.0.9:5431:mimic:wrong:no",
  "*:5431:mimic:carol:wildpw"
), tmp)
pick <- lookup_pgpass(tmp)
check("skips bad host and matches wildcard record", !is.null(pick) &&
        pick$user == "carol" && pick$password == "wildpw")

# 5b. First matching record is returned in file order.
tmp2 <- tempfile()
writeLines(c(
  "127.0.0.1:5431:mimic:exact:exactpw",
  "*:5431:mimic:carol:wildpw"
), tmp2)
pick2 <- lookup_pgpass(tmp2)
check("first matching record returned in order", !is.null(pick2) &&
        pick2$user == "exact" && pick2$password == "exactpw")

# 5c. A database that is neither target db is rejected.
tmp3 <- tempfile()
writeLines("127.0.0.1:5431:otherdb:any:no", tmp3)
check("non-target database rejected", is.null(lookup_pgpass(tmp3)))

# 5d. With PGUSER set, a non-matching user is skipped.
PGUSER <- "alice"
tmp4 <- tempfile()
writeLines(c(
  "127.0.0.1:5431:mimic:bob:wronguser",
  "127.0.0.1:5431:mimic:alice:correctuser"
), tmp4)
pick4 <- lookup_pgpass(tmp4)
check("user mismatch rejected and matching user selected", !is.null(pick4) &&
        pick4$user == "alice" && pick4$password == "correctuser")

# 5e. A missing file returns NULL.
check("missing pgpass file returns NULL", is.null(lookup_pgpass(tempfile())))

# Restore the original connection variables and clean up.
PGHOST <- saved$PGHOST; PGPORT <- saved$PGPORT; MIMIC_DB <- saved$MIMIC_DB
EICU_DB <- saved$EICU_DB; PGUSER <- saved$PGUSER
unlink(c(tmp, tmp2, tmp3, tmp4))

# 6. Malformed line (<5 fields) yields NULL
check("malformed line returns NULL", is.null(parse_pgpass_line("a:b:c")))

if (fails > 0L) stop(fails, " pgpass test(s) failed")
cat("All pgpass tests passed.\n")
