Sys.setenv(FRONTMED_SKIP_MAIN = "1")
args <- commandArgs(FALSE)
this <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]), winslash = "/")
source(file.path(dirname(this), "..", "R", "21_reclassification_winsorization.R"))
y <- c(0L, 0L, 1L, 1L)
old <- rep(.5, 4); better <- c(.1, .2, .8, .9)
zero <- reclassification_metrics(y, old, old)
good <- reclassification_metrics(y, old, better)
stopifnot(all(zero == 0), good[["NRI_continuous"]] == 2,
          abs(good[["IDI"]] - .7) < 1e-12,
          all(abs(good + reclassification_metrics(y, better, old)) < 1e-12),
          auc_rank(y, better) == 1, auc_rank(y, 1 - better) == 0,
          auc_rank(y, old) == .5)
set.seed(9)
subject <- rep(seq_len(60), each = 2)
yy <- rep(rep(0:1, each = 30), each = 2)
f <- make_folds(yy, subject)
stopifnot(all(vapply(split(f, subject), function(x) length(unique(x)), integer(1)) == 1L))
x <- c(1:99, 10000)
stopifnot(length(winsorize(x, .01)) == length(x), max(winsorize(x, .01)) < max(x),
          identical(winsorize(x, 0), x))
set.seed(31)
d <- data.frame(death = rbinom(300, 1, .15), age = rnorm(300, 60, 10),
  sex = factor(sample(c("F", "M"), 300, TRUE)),
  adm_prediction = factor(sample(c("emergency", "elective", "other"), 300, TRUE)),
  charlson = sample(0:8, 300, TRUE))
X <- model.matrix(~ age + sex + adm_prediction + charlson, d)
fast <- glm.fit(X, d$death, family = binomial())
ref <- glm(death ~ age + sex + adm_prediction + charlson, data = d, family = binomial())
stopifnot(max(abs(plogis(drop(X %*% fast$coefficients)) - predict(ref, type = "response"))) < 1e-12)
cat("PASS: NRI/IDI identity/direction, AUC ties/direction, grouped folds, winsorization, glm.fit equivalence\n")
