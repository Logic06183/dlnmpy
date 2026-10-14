################################################################################
# Reference values for examples/statmed2014.py
#
#   Gasparrini A. Modeling exposure-lag-response associations with distributed
#   lag non-linear models. Statistics in Medicine 2014;33(5):881-899.
#   doi:10.1002/sim.5963
#
# Runs the author's own scripts 01.prep.R to 06.additional.R, unchanged, from a
# clone of https://github.com/gasparrini/2014_gasparrini_StatMed_Rcodedata at
# the pinned commit below, then writes every quantity worth comparing to
# tests/reference/statmed2014.json. Only their pdf() calls are redirected, to a
# temporary directory. The simulation study (07.simcall.R, 08.simres.R) is not
# run: it is random by design and takes hours. About a minute.
#
# Usage, from the repository root:
#
#   git clone https://github.com/gasparrini/2014_gasparrini_StatMed_Rcodedata \
#       /tmp/repro/statmed2014
#   git -C /tmp/repro/statmed2014 checkout 3ae639c53023b74fd1c332aa0bae3ea83b9dc594
#   Rscript tools/reproduce_statmed2014.R [/tmp/repro/statmed2014] [out.json]
#
# Needs survival, dlnm, splines, foreign, xtable and jsonlite.
################################################################################

PINNED_SHA <- "3ae639c53023b74fd1c332aa0bae3ea83b9dc594"

args <- commandArgs(trailingOnly = TRUE)
src <- normalizePath(if (length(args) >= 1) args[1] else "/tmp/repro/statmed2014")
out <- if (length(args) >= 2) args[2] else "tests/reference/statmed2014.json"
out <- file.path(normalizePath(dirname(out), mustWork = FALSE), basename(out))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

sha <- system2("git", c("-C", shQuote(src), "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(sha, PINNED_SHA))
  stop("authors' repository is at ", sha, ", expected ", PINNED_SHA)
csvsha <- sub(" .*", "", system2("sha256sum", shQuote(file.path(src, "uminers.csv")),
  stdout = TRUE))

suppressPackageStartupMessages(library(jsonlite))

# the scripts write their figures with pdf("name.pdf"); send them to a
# temporary directory so the clone stays as checked out
figdir <- file.path(tempdir(), "statmed2014_figures")
dir.create(figdir, showWarnings = FALSE)
pdf <- function(file, ...) grDevices::pdf(file.path(figdir, file), ...)

owd <- setwd(src)
t0 <- proc.time()[["elapsed"]]
for (f in c("01.prep.R", "02.descriptive.R", "03.dlm.R", "04.dlnm.R",
  "05.prediction.R", "06.additional.R")) {
  cat("sourcing", f, "\n")
  suppressMessages(source(f, local = globalenv()))
}
graphics.off()
setwd(owd)
runtime <- proc.time()[["elapsed"]] - t0

################################################################################
# helpers

um <- function(x) { x <- unname(x); if (is.matrix(x)) dimnames(x) <- NULL; x }

# coxph stops when the log-likelihood changes by less than eps = 1e-9, relative.
# For model 8 that leaves the coefficients 1.6e-6 standard errors and the vcov
# 3.3e-6 short of a refit at eps = 1e-14 (models 5 and 7: ~3e-9; the rest
# ~1e-11). Each fit is therefore stored twice: as the authors' scripts leave it
# ("published") and refitted to convergence ("converged"), so the Python side
# can be measured against the maximum rather than against R's stopping rule.
# The warning that toler.chol is not below eps is harmless here: no column of
# these designs is close to aliased.
tight <- function(m) suppressWarnings(
  update(m, control = coxph.control(eps = 1e-14, iter.max = 100)))

fmod1 <- function(m) list(
  coef = um(coef(m)), names = names(coef(m)), vcov = um(vcov(m)),
  se = um(sqrt(diag(vcov(m)))),
  loglik = m$loglik[2], loglik_null = m$loglik[1],
  AIC = AIC(m), BIC = BIC(m), n = m$n, nevent = m$nevent, ties = m$method,
  iter = m$iter)
fmod <- function(m) list(published = fmod1(m), converged = fmod1(tight(m)))

fpred <- function(p) list(
  predvar = um(p$predvar), lag = um(p$lag), bylag = p$bylag,
  cen = if (is.null(p$cen)) NULL else p$cen,
  coef = um(p$coefficients), vcov = um(p$vcov),
  matfit = um(p$matfit), matse = um(p$matse),
  allfit = um(p$allfit), allse = um(p$allse),
  allRRfit = um(p$allRRfit), allRRlow = um(p$allRRlow), allRRhigh = um(p$allRRhigh),
  link = p$model.link)

fcb <- function(cb) list(dim = dim(cb), df = um(attr(cb, "df")),
  range = um(attr(cb, "range")), colsums = um(colSums(cb)),
  rows = um(cb[ROWS, , drop = FALSE]))

fred <- function(r) list(coef = um(coef(r)), vcov = um(vcov(r)),
  predvar = um(r$predvar), fit = um(r$fit), se = um(r$se))

flr <- function(m1, m2) {
  ll <- abs(m1$loglik - m2$loglik)[2]
  df <- abs(length(coef(m1)) - length(coef(m2)))
  list(stat = 2 * ll, df = df, pvalue = 1 - pchisq(2 * ll, df))
}

# a few fixed rows of the design matrices, spread over the data
ROWS <- unique(round(seq(1, nrow(mainsub), length.out = 25)))

################################################################################
# data preparation (01.prep.R)

prep <- list(
  n_uminers = nrow(uminers), n_main = nrow(main),
  n_cases = sum(main$ind), n_ftime = length(ftime), ftime = ftime,
  # mainspl as 01.prep.R leaves it: only the rows ending at a failure time
  n_mainspl = nrow(mainspl), n_mainsub = nrow(mainsub),
  # the sampled risk sets, row by row: this pins down R's set.seed()/sample()
  mainsub_record = um(mainsub$record), mainsub_riskset = um(mainsub$riskset),
  n_mainsub_cases = sum(mainsub$ind), mainsub_agest_sum = sum(mainsub$agest),
  rows = ROWS,
  rexp_rowsums = um(rowSums(rexp)), sexp_rowsums = um(rowSums(sexp)),
  rexp_colsums = um(colSums(rexp)), sexp_colsums = um(colSums(sexp)),
  Qx_colsums = um(colSums(Qx)), Qz_colsums = um(colSums(Qz)),
  Qx_rows = um(Qx[ROWS, ]), Qz_rows = um(Qz[ROWS, ]),
  caltime_sum = sum(caltime), caltime_rows = um(caltime[ROWS]),
  cbz = fcb(cbz))

################################################################################
# Table I (02.descriptive.R)

table1_out <- um(table1)

################################################################################
# Table II and the models behind Figures 1-3 (03.dlm.R, 04.dlnm.R)

models <- list(); preds <- list(); cbs <- list()
for (i in 1:8) {
  models[[paste0("model", i)]] <- fmod(get(paste0("model", i)))
  preds[[paste0("pred", i)]] <- fpred(get(paste0("pred", i)))
  cbs[[paste0("cbx", i)]] <- fcb(get(paste0("cbx", i)))
}
table2 <- list(rows = um(table2), df = as.numeric(table2[, "df"]),
  AIC = sapply(1:8, function(i) AIC(get(paste0("model", i)))),
  BIC = sapply(1:8, function(i) BIC(get(paste0("model", i)))))

# HR for 100 WLM/year at lags 10-12 (text of section 5)
lag1012 <- lapply(list(pred1 = pred1, pred4 = pred4, pred8 = pred8), function(p)
  list(fit = um(p$matRRfit["100", paste0("lag", 10:12)]),
    low = um(p$matRRlow["100", paste0("lag", 10:12)]),
    high = um(p$matRRhigh["100", paste0("lag", 10:12)])))

lrtests <- list(m3_m1 = flr(model3, model1), m3_m4 = flr(model3, model4),
  m8_m4 = flr(model8, model4), m8_m5 = flr(model8, model5),
  m8_m7 = flr(model8, model7))

# peak lag at each exposure level, model 8
peak8 <- as.numeric(sub("lag", "",
  with(pred8, apply(matRRfit, 1, function(x) colnames(matRRfit)[which.max(x)]))))

################################################################################
# Table III and Figure 5 (05.prediction.R)

tab3pred <- list()
for (i in c(1, 4, 5, 8)) {
  # as in the authors' loop: crosspred matches coefficients by the basis name
  cbx <- get(paste0("cbx", i)); model <- get(paste0("model", i))
  tab3pred[[paste0("model", i)]] <- fpred(crosspred(cbx, model, at = hist, cen = 0))
}
predictions <- list(hist = um(hist), table3 = um(table3), tab3 = tab3pred,
  dynhist = um(dynhist),
  predhist4 = fpred(predhist4), predhist5 = fpred(predhist5),
  predhist8 = fpred(predhist8))

################################################################################
# Additional analyses (06.additional.R)

additional <- list(
  rc = list(cb = fcb(cbxrc), model = fmod(modelrc), pred = fpred(predrc)),
  sub = list(n = sum(sub), subjects = length(unique(mainsub[sub, ]$record)),
    model4 = fmod(model4sub), pred4 = fpred(pred4sub),
    model8 = fmod(model8sub), pred8 = fpred(pred8sub)),
  log = list(cb = fcb(cbxlog), model = fmod(modellog), pred = fpred(predlog)),
  knots = list(
    model8b = fmod(model8b), pred8b = fpred(pred8b),
    model8c = fmod(model8c), pred8c = fpred(pred8c),
    model8d = fmod(model8d), pred8d = fpred(pred8d)))

################################################################################
# Reductions of model 8, not in the paper: the overall cumulative
# exposure-response, the lag-response at 100 WLM/year (Figure 3, left) and the
# exposure-response at lag 15 (Figure 3, right) as one-dimensional bases. They
# are taken from the converged refit: the overall reduction sums 39 lags, and
# the cancellation magnifies the 3.3e-6 convergence error of the published fit
# to 2e-4 in units of the reduced standard errors.
model8c <- tight(model8)
reductions <- list(
  overall = fred(crossreduce(cbx8, model8c, cen = 0)),
  var100 = fred(crossreduce(cbx8, model8c, type = "var", value = 100, cen = 0)),
  lag15 = fred(crossreduce(cbx8, model8c, type = "lag", value = 15, cen = 0)))

################################################################################

versions <- list(
  R = R.version.string, platform = R.version$platform,
  dlnm = as.character(packageVersion("dlnm")),
  survival = as.character(packageVersion("survival")),
  splines = as.character(packageVersion("splines")),
  jsonlite = as.character(packageVersion("jsonlite")),
  RNGkind = RNGkind(),
  authors_repo = "https://github.com/gasparrini/2014_gasparrini_StatMed_Rcodedata",
  authors_commit = sha, uminers_sha256 = csvsha,
  date = format(Sys.time(), "%Y-%m-%d"),
  runtime_seconds = runtime)

res <- list(versions = versions, prep = prep, table1 = table1_out,
  table2 = table2, models = models, preds = preds, cbs = cbs,
  lag1012 = lag1012, lrtests = lrtests, peak8 = peak8,
  predictions = predictions, additional = additional, reductions = reductions)

writeLines(toJSON(res, digits = I(17), matrix = "rowmajor", auto_unbox = TRUE,
  na = "null", null = "null", pretty = FALSE), out)
cat("wrote", out, sprintf("(%.1f s in the authors' scripts)\n", runtime))
