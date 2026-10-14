# Reference values for examples/biomet2017.py: the penalised models of
#
#   Gasparrini A, Scheipl F, Armstrong B, Kenward MG. A penalized framework for
#   distributed lag non-linear models. Biometrics 2017;73(3):938-948.
#   doi:10.1111/biom.12645
#
# run with the authors' public code and data, pinned to one commit of
# https://github.com/gasparrini/2017_gasparrini_Biomet_Rcodedata
#
# Usage, from the repository root:
#   Rscript tools/reproduce_biomet2017.R [path/to/authors/clone]
# The clone defaults to /tmp/repro/biomet2017 and is created (and checked out at
# the pinned commit) if it does not exist. Writes tests/reference/biomet2017.json.
#
# What is run, and how it departs from the authors' scripts:
#
# * Example 1 (London, ex1_01.ext.R): the three penalised GAMs with the
#   "external" method, gam(..., paraPen=list(cb=cbPen(cb)), method="REML"),
#   copied line for line: gam1 (default penalties), gam2 (doubly varying lag
#   penalty, eq. 8b and 7a), gam3 (double-threshold exposure-response,
#   penalised lag). The two unpenalised GLMs and the Q-AIC knot search are not
#   rerun; they are not part of the penalised framework.
# * Example 2 (Colorado uranium miners): ex2_01.prep.R is sourced unchanged.
#   The penalised GAM of ex2_02.ext.R is run as written. The authors do not set
#   'method', so mgcv uses its default for a known-scale family, UBRE. dlnmpy has
#   no UBRE criterion, so the same model is also fitted with method="REML"
#   (an addition, not in the authors' code) to give dlnmpy's own smoothing
#   parameter selection something to be compared with. The unpenalised GLM is
#   not rerun.
# * The "internal" method of ex1_01.int.R / ex2_02.int.R, s(Q, L, bs="cb"), is
#   skipped: dlnmpy has no smooth-constructor interface. The authors' paper
#   reports that the two methods give the same fits.
# * The simulation study (simul/) is not rerun.
# * Plots are not drawn.

suppressPackageStartupMessages({
  library(dlnm); library(mgcv); library(splines); library(jsonlite)
})

SHA <- "531c5c8bec3d477a32e5fd66e7f37f1e4fdead1d"
URL <- "https://github.com/gasparrini/2017_gasparrini_Biomet_Rcodedata"
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "/tmp/repro/biomet2017"
if (!dir.exists(src)) {
  dir.create(dirname(src), recursive = TRUE, showWarnings = FALSE)
  stopifnot(system2("git", c("clone", "-q", URL, src)) == 0)
}
stopifnot(system2("git", c("-C", src, "checkout", "-q", SHA)) == 0)
head <- system2("git", c("-C", src, "rev-parse", "HEAD"), stdout = TRUE)
if (head != SHA) stop("authors' clone is at ", head, ", expected ", SHA)

sha256 <- function(f) strsplit(system2("sha256sum", f, stdout = TRUE), " ")[[1]][1]
outfile <- "tests/reference/biomet2017.json"
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
t0 <- proc.time()[["elapsed"]]

# ---------------------------------------------------------------------------
# helpers: everything dlnmpy is compared on
penout <- function(pen) {
  mats <- pen[setdiff(names(pen), c("rank", "sp"))]
  list(names = names(mats), mats = lapply(mats, unname), rank = unname(pen$rank))
}
gamout <- function(m, cbname, pen, pred, elapsed) {
  cf <- coef(m)
  ind <- grep(paste0("^", cbname), names(cf))
  list(
    method = m$method, family = m$family$family, converged = m$converged,
    outer_conv = if (!is.null(m$outer.info$conv)) m$outer.info$conv else NA,
    sp = unname(m$sp), score = as.numeric(m$gcv.ubre),
    scale = m$sig2, reml_scale = if (is.null(m$reml.scale)) NA else m$reml.scale,
    deviance = deviance(m),
    coef_names = names(cf), coef = unname(cf), se = unname(sqrt(diag(m$Vp))),
    cb_index = ind, vcov_cb = unname(m$Vp[ind, ind]),
    edf = unname(m$edf), edf_cb = sum(m$edf[ind]),
    penalty = penout(pen),
    pred = list(cen = if (is.null(pred$cen)) NA else pred$cen, predvar = pred$predvar,
                allfit = unname(pred$allfit), allse = unname(pred$allse)),
    seconds = elapsed)
}
timed <- function(expr) {
  t <- system.time(v <- expr)[["elapsed"]]
  list(value = v, seconds = t)
}

out <- list()

# ---------------------------------------------------------------------------
# EXAMPLE 1: London 1993-2006, ex1_01.ext.R (penalised part)
london <- read.csv(file.path(src, "example1", "london.csv"))

# GAM WITH DEFAULT PENALTIES
cbgam1 <- crossbasis(london$tmean, lag = 25, argvar = list(fun = "ps", df = 9),
  arglag = list(fun = "ps", df = 10))
cbgam1Pen <- cbPen(cbgam1)
f <- timed(gam(death ~ cbgam1 + ns(time, 10 * 14) + dow, family = quasipoisson(), london,
  paraPen = list(cbgam1 = cbgam1Pen), method = "REML"))
gam1 <- f$value
predslgam1 <- crosspred(cbgam1, gam1, by = 0.2, bylag = 0.2, cen = 20)
cat("gam1  sp", signif(gam1$sp, 6), " REML", gam1$gcv.ubre, " edf", sum(gam1$edf[2:91]),
    " ", f$seconds, "s\n")
out$ex1_gam1 <- gamout(gam1, "cbgam1", cbgam1Pen, predslgam1, f$seconds)

# GAM WITH DOUBLY VARYING PENALTY ON THE LAG
cbgam2 <- crossbasis(london$tmean, lag = 25, argvar = list(fun = "ps", df = 9),
  arglag = list(fun = "ps", df = 10, fx = T))
C <- do.call("onebasis", c(list(x = 0:25, fun = "ps", df = 10, intercept = T)))
D <- diff(diag(25 + 1), diff = 2)
P <- diag((seq(0, 25 - 2))^2)
Slag1 <- t(C) %*% t(D) %*% P %*% D %*% C
Slag2 <- diag(rep(0:1, c(6, 4)))
cbgam2Pen <- cbPen(cbgam2, addSlag = list(Slag1, Slag2))
f <- timed(gam(death ~ cbgam2 + ns(time, 10 * 14) + dow, family = quasipoisson(), london,
  paraPen = list(cbgam2 = cbgam2Pen), method = "REML"))
gam2 <- f$value
predslgam2 <- crosspred(cbgam2, gam2, by = 0.2, bylag = 0.2, cen = 20)
cat("gam2  sp", signif(gam2$sp, 6), " REML", gam2$gcv.ubre, " edf", sum(gam2$edf[2:91]),
    " ", f$seconds, "s\n")
out$ex1_gam2 <- gamout(gam2, "cbgam2", cbgam2Pen, predslgam2, f$seconds)
out$ex1_gam2$addslag <- list(Slag1 = unname(Slag1), Slag2 = unname(Slag2))

# GAM WITH ONE DIMENSION UNPENALIZED (DOUBLE THRESHOLD)
cbgam3 <- crossbasis(london$tmean, lag = 25, argvar = list(fun = "thr", thr = c(17, 21)),
  arglag = list(fun = "ps", df = 10))
cbgam3Pen <- cbPen(cbgam3, addSlag = list(Slag1, Slag2))
f <- timed(gam(death ~ cbgam3 + ns(time, 10 * 14) + dow, family = quasipoisson(), london,
  paraPen = list(cbgam3 = cbgam3Pen), method = "REML"))
gam3 <- f$value
predslgam3 <- crosspred(cbgam3, gam3, by = 0.2, bylag = 0.2)
cat("gam3  sp", signif(gam3$sp, 6), " REML", gam3$gcv.ubre, " edf", sum(gam3$edf[2:21]),
    " ", f$seconds, "s\n")
out$ex1_gam3 <- gamout(gam3, "cbgam3", cbgam3Pen, predslgam3, f$seconds)

# ---------------------------------------------------------------------------
# EXAMPLE 2: Colorado Plateau uranium miners
# the authors' preparation script, unchanged (it reads uminers.csv from its
# own directory, defines main, Qx, Qz and cbz)
e2 <- new.env()
f <- timed(sys.source(file.path(src, "example2", "ex2_01.prep.R"), envir = e2, chdir = TRUE))
cat("ex2 prep", f$seconds, "s\n")
main <- e2$main; Qx <- e2$Qx; Qz <- e2$Qz; cbz <- e2$cbz
out$ex2_data <- list(nrow = nrow(main), sum_out = sum(main$out), sum_age = sum(main$age),
  sum_year = sum(main$year), colsum_Qx = unname(colSums(Qx)), colsum_Qz = unname(colSums(Qz)),
  sumsq_Qx = sum(Qx^2), sumsq_Qz = sum(Qz^2), prep_seconds = f$seconds)

# GAM (ex2_02.ext.R)
vk <- c(min(Qx), logknots(Qx, nk = 9), max(Qx))
cbgam <- crossbasis(Qx, lag = c(2, 40), argvar = list(fun = "cr", knots = vk),
  arglag = list(fun = "cr", int = F))
C <- do.call("onebasis", c(list(x = 2:40, fun = "cr", df = 10, intercept = F)))
P <- diag(rep(0:1, c(30 - 2 + 1, 40 - 30)))
Slag2 <- t(C) %*% P %*% C
cbgamPen <- cbPen(cbgam, addSlag = Slag2)
out$ex2_spec <- list(vk = vk, Slag2 = unname(Slag2))

# as published: no 'method', so mgcv's default (UBRE for Poisson)
f <- timed(gam(out ~ cbgam + cbz + bs(age, df = 5) + year, data = main, family = poisson,
  paraPen = list(cbgam = cbgamPen)))
gam_ubre <- f$value
predgam <- crosspred(cbgam, gam_ubre, at = 0:25 * 10, cen = 0)
cat("ex2 UBRE  sp", signif(gam_ubre$sp, 6), " UBRE", gam_ubre$gcv.ubre, " edf",
    sum(gam_ubre$edf[2:101]), " ", f$seconds, "s\n")
out$ex2_gam_ubre <- gamout(gam_ubre, "cbgam", cbgamPen, predgam, f$seconds)

# addition: the same model with REML, for dlnmpy's own selection
f <- timed(gam(out ~ cbgam + cbz + bs(age, df = 5) + year, data = main, family = poisson,
  paraPen = list(cbgam = cbgamPen), method = "REML"))
gam_reml <- f$value
predgamr <- crosspred(cbgam, gam_reml, at = 0:25 * 10, cen = 0)
cat("ex2 REML  sp", signif(gam_reml$sp, 6), " REML", gam_reml$gcv.ubre, " edf",
    sum(gam_reml$edf[2:101]), " ", f$seconds, "s\n")
out$ex2_gam_reml <- gamout(gam_reml, "cbgam", cbgamPen, predgamr, f$seconds)

# ---------------------------------------------------------------------------
out$versions <- list(
  R = R.version.string, dlnm = as.character(packageVersion("dlnm")),
  mgcv = as.character(packageVersion("mgcv")), jsonlite = as.character(packageVersion("jsonlite")),
  authors_repo = URL, authors_commit = SHA,
  data_sha256 = list(london.csv = sha256(file.path(src, "example1", "london.csv")),
                     uminers.csv = sha256(file.path(src, "example2", "uminers.csv"))),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  total_seconds = proc.time()[["elapsed"]] - t0)
write_json(out, outfile, digits = NA, auto_unbox = TRUE, na = "null")
cat("written", outfile, "in", round(out$versions$total_seconds, 1), "s\n")
