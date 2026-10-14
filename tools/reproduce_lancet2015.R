################################################################################
# Reference values for examples/lancet2015_vs_r.py
#
#   Gasparrini A, Guo Y, Hashizume M, et al. Mortality risk attributable to
#   high and low ambient temperature: a multicountry observational study.
#   The Lancet 2015;386(9991):369-375. doi:10.1016/S0140-6736(14)62114-0
#
# Runs the authors' own code (00.prepdata.R to 04.tables.R and attrdl.R,
# unmodified) from a clone of
# https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata pinned to
# commit 7c3e03c425529d2f210b1659a01bbc582448b9c0, on its England & Wales
# data, and writes every stage intermediate to tests/reference/lancet2015.json:
#
#   1. per-region reduced coefficients and covariance (overall cumulative)
#   2. the multivariate meta-regression on mean and range of temperature:
#      coefficients, covariance, Psi, log-likelihood; BLUPs and their vcov
#   3. minimum-mortality percentile and temperature per region
#   4. attributable deaths (forward) per region and in total, with the
#      empirical intervals, plus the standard normal deviates attrdl() drew,
#      so that Python can replay them exactly
#
# Two chains are recorded.
#
#   "published": the authors' code as written: glm() at its default
#     epsilon = 1e-8, mvmeta() for the second stage.
#   "converged" (not in the authors' code): the same first-stage models
#     refitted with glm.control(epsilon = 1e-14), and the second stage fitted
#     with mixmeta(), the successor of mvmeta that dlnmpy's companion
#     mixmetapy ports. Stages 3 and 4 are then rerun from these BLUPs, stage 4
#     by sourcing the authors' 03.attr.R again with the same seed, so the
#     same deviates are drawn.
#
# The cross-fits (mvmeta on the converged inputs, mixmeta on the published
# ones) are recorded too, to show how far two R implementations of the same
# REML model sit apart on this flat likelihood.
#
# Usage, from the repository root:
#
#   git clone https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata /tmp/repro/lancet2015
#   git -C /tmp/repro/lancet2015 checkout 7c3e03c425529d2f210b1659a01bbc582448b9c0
#   Rscript tools/reproduce_lancet2015.R [/tmp/repro/lancet2015]
#
# Needs dlnm, mvmeta, mixmeta, tsModel, jsonlite. About 40 seconds.
################################################################################

SHA <- "7c3e03c425529d2f210b1659a01bbc582448b9c0"
SEED <- 20150725 # attrdl() draws with rnorm() and the authors set no seed

args <- commandArgs(trailingOnly = TRUE)
src <- normalizePath(if (length(args)) args[1] else "/tmp/repro/lancet2015")
out <- file.path(getwd(), "tests", "reference", "lancet2015.json")
if (!file.exists(file.path(getwd(), "src", "dlnmpy")))
  stop("run from the dlnmpy repository root")

head <- system2("git", c("-C", shQuote(src), "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(head, SHA)) stop("authors' repository is at ", head, ", expected ", SHA)
dirty <- system2("git", c("-C", shQuote(src), "status", "--porcelain", "--untracked-files=no"),
  stdout = TRUE)
if (length(dirty)) stop("authors' scripts are modified:\n", paste(dirty, collapse = "\n"))

suppressPackageStartupMessages({library(jsonlite); library(mixmeta)})

# attrdl() draws its deviates as X <- matrix(rnorm(k * nsim), nsim). attrdl.R
# is sourced into the global environment, so a global rnorm() is found first:
# this one records each call's deviates without changing the random stream.
.Zrec <- list()
rnorm <- function(n, mean = 0, sd = 1) {
  z <- stats::rnorm(n, mean, sd)
  .Zrec[[length(.Zrec) + 1L]] <<- z
  z
}

################################################################################
# THE AUTHORS' CODE, AS PUBLISHED

owd <- setwd(src)
invisible(capture.output({
  source("00.prepdata.R")
  source("01.firststage.R")
  source("02.secondstage.R")
  set.seed(SEED)
  source("03.attr.R")
  source("04.tables.R")
}))
setwd(owd)
if (length(.Zrec) != 3L * length(dlist)) stop("expected ", 3 * length(dlist), " rnorm calls")
Zpub <- .Zrec

# keep the published chain before anything below overwrites the globals
pub <- list(coef = unname(coef), vcov = unname(vcov), mv = mv, blup = blup,
  minperc = unname(minperccity), mintemp = unname(mintempcity),
  an = unname(matsim), anlow = unname(ancitylow), anhigh = unname(ancityhigh),
  antot = unname(antot), antotlow = unname(antotlow), antothigh = unname(antothigh),
  af = unname(afcity), aflow = unname(afcitylow), afhigh = unname(afcityhigh),
  aftot = unname(aftot), aftotlow = unname(aftotlow), aftothigh = unname(aftothigh),
  arraysim = arraysim, minperccountry = minperccountry, totdeath = unname(totdeath))

# check that the recorded deviates are the ones attrdl() used: rebuild a few
# draws per call from them and push each through attrdl() as a point estimate
checkdraws <- function(Zlist, blup, mintemp, arraysim, nchk = 3) {
  worst <- 0
  for (i in seq(dlist)) {
    data <- dlist[[i]]
    argvar <- list(x = data$tmean, fun = varfun, knots = quantile(data$tmean,
      varper / 100, na.rm = T), degree = vardegree)
    cb <- crossbasis(data$tmean, lag = lag, argvar = argvar,
      arglag = list(knots = logknots(lag, lagnk)))
    e <- eigen(blup[[i]]$vcov)
    for (j in 1:3) {
      Z <- matrix(Zlist[[3 * (i - 1) + j]], nsim)
      cs <- blup[[i]]$blup + e$vectors %*% diag(sqrt(e$values), 5) %*% t(Z)
      rg <- list(NULL, c(-100, mintemp[i]), c(mintemp[i], 100))[[j]]
      for (s in seq(nchk)) {
        a <- attrdl(data$tmean, cb, data$death, coef = cs[, s], vcov = blup[[i]]$vcov,
          type = "an", dir = "forw", cen = mintemp[i], range = rg)
        worst <- max(worst, abs(a - arraysim[i, j, s]) / abs(arraysim[i, j, s]))
      }
    }
  }
  worst
}
chk <- checkdraws(Zpub, pub$blup, pub$mintemp, pub$arraysim)
if (chk > 1e-12) stop("recorded deviates do not reproduce attrdl() draws: ", chk)

################################################################################
# FIRST STAGE RUN TO CONVERGENCE (not in the authors' code)
#
# The authors' loop with one change: glm.control(epsilon = 1e-14). glm()'s
# default epsilon = 1e-8 stops IRLS about 1e-7 (relative) short of the
# optimum in the reduced coefficients; statsmodels iterates further.

coefc <- matrix(NA, nrow(cities), length(varper) + vardegree, dimnames = list(cities$city))
vcovc <- vector("list", nrow(cities)); names(vcovc) <- cities$city
iterpub <- iterc <- dispc <- disppub <- rep(NA, nrow(cities))
for (i in seq(length(dlist))) {
  data <- dlist[[i]]
  argvar <- list(fun = varfun, knots = quantile(data$tmean, varper / 100, na.rm = T),
    degree = vardegree)
  cb <- crossbasis(data$tmean, lag = lag, argvar = argvar,
    arglag = list(knots = logknots(lag, lagnk)))
  model0 <- glm(formula, data, family = quasipoisson, na.action = "na.exclude")
  model <- glm(formula, data, family = quasipoisson, na.action = "na.exclude",
    control = glm.control(epsilon = 1e-14, maxit = 100))
  cen <- mean(data$tmean, na.rm = T)
  red0 <- crossreduce(cb, model0, cen = cen)
  if (max(abs(coef(red0) - pub$coef[i, ])) > 0) stop("published first stage not reproduced")
  red <- crossreduce(cb, model, cen = cen)
  coefc[i, ] <- coef(red)
  vcovc[[i]] <- vcov(red)
  iterpub[i] <- model0$iter; iterc[i] <- model$iter
  disppub[i] <- summary(model0)$dispersion; dispc[i] <- summary(model)$dispersion
}

################################################################################
# SECOND STAGE: mvmeta (authors) and mixmeta, on both sets of inputs

avgtmean <- sapply(dlist, function(x) mean(x$tmean, na.rm = T))
rangetmean <- sapply(dlist, function(x) diff(range(x$tmean, na.rm = T)))

# flat vector and vcov in predictor-major order (all intercepts, then all
# avgtmean slopes, then all rangetmean slopes), which is mixmeta's and dlnmpy's;
# mvmeta's coef() and vcov() are outcome-major
metafit <- function(fit) {
  p <- fit$dim$p; k <- fit$dim$k
  v <- vcov(fit)
  b <- coef(fit)
  if (inherits(fit, "mvmeta")) {
    B <- fit$coefficients                                   # (p, k)
    perm <- as.vector(t(matrix(seq(p * k), p, k)))
    v <- v[perm, perm]
    b <- b[perm]
  } else {
    B <- matrix(b, p, k, byrow = TRUE)
  }
  if (!identical(unname(as.vector(t(B))), unname(b))) stop("unexpected coefficient order")
  bl <- blup(fit, vcov = TRUE)
  list(package = class(fit)[1], coef = unname(B), vcov = unname(v), Psi = unname(fit$Psi),
    logLik = as.numeric(logLik(fit)), converged = fit$converged, niter = fit$niter,
    blup = unname(t(sapply(bl, `[[`, "blup"))),
    blupvcov = unname(lapply(bl, function(b) unname(b$vcov))))
}

mvpub <- pub$mv                                                     # authors' fit
mxpub <- mixmeta(pub$coef ~ avgtmean + rangetmean, pub$vcov, method = "reml")
mvconv <- mvmeta(coefc ~ avgtmean + rangetmean, vcovc, data = cities)
mxconv <- mixmeta(coefc ~ avgtmean + rangetmean, vcovc, method = "reml")
# and the same with the optimiser's relative tolerance at 1e-14 in place of
# its default sqrt(.Machine$double.eps) = 1.5e-8: on this flat surface the
# default stops ~1e-8 short in the log-likelihood and ~1e-5 in the parameters
mxtight <- mixmeta(coefc ~ avgtmean + rangetmean, vcovc, method = "reml",
  control = list(reltol = 1e-14, maxiter = 1000))

# same fit by mvmeta as the authors' code produced: the BLUPs it used
if (max(abs(t(sapply(pub$blup, `[[`, "blup")) - unname(blup(mvpub)))) > 0)
  stop("mvmeta BLUPs not reproduced")

################################################################################
# STAGES 3 AND 4 FROM THE CONVERGED CHAIN (mixmeta BLUPs)
#
# Stage 3 is the authors' loop from 02.secondstage.R; stage 4 is their
# 03.attr.R sourced again, unchanged, with the globals it reads (blup,
# mintempcity) set to the converged chain and the same seed.

blup <- blup(mxconv, vcov = TRUE)
names(blup) <- cities$city
minperccity <- mintempcity <- rep(NA, length(dlist))
names(mintempcity) <- names(minperccity) <- cities$city
for (i in seq(length(dlist))) {
  data <- dlist[[i]]
  predvar <- quantile(data$tmean, 1:99 / 100, na.rm = T)
  argvar <- list(x = predvar, fun = varfun,
    knots = quantile(data$tmean, varper / 100, na.rm = T), degree = vardegree,
    Bound = range(data$tmean, na.rm = T))
  bvar <- do.call(onebasis, argvar)
  minperccity[i] <- (1:99)[which.min((bvar %*% blup[[i]]$blup))]
  mintempcity[i] <- quantile(data$tmean, minperccity[i] / 100, na.rm = T)
}

.Zrec <- list()
owd <- setwd(src)
invisible(capture.output({set.seed(SEED); source("03.attr.R")}))
setwd(owd)
Zconv <- .Zrec
if (!identical(Zconv, Zpub)) stop("the converged chain did not draw the same deviates")
chk2 <- checkdraws(Zconv, blup, mintempcity, arraysim)
if (chk2 > 1e-12) stop("recorded deviates do not reproduce attrdl() draws (converged): ", chk2)
rm(rnorm)

conv <- list(minperc = unname(minperccity), mintemp = unname(mintempcity),
  an = unname(matsim), anlow = unname(ancitylow), anhigh = unname(ancityhigh),
  antot = unname(antot), antotlow = unname(antotlow), antothigh = unname(antothigh),
  af = unname(afcity), aflow = unname(afcitylow), afhigh = unname(afcityhigh),
  aftot = unname(aftot), aftotlow = unname(aftotlow), aftothigh = unname(aftothigh),
  minperccountry = median(minperccity))

################################################################################

attrblock <- function(a, bl) c(a[c("minperc", "mintemp", "an", "anlow", "anhigh", "antot",
  "antotlow", "antothigh", "af", "aflow", "afhigh", "aftot", "aftotlow", "aftothigh",
  "minperccountry")],
  # the square root of each BLUP vcov that attrdl() used, for the replay
  list(eigen_vectors = unname(lapply(bl, function(b) unname(eigen(b$vcov)$vectors)))))

Z <- lapply(seq(dlist), function(i) lapply(1:3, function(j)
  matrix(Zpub[[3 * (i - 1) + j]], nsim)))

res <- list(
  versions = list(
    R = R.version.string,
    dlnm = as.character(packageVersion("dlnm")),
    mvmeta = as.character(packageVersion("mvmeta")),
    mixmeta = as.character(packageVersion("mixmeta")),
    tsModel = as.character(packageVersion("tsModel")),
    authors_repo = "https://github.com/gasparrini/2015_gasparrini_Lancet_Rcodedata",
    authors_commit = SHA,
    regEngWales_csv_md5 = unname(tools::md5sum(file.path(src, "regEngWales.csv"))),
    attrdl_R_md5 = unname(tools::md5sum(file.path(src, "attrdl.R")))),
  settings = list(regions = cities$city, citynames = cities$cityname,
    varper = varper, vardegree = vardegree, lag = lag, lagnk = lagnk, dfseas = dfseas,
    nsim = nsim, seed = SEED, totdeath = pub$totdeath,
    avgtmean = unname(avgtmean), rangetmean = unname(rangetmean),
    knots = unname(t(sapply(dlist, function(x) quantile(x$tmean, varper / 100, na.rm = T)))),
    cen = unname(avgtmean)),
  first_stage = list(
    published = list(coef = pub$coef, vcov = pub$vcov, iter = iterpub, dispersion = disppub),
    converged = list(epsilon = 1e-14, coef = unname(coefc), vcov = unname(vcovc),
      iter = iterc, dispersion = dispc)),
  second_stage = list(
    published_mvmeta = metafit(mvpub),
    published_mixmeta = metafit(mxpub),
    converged_mvmeta = metafit(mvconv),
    converged_mixmeta = metafit(mxconv),
    converged_mixmeta_reltol14 = metafit(mxtight),
    # REML log-likelihoods of mvmeta and mixmeta differ by a constant
    loglik_constant_mvmeta_minus_mixmeta = as.numeric(logLik(mvconv) - logLik(mxconv))),
  attributable = list(
    published = attrblock(pub, pub$blup),
    converged = attrblock(conv, blup)),
  normals = Z)

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(res, digits = NA, auto_unbox = TRUE, null = "null", na = "null",
  pretty = FALSE), out)
cat("wrote", out, "\n")
cat(sprintf("published: MMP %.1f, AF %% total %.2f (%.2f-%.2f) cold %.2f heat %.2f\n",
  pub$minperccountry, pub$aftot[1], pub$aftotlow[1], pub$aftothigh[1], pub$aftot[2], pub$aftot[3]))
cat(sprintf("converged: MMP %.1f, AF %% total %.2f (%.2f-%.2f) cold %.2f heat %.2f\n",
  conv$minperccountry, conv$aftot[1], conv$aftotlow[1], conv$aftothigh[1], conv$aftot[2],
  conv$aftot[3]))
cat(sprintf("glm iterations published %s, converged %s\n",
  paste(iterpub, collapse = " "), paste(iterc, collapse = " ")))
cat(sprintf("logLik mvmeta %.10f mixmeta %.10f (published inputs); %.10f %.10f (converged); %.10f (reltol 1e-14)\n",
  logLik(mvpub), logLik(mxpub), logLik(mvconv), logLik(mxconv), logLik(mxtight)))
