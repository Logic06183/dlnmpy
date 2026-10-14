################################################################################
# Reference values for examples/leone2014.py
#
#   Gasparrini A, Leone M. Attributable risk from distributed lag models.
#   BMC Medical Research Methodology 2014;14:55. doi:10.1186/1471-2288-14-55
#
# Runs the authors' own code (01.runmodel.R and attrdl.R, unmodified) from a
# clone of https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata
# pinned to commit 68a8d444f85c673868bdd1511619dac93af6f563, then evaluates
# every attributable quantity in Tables 1 and 2 of the paper (London) plus the
# calls in 02.results.R and 03.graphs.R, and writes them to
# tests/reference/leone2014.json together with the standard normal draws
# attrdl() used for the empirical intervals, so that Python can replay them.
#
# Usage, from the repository root:
#
#   git clone https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata /tmp/repro/leone2014
#   git -C /tmp/repro/leone2014 checkout 68a8d444f85c673868bdd1511619dac93af6f563
#   Rscript tools/reproduce_leone2014.R [/tmp/repro/leone2014]
#
# Needs dlnm, tsModel, jsonlite (and splines, foreign, loaded by 01.runmodel.R).
################################################################################

SHA <- "68a8d444f85c673868bdd1511619dac93af6f563"
NSIM <- 5000     # as in the paper ("5,000 samples"); 02.results.R uses 1000
SEED <- 20140421 # attrdl() draws with rnorm() and no seed; we fix one

args <- commandArgs(trailingOnly = TRUE)
src <- normalizePath(if (length(args)) args[1] else "/tmp/repro/leone2014")
out <- file.path(getwd(), "tests", "reference", "leone2014.json")
if (!file.exists(file.path(getwd(), "src", "dlnmpy")))
  stop("run from the dlnmpy repository root")

head <- system2("git", c("-C", shQuote(src), "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(head, SHA)) stop("authors' repository is at ", head, ", expected ", SHA)

suppressPackageStartupMessages(library(jsonlite))

# the authors' scripts read london.csv and source attrdl.R from their own folder
owd <- setwd(src)
invisible(capture.output(source("01.runmodel.R", echo = FALSE)))
source("attrdl.R")
setwd(owd)

################################################################################
# MODEL

ind <- grep("^cb", names(coef(model)))
cbcoef <- unname(coef(model)[ind])
cbvcov <- unname(vcov(model)[ind, ind])
k <- length(cbcoef)

# the 1st and 99th percentiles are the cut-offs for mild and extreme
perc <- unname(perc)

################################################################################
# ATTRIBUTABLE RISK: every cell of Tables 1 and 2 (London), backward and
# forward, number and fraction, with empirical CIs from NSIM draws

ranges <- list(
  overall      = NULL,
  cold         = c(-100, cen),
  heat         = c(cen, 100),
  extreme_cold = c(-100, perc[1]),
  mild_cold    = c(perc[1], cen),
  mild_heat    = c(cen, perc[2]),
  extreme_heat = c(perc[2], 100))

# attrdl() draws X <- matrix(rnorm(k * nsim), nsim) as its first use of the
# RNG, so after set.seed(SEED) these are exactly the deviates it uses. Every
# simulated call below starts from the same seed (common random numbers).
set.seed(SEED)
Z <- matrix(rnorm(k * NSIM), NSIM)

attr <- list()
for (comp in names(ranges)) {
  for (dir in c("back", "forw")) {
    rg <- ranges[[comp]]
    an <- attrdl(lndn$tmean, cb, lndn$death, model, type = "an", dir = dir,
      cen = cen, range = rg)
    af <- attrdl(lndn$tmean, cb, lndn$death, model, type = "af", dir = dir,
      cen = cen, range = rg)
    set.seed(SEED)
    afsim <- attrdl(lndn$tmean, cb, lndn$death, model, type = "af", dir = dir,
      cen = cen, range = rg, sim = TRUE, nsim = NSIM)
    set.seed(SEED)
    ansim <- attrdl(lndn$tmean, cb, lndn$death, model, type = "an", dir = dir,
      cen = cen, range = rg, sim = TRUE, nsim = NSIM)
    attr[[paste(comp, dir, sep = "_")]] <- list(
      component = comp, dir = dir, range = rg,
      an = an, af = af,
      an_ci = unname(quantile(ansim, c(0.025, 0.975))),
      af_ci = unname(quantile(afsim, c(0.025, 0.975))))
  }
}

# check that the exported deviates are the ones attrdl() used: rebuild one
# set of simulated fractions from Z with the code of attrdl() and compare
chk <- local({
  e <- eigen(cbvcov)
  coefsim <- cbcoef + e$vectors %*% diag(sqrt(e$values), k) %*% t(Z)
  set.seed(SEED)
  r <- attrdl(lndn$tmean, cb, lndn$death, model, cen = cen, sim = TRUE, nsim = NSIM)
  at <- tsModel::Lag(lndn$tmean, 0:25)
  X <- dlnm:::mkXpred("cb", cb, at, seq(nrow(at)), dlnm:::seqlag(c(0, 25)), cen)
  Xall <- Reduce(`+`, lapply(0:25, function(i) X[seq(nrow(at)) + nrow(at) * i, ]))
  m <- (1 - exp(-Xall %*% coefsim)) * lndn$death
  ok <- !is.na(m[, 1])
  max(abs(colSums(m[ok, ]) / sum(lndn$death[ok]) - r))
})
if (chk > 1e-12) stop("exported deviates do not reproduce attrdl() draws: ", chk)

################################################################################
# DAILY SERIES (tot = FALSE): cold as in 02.results.R, heat as in 03.graphs.R

daily <- list(
  cold_back = attrdl(lndn$tmean, cb, lndn$death, model, tot = FALSE, type = "an",
    cen = cen, range = c(-100, cen)),
  cold_forw = attrdl(lndn$tmean, cb, lndn$death, model, tot = FALSE, type = "an",
    dir = "forw", cen = cen, range = c(-100, cen)),
  heat_back = attrdl(lndn$tmean, cb, lndn$death, model, tot = FALSE, type = "an",
    cen = cen, range = c(cen, 100)),
  heat_forw = attrdl(lndn$tmean, cb, lndn$death, model, tot = FALSE, type = "an",
    dir = "forw", cen = cen, range = c(cen, 100)))
# the leading (backward) or trailing (forward) 25 days are NA, written as null

################################################################################
# THE SAME MODEL RUN TO CONVERGENCE (not in the authors' code)
#
# glm() stops after 3 IRLS iterations at its default epsilon = 1e-8, and its
# reported dispersion then uses the working weights of the previous iteration,
# 6.5e-6 (relative) away from the Pearson statistic of the final fit. The
# same model refitted with epsilon = 1e-15 separates porting differences from
# that convergence tolerance. Only point estimates are recomputed.

model_conv <- update(model, control = glm.control(epsilon = 1e-15, maxit = 100))
indc <- grep("^cb", names(coef(model_conv)))
conv_attr <- list()
for (comp in names(ranges)) {
  for (dir in c("back", "forw")) {
    rg <- ranges[[comp]]
    conv_attr[[paste(comp, dir, sep = "_")]] <- list(component = comp, dir = dir,
      an = attrdl(lndn$tmean, cb, lndn$death, model_conv, type = "an", dir = dir,
        cen = cen, range = rg),
      af = attrdl(lndn$tmean, cb, lndn$death, model_conv, type = "af", dir = dir,
        cen = cen, range = rg))
  }
}
cp2conv <- crosspred(cb, model_conv, at = perc, bylag = 0.2, cen = cen)

################################################################################
# PREDICTIONS (01.runmodel.R): overall curve, RR at the percentiles, and the
# lag-response at 25, 27, 29 C (Figure 3)

pred <- list(
  cp1 = list(at = unname(cp1$predvar), allRRfit = unname(cp1$allRRfit),
    allRRlow = unname(cp1$allRRlow), allRRhigh = unname(cp1$allRRhigh)),
  cp2 = list(at = unname(cp2$predvar), allRRfit = unname(cp2$allRRfit),
    allRRlow = unname(cp2$allRRlow), allRRhigh = unname(cp2$allRRhigh)),
  cp3 = list(at = unname(cp3$predvar), bylag = 0.2,
    matRRfit = unname(cp3$matRRfit), matRRlow = unname(cp3$matRRlow),
    matRRhigh = unname(cp3$matRRhigh)))

################################################################################

res <- list(
  versions = list(
    R = R.version.string,
    dlnm = as.character(packageVersion("dlnm")),
    mvmeta = tryCatch(as.character(packageVersion("mvmeta")), error = function(e) NULL),
    mixmeta = tryCatch(as.character(packageVersion("mixmeta")), error = function(e) NULL),
    tsModel = as.character(packageVersion("tsModel")),
    authors_repo = "https://github.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata",
    authors_commit = SHA,
    london_csv_md5 = unname(tools::md5sum(file.path(src, "london.csv"))),
    attrdl_R_md5 = unname(tools::md5sum(file.path(src, "attrdl.R"))),
    note = "mvmeta/mixmeta are not used by this single-city analysis; recorded for completeness"),
  settings = list(maxlag = maxlag, cen = cen, perc = perc, varknots = unname(vk),
    lagknots = unname(lk), nsim = NSIM, seed = SEED,
    cen_percentile = sum(lndn$tmean < cen, na.rm = TRUE) / sum(!is.na(lndn$tmean)),
    ndeaths = sum(lndn$death)),
  model = list(coef = cbcoef, vcov = cbvcov, dispersion = summary(model)$dispersion,
    deviance = deviance(model), df_residual = df.residual(model),
    iter = model$iter, converged = model$converged,
    pearson_dispersion = sum((model$y - fitted(model))^2 / fitted(model)) / df.residual(model),
    # attrdl()'s square root of vcov, for replaying the draws independently
    # of the sign convention of a particular eigen solver
    eigen_vectors = eigen(cbvcov)$vectors, eigen_values = eigen(cbvcov)$values),
  model_converged = list(epsilon = 1e-15, iter = model_conv$iter,
    coef = unname(coef(model_conv)[indc]), vcov = unname(vcov(model_conv)[indc, indc]),
    dispersion = summary(model_conv)$dispersion,
    attributable = unname(conv_attr),
    cp2 = list(allRRfit = unname(cp2conv$allRRfit), allRRlow = unname(cp2conv$allRRlow),
      allRRhigh = unname(cp2conv$allRRhigh))),
  attributable = unname(attr),
  daily = daily,
  predictions = pred,
  normals = Z)

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(res, digits = NA, auto_unbox = TRUE, null = "null", na = "null", pretty = FALSE), out)
cat("wrote", out, "\n")
cat(sprintf("%-14s %-5s  AF%% %6.2f (%5.2f-%5.2f)  AN %8.0f\n",
  sapply(attr, `[[`, "component"), sapply(attr, `[[`, "dir"),
  100 * sapply(attr, `[[`, "af"), 100 * sapply(attr, function(a) a$af_ci[1]),
  100 * sapply(attr, function(a) a$af_ci[2]), sapply(attr, `[[`, "an")), sep = "")
