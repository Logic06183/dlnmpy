################################################################################
# Reference values for examples/vicedo2019.py
#
#   Vicedo-Cabrera AM, Sera F, Gasparrini A. Hands-on tutorial on a modeling
#   framework for projections of climate change impacts on health.
#   Epidemiology 2019;30(3):321-329. doi:10.1097/EDE.0000000000000982
#
# Runs the authors' own scripts, unmodified, from a clone of
# https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata pinned to
# commit 313606cd0fbeec3cd52e72cd4fd5956faece9c48, then writes the quantities
# the Python translation is compared against to tests/reference/vicedo2019.json.
#
# Usage, from the repository root:
#
#   git clone https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata \
#       /tmp/repro/vicedo2019
#   git -C /tmp/repro/vicedo2019 checkout 313606cd0fbeec3cd52e72cd4fd5956faece9c48
#   Rscript tools/reproduce_vicedo2019.R [authors_dir] [out_json]
#
# Needs dlnm, tsModel, MASS and jsonlite. Takes about a minute.
#
# Two things are done around the authors' code, neither of which changes a
# number it computes:
#   - 03DownscalingCalibration.R calls source("fhempel.R") but the file in the
#     repository is fhempel.r, which fails on a case-sensitive file system. A
#     `source` defined in the evaluation environment falls back to the
#     lower-case name.
#   - 04_05_06ExtCurveProjUncert.R ends its loop with rm(ansim). The script is
#     evaluated one expression at a time and ansim is copied before that rm, so
#     the per-GCM projections can be written out too.
# Plots go to a null device.
################################################################################

SHA <- "313606cd0fbeec3cd52e72cd4fd5956faece9c48"
args <- commandArgs(trailingOnly = TRUE)
src <- normalizePath(if (length(args) >= 1) args[1] else "/tmp/repro/vicedo2019")
out <- if (length(args) >= 2) args[2] else "tests/reference/vicedo2019.json"
out <- file.path(normalizePath(dirname(out), mustWork = FALSE), basename(out))

suppressPackageStartupMessages({
  library(jsonlite); library(dlnm); library(MASS); library(tsModel)
})

# ---- the clone must be the pinned commit, with nothing changed ---------------
head <- system2("git", c("-C", shQuote(src), "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(head, SHA)) stop("authors' repository is at ", head, ", expected ", SHA)
dirty <- system2("git", c("-C", shQuote(src), "status", "--porcelain"), stdout = TRUE)
if (length(dirty)) stop("authors' repository has local changes:\n", paste(dirty, collapse = "\n"))

# ---- run the authors' scripts in one environment -----------------------------
env <- new.env(parent = globalenv())
env$source <- function(file, ...) {
  if (!file.exists(file) && file.exists(sub("\\.R$", ".r", file)))
    file <- sub("\\.R$", ".r", file)
  base::source(file, local = env, ...)
}

t0 <- proc.time()[["elapsed"]]
owd <- setwd(src)
pdf(NULL)
for (f in c("01EstimationERassociation.R", "02ProjExpHOutSeries.R",
            "03DownscalingCalibration.R")) {
  cat("running", f, "\n")
  base::source(f, local = env)
}
cat("running 04_05_06ExtCurveProjUncert.R\n")
for (e in parse("04_05_06ExtCurveProjUncert.R")) {
  if (identical(e, quote(rm(ansim)))) assign("ansim_kept", env$ansim, envir = env)
  eval(e, env)
}
cat("\nrunning 07DemChangesAdapt.R\n")
base::source("07DemChangesAdapt.R", local = env)
dev.off()
setwd(owd)
cat(sprintf("authors' scripts ran in %.0f s\n", proc.time()[["elapsed"]] - t0))

# ---- the normal deviates behind each mvrnorm() call ---------------------------
# The script draws coefsim <- mvrnorm(nsim, coef, vcov) after
# set.seed(13041975 + j), the same for both RCPs. mvrnorm turns
# Z <- matrix(rnorm(nsim * k), nsim) into coef + V^(1/2) t(Z) with the eigen
# square root; regenerate Z and check it reproduces the draws exactly.
k <- length(env$coef)
normals <- coefsim_head <- eigvecs <- list()
for (j in seq(env$gcm)) {
  set.seed(13041975 + j)
  Z <- matrix(rnorm(env$nsim * k), env$nsim)
  set.seed(13041975 + j)
  cs <- mvrnorm(env$nsim, env$coef, env$vcov)
  eS <- eigen(env$vcov, symmetric = TRUE)
  rebuilt <- t(drop(env$coef) + eS$vectors %*% diag(sqrt(pmax(eS$values, 0)), k) %*% t(Z))
  stopifnot(max(abs(rebuilt - cs)) == 0)
  normals[[names(env$gcm)[j]]] <- Z
  coefsim_head[[names(env$gcm)[j]]] <- cs[1:10, ]
  eigvecs[[names(env$gcm)[j]]] <- eS$vectors
}

# ---- quantities to compare -----------------------------------------------------
with(env, {
  cbcoef <- function(model, nm = "cb") {
    i <- grep(paste0("^", nm), names(coef(model)))
    list(coef = unname(coef(model)[i]), vcov = unname(vcov(model)[i, i]))
  }
  arr <- function(a) list(dim = dim(a), dimnames = dimnames(a), values = as.vector(a))

  ak <- ansim_kept
  # per-GCM point estimates (absolute AN) and per-GCM 95% empirical intervals
  gcm_est <- ak[, , "abs", , , 1]
  gcm_ci <- apply(ak[, , "abs", , , -1], 1:4, quantile, c(0.025, 0.975))

  # bias-correction parameters behind the calibrated series
  corr <- lapply(list(rcp4p5 = "lndn_rcp4p5.csv", rcp8p5 = "lndn_rcp8p5.csv"), function(f) {
    raw <- read.csv(file.path(src, f), colClasses = c(date = "Date"))
    lapply(fhempel(obs[c("date", "tmean")], raw, output = "correction"),
      function(x) list(add = unname(x[, "add"]), mult = unname(x[, "mult"])))
  })
  rows <- seq(1, nrow(rcp4p5cal), by = 100)
  calib <- lapply(list(rcp4p5 = rcp4p5cal, rcp8p5 = rcp8p5cal), function(d)
    list(rows = rows, colsum = unname(colSums(d[, -1])), colsumsq = unname(colSums(d[, -1]^2)),
      sample = lapply(as.list(d[rows, -1]), unname)))

  ref <<- list(
    versions = list(
      R = R.version.string,
      dlnm = as.character(packageVersion("dlnm")),
      MASS = as.character(packageVersion("MASS")),
      tsModel = as.character(packageVersion("tsModel")),
      splines = as.character(packageVersion("splines")),
      jsonlite = as.character(packageVersion("jsonlite")),
      authors_repo = "https://github.com/gasparrini/2019_vicedo-cabrera_Epidem_Rcodedata",
      authors_commit = SHA,
      generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    nsim = nsim,
    # 01 exposure-response
    stage1 = list(
      ns = cbcoef(m), lin = cbcoef(m1), thr = cbcoef(m2),
      seas_df = round(8 * length(obs$date) / 365.25),
      dispersion = summary(m)$dispersion,
      deviance = deviance(m), nobs = nobs(m)),
    cp_cen19 = list(predvar = cp$predvar, allRRfit = unname(cp$allRRfit)),
    mmt = cen,
    pred = list(predvar = pred$predvar, allRRfit = unname(pred$allRRfit),
      allRRlow = unname(pred$allRRlow), allRRhigh = unname(pred$allRRhigh)),
    pred_lin = list(predvar = pred1$predvar, allRRfit = unname(pred1$allRRfit)),
    pred_thr = list(predvar = pred2$predvar, allRRfit = unname(pred2$allRRfit)),
    # 02 projected mortality
    deathdoy = unname(as.vector(deathdoy)),
    deathperiod = deathperiod,
    # 03 calibration
    calibration = list(correction = corr, series = calib),
    # 04-06 reduced coefficients, draws, projections
    red = list(coef = unname(coef), vcov = unname(vcov)),
    normals = lapply(normals, function(Z) as.vector(Z)),   # column-major (nsim, k)
    coefsim_head = lapply(coefsim_head, function(x) as.vector(x)),
    eigvecs = lapply(eigvecs, function(x) as.vector(x)),   # column-major (k, k)
    gcm_est = arr(gcm_est),
    gcm_ci = arr(gcm_ci),
    anabs = arr(anabs), anrel = arr(anrel), afabs = arr(afabs), afrel = arr(afrel),
    # Figure 5
    fig5 = list(predvar = pred_exp$predvar, allRRfit = unname(pred_exp$allRRfit),
      allRRlow = unname(pred_exp$allRRlow), allRRhigh = unname(pred_exp$allRRhigh),
      anhist = sum(anhist), anproj = sum(anproj)),
    # 07 age-specific curves and adaptation
    age = list(
      mmt = c(cen0_64, cen65_74, cen75_84, cen85plus),
      coef = list(cbcoef(m0_64)$coef, cbcoef(m65_74)$coef, cbcoef(m75_84)$coef,
        cbcoef(m85plus)$coef)),
    adapt = unname(exp(pred$allfit[ind] * 0.7))
  )
})

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(ref, digits = NA, auto_unbox = TRUE), out)
cat("wrote", out, "\n")
cat("change in AF (%) 2010-19 to 2090-99, RCP8.5:\n")
print(round(rbind(heat = env$afrel["2090-99", , "heat", "RCP8.5"],
  cold = env$afrel["2090-99", , "cold", "RCP8.5"]), 2))
