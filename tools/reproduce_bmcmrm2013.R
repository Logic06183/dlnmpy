################################################################################
# Reference values for examples/bmcmrm2013_vs_r.py
#
#   Gasparrini A, Armstrong B. Reducing and meta-analysing estimates from
#   distributed lag non-linear models. BMC Medical Research Methodology
#   2013;13:1. doi:10.1186/1471-2288-13-1
#
# Runs the authors' own code (01.prep.R to 06.metareg.R, unmodified) from a
# clone of https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata
# pinned to commit d542c9314901b00d89dae6d5000a599b7700d6da, on its England &
# Wales data, and writes every stage intermediate to
# tests/reference/bmcmrm2013.json:
#
#   first stage   the five reductions per region (overall cumulative for the
#                 three lag models; predictor-specific at 22 C and 0 C for
#                 the main one), coefficients and covariance; QAIC per region
#   second stage  the five pooled fits (coef, vcov, Psi, logLik, Q test, I2),
#                 the fixed-effects fit at 22 C, the three meta-regressions
#                 on latitude with their predictions at the 25th and 75th
#                 percentiles of latitude and Wald tests
#   predictions   the pooled curves (cpall, cpall2, cpall3, cphot, cpcold,
#                 cphot2 and the six latitude curves), the minimum-mortality
#                 temperature, RR at 22 C and 0 C
#
# Three chains are recorded.
#
#   "published": the authors' code as written: glm() at its default
#     epsilon = 1e-8, mvmeta() for the second stage.
#   "converged" (not in the authors' code): 03.firststage.R sourced again,
#     unchanged, with glm() masked by a wrapper that adds
#     glm.control(epsilon = 1e-14); the second stage then fitted with
#     mixmeta(), the successor of mvmeta that mixmetapy ports, at its
#     default control, and the predictions repeated from those fits.
#   "reltol14": the converged first stage with mixmeta() at
#     control = list(reltol = 1e-14), its optimum to within machine
#     precision (second stage only).
#
# Usage, from the repository root:
#
#   git clone https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata /tmp/repro/bmcmrm2013
#   git -C /tmp/repro/bmcmrm2013 checkout d542c9314901b00d89dae6d5000a599b7700d6da
#   Rscript tools/reproduce_bmcmrm2013.R [/tmp/repro/bmcmrm2013]
#
# Needs dlnm, mvmeta, mixmeta, jsonlite. About a minute. The authors'
# figures are written to a temporary directory, not into their clone.
################################################################################

SHA <- "d542c9314901b00d89dae6d5000a599b7700d6da"

args <- commandArgs(trailingOnly = TRUE)
src <- normalizePath(if (length(args)) args[1] else "/tmp/repro/bmcmrm2013")
out <- file.path(getwd(), "tests", "reference", "bmcmrm2013.json")
if (!file.exists(file.path(getwd(), "src", "dlnmpy")))
  stop("run from the dlnmpy repository root")

head <- system2("git", c("-C", shQuote(src), "rev-parse", "HEAD"), stdout = TRUE)
if (!identical(head, SHA)) stop("authors' repository is at ", head, ", expected ", SHA)
dirty <- system2("git", c("-C", shQuote(src), "status", "--porcelain", "--untracked-files=no"),
  stdout = TRUE)
if (length(dirty)) stop("authors' scripts are modified:\n", paste(dirty, collapse = "\n"))

suppressPackageStartupMessages({library(jsonlite); library(mixmeta)})

figdir <- tempfile("bmcmrm2013_figures")
dir.create(figdir)
runscript <- function(f, wd = figdir) {
  owd <- setwd(wd); on.exit(setwd(owd))
  invisible(capture.output(source(file.path(src, f), echo = FALSE)))
}

################################################################################
# THE AUTHORS' CODE, AS PUBLISHED

runscript("01.prep.R", wd = src)          # reads regEngWales.csv from the clone
for (f in c("02.singleregion.R", "03.firststage.R", "04.secondstage.R",
  "05.mainresults.R", "06.metareg.R")) runscript(f)
graphics.off()

RED <- c(all = "all", all2 = "all2", all3 = "all3", hot = "hot", cold = "cold")
reductions <- function() lapply(RED, function(k) list(
  coef = unname(get(paste0("y", k))),
  vcov = unname(lapply(get(paste0("S", k)), unname))))

qt <- function(fit) {
  q <- qtest(fit)
  list(Q = unname(q$Q), df = unname(q$df), pvalue = unname(q$pvalue),
    I2 = unname(((q$Q - q$df) / q$Q)[1] * 100))
}

# flat coefficients and vcov in predictor-major order (mixmeta's and
# mixmetapy's); mvmeta's are outcome-major. Same for intercept-only fits.
metafit <- function(fit) {
  p <- fit$dim$p; k <- fit$dim$k
  v <- vcov(fit)
  b <- coef(fit)
  if (inherits(fit, "mvmeta")) {
    B <- matrix(fit$coefficients, p, k)
    perm <- as.vector(t(matrix(seq(p * k), p, k)))
    v <- v[perm, perm]
    b <- b[perm]
  } else {
    B <- matrix(b, p, k, byrow = TRUE)
  }
  if (!identical(unname(as.vector(t(B))), unname(b))) stop("unexpected coefficient order")
  Psi <- if (is.null(fit$Psi)) matrix(0, k, k) else unname(fit$Psi)
  c(list(package = class(fit)[1], method = fit$method, coef = unname(B),
    vcov = unname(v), Psi = Psi, logLik = as.numeric(logLik(fit)),
    converged = if (is.null(fit$converged)) TRUE else fit$converged), list(qtest = qt(fit)))
}

pred <- function(cp) list(predvar = unname(cp$predvar), allRRfit = unname(cp$allRRfit),
  allRRlow = unname(cp$allRRlow), allRRhigh = unname(cp$allRRhigh))

latpred <- function(p) lapply(unname(p), function(x) list(fit = unname(x$fit),
  vcov = unname(x$vcov)))

# everything from the second stage on, from the globals the authors' scripts
# leave behind (or that the converged chain below recreates)
secondstage <- function() list(
  pooled = lapply(c(all = "mvall", all2 = "mvall2", all3 = "mvall3", hot = "mvhot",
    cold = "mvcold"), function(n) metafit(get(n))),
  fixed_hot = metafit(mvhot2),
  metareg = list(all = metafit(mvalllat), hot = metafit(mvhotlat), cold = metafit(mvcoldlat)),
  wald_lat = c(all = fwald(mvalllat, "lat"), hot = fwald(mvhotlat, "lat"),
    cold = fwald(mvcoldlat, "lat")),
  latpred = list(val = unname(val), all = latpred(predall), hot = latpred(predhot),
    cold = latpred(predcold)),
  curves = list(cpall = pred(cpall), cpall2 = pred(cpall2), cpall3 = pred(cpall3),
    cphot = pred(cphot), cpcold = pred(cpcold), cphot2 = pred(cphot2),
    cpalllat25 = pred(cpalllat25), cpalllat75 = pred(cpalllat75),
    cphotlat25 = pred(cphotlat25), cphotlat75 = pred(cphotlat75),
    cpcoldlat25 = pred(cpcoldlat25), cpcoldlat75 = pred(cpcoldlat75)),
  mmt = cpall$predvar[which.min(cpall$allRRfit)],
  rr22 = unname(with(cpall, cbind(allRRfit, allRRlow, allRRhigh)["22", ])),
  rr0 = unname(with(cpall, cbind(allRRfit, allRRlow, allRRhigh)["0", ])))

published <- list(
  first_stage = list(reductions = reductions(),
    qaic = list(m1 = qaic, m2 = qaic2, m3 = qaic3),
    qaic_total = c(m1 = sum(qaic), m2 = sum(qaic2), m3 = sum(qaic3))),
  second_stage = secondstage())

settings <- list(regions = regions, bound = unname(bound), varknots = unname(varknots),
  lagknots = unname(lagknots), lat = lat, cen = 17,
  arglag_cb = lapply(attr(cb, "arglag"), unname),
  mmt_percentile = round(sum(regEngWales$tmean < 17.1) / nrow(regEngWales) * 100, 1))

################################################################################
# FIRST STAGE RUN TO CONVERGENCE: 03.firststage.R again, unchanged, with
# glm() masked so that every fit gets glm.control(epsilon = 1e-14)

glm <- function(...) stats::glm(..., control = glm.control(epsilon = 1e-14, maxit = 100))
runscript("03.firststage.R")
rm(glm)
converged_first <- list(epsilon = 1e-14, reductions = reductions(),
  qaic = list(m1 = qaic, m2 = qaic2, m3 = qaic3),
  qaic_total = c(m1 = sum(qaic), m2 = sum(qaic2), m3 = sum(qaic3)))

################################################################################
# SECOND STAGE WITH mixmeta, on the converged reductions
#
# The calls of 04.secondstage.R, 05.mainresults.R and 06.metareg.R with
# mvmeta() replaced by mixmeta(); bvar and blag are unchanged.

mixchain <- function(control = list()) {
  e <- globalenv()
  mm <- function(y, S, f = ~1) {
    y <- get(y, e); S <- get(S, e)
    form <- update(as.formula(paste("y", "~ .")), f)
    mixmeta(form, S, method = "reml", control = control)
  }
  assign("mvall", mm("yall", "Sall"), e)
  assign("mvall2", mm("yall2", "Sall2"), e)
  assign("mvall3", mm("yall3", "Sall3"), e)
  assign("mvhot", mm("yhot", "Shot"), e)
  assign("mvcold", mm("ycold", "Scold"), e)
  assign("mvhot2", mixmeta(yhot ~ 1, Shot, method = "fixed"), e)
  assign("mvalllat", mm("yall", "Sall", ~ lat), e)
  assign("mvhotlat", mm("yhot", "Shot", ~ lat), e)
  assign("mvcoldlat", mm("ycold", "Scold", ~ lat), e)
}

mixchain()
# the prediction calls of the authors' scripts, verbatim
cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
  model.link="log",by=0.1,from=bound[1],to=bound[2],cen=17)
cpall2 <- crosspred(bvar,coef=coef(mvall2),vcov=vcov(mvall2),
  model.link="log",by=0.1,from=bound[1],to=bound[2],cen=17)
cpall3 <- crosspred(bvar,coef=coef(mvall3),vcov=vcov(mvall3),
  model.link="log",by=0.1,from=bound[1],to=bound[2],cen=17)
cphot <- crosspred(blag,coef=coef(mvhot),vcov=vcov(mvhot),
  model.link="log",at=0:210/10,cen=17)
cpcold <- crosspred(blag,coef=coef(mvcold),vcov=vcov(mvcold),
  model.link="log",at=0:210/10,cen=17)
cphot2 <- crosspred(blag,coef=coef(mvhot2),vcov=vcov(mvhot2),
  model.link="log",at=0:210/10)
predall <- predict(mvalllat,data.frame(lat=val),vcov=T)
cpalllat25 <- crosspred(bvar,coef=predall[[1]]$fit,vcov=predall[[1]]$vcov,
  model.link="log",by=0.2,cen=17)
cpalllat75 <- crosspred(bvar,coef=predall[[2]]$fit,vcov=predall[[2]]$vcov,
  model.link="log",by=0.2,cen=17)
predhot <- predict(mvhotlat,data.frame(lat=val),vcov=T)
cphotlat25 <- crosspred(blag,coef=predhot[[1]]$fit,vcov=predhot[[1]]$vcov,
  model.link="log",by=0.2)
cphotlat75 <- crosspred(blag,coef=predhot[[2]]$fit,vcov=predhot[[2]]$vcov,
  model.link="log",by=0.2)
predcold <- predict(mvcoldlat,data.frame(lat=val),vcov=T)
cpcoldlat25 <- crosspred(blag,coef=predcold[[1]]$fit,vcov=predcold[[1]]$vcov,
  model.link="log",by=0.2)
cpcoldlat75 <- crosspred(blag,coef=predcold[[2]]$fit,vcov=predcold[[2]]$vcov,
  model.link="log",by=0.2)
converged <- list(first_stage = converged_first, second_stage = secondstage())

# mvmeta and mixmeta report REML log-likelihoods that differ by a constant
# depending only on the dimensions; record it for each fit, both packages run
# to their optimum so that the difference is the constant alone
llconst <- local({
  tc <- list(reltol = 1e-14, maxiter = 1000)
  one <- function(y, f, method = "reml") {
    Y <- get(y); S <- get(sub("^y", "S", y))
    form <- update(Y ~ 1, f)
    as.numeric(logLik(mvmeta(form, S, method = method, control = tc)) -
      logLik(mixmeta(form, S, method = method, control = tc)))
  }
  list(pooled = as.list(sapply(c(all = "yall", all2 = "yall2", all3 = "yall3", hot = "yhot",
      cold = "ycold"), one, f = ~ 1)),
    fixed_hot = one("yhot", ~ 1, "fixed"),
    metareg = as.list(sapply(c(all = "yall", hot = "yhot", cold = "ycold"), one, f = ~ lat)))
})

mixchain(control = list(reltol = 1e-14, maxiter = 1000))
tight <- list(
  pooled = lapply(c(all = "mvall", all2 = "mvall2", all3 = "mvall3", hot = "mvhot",
    cold = "mvcold"), function(n) metafit(get(n))),
  metareg = list(all = metafit(mvalllat), hot = metafit(mvhotlat), cold = metafit(mvcoldlat)))

################################################################################

res <- list(
  versions = list(
    R = R.version.string,
    dlnm = as.character(packageVersion("dlnm")),
    mvmeta = as.character(packageVersion("mvmeta")),
    mixmeta = as.character(packageVersion("mixmeta")),
    authors_repo = "https://github.com/gasparrini/2013_gasparrini_BMCmrm_Rcodedata",
    authors_commit = SHA,
    regEngWales_csv_md5 = unname(tools::md5sum(file.path(src, "regEngWales.csv")))),
  settings = settings,
  published = published,
  converged = converged,
  reltol14 = tight,
  loglik_constant_mvmeta_minus_mixmeta = llconst)

dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(res, digits = NA, auto_unbox = TRUE, null = "null", na = "null",
  pretty = FALSE), out)
unlink(figdir, recursive = TRUE)
cat("wrote", out, "\n")
for (ch in c("published", "converged")) {
  s2 <- res[[ch]]$second_stage
  cat(sprintf("%-9s MMT %.1f  RR22 %.3f (%.3f-%.3f)  RR0 %.3f (%.3f-%.3f)  I2 %.1f/%.1f/%.1f  QAIC %.1f %.1f %.1f\n",
    ch, s2$mmt, s2$rr22[1], s2$rr22[2], s2$rr22[3], s2$rr0[1], s2$rr0[2], s2$rr0[3],
    s2$pooled$all$qtest$I2, s2$pooled$hot$qtest$I2, s2$pooled$cold$qtest$I2,
    res[[ch]]$first_stage$qaic_total[1], res[[ch]]$first_stage$qaic_total[2],
    res[[ch]]$first_stage$qaic_total[3]))
}
cat("logLik, pooled overall: mvmeta", published$second_stage$pooled$all$logLik,
  " mixmeta", converged$second_stage$pooled$all$logLik, " reltol 1e-14", tight$pooled$all$logLik, "\n")
