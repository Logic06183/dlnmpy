# Reference values for dlnmpy.attribution from Gasparrini's published R
# functions attrdl.R (2014 BMC MRM) and findmin.R (2017 Epidemiology),
# sourced from his GitHub repositories.
suppressPackageStartupMessages({ library(dlnm); library(splines); library(jsonlite) })
source("https://raw.githubusercontent.com/gasparrini/2014_gasparrini_BMCmrm_Rcodedata/master/attrdl.R")
source("https://raw.githubusercontent.com/gasparrini/2017_tobias_Epidem_Rcodedata/master/findmin.R")
set.seed(20260830)
outdir <- "tests/fixtures"
d <- chicagoNMMAPS
cb <- crossbasis(d$temp, lag = 21, argvar = list(fun = "bs", degree = 2, knots = quantile(d$temp, c(10, 75, 90) / 100)),
  arglag = list(knots = logknots(21, 3)))
m <- glm(death ~ cb + ns(time, 7 * 14) + dow, family = quasipoisson(), d, control = glm.control(epsilon = 1e-14))
cf <- coef(m); ind <- grep("cb[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}", names(cf))
coef <- unname(cf[ind]); vcov <- unname(vcov(m)[ind, ind])
# findmin on the 1st-99th percentile range, step 0.1
q <- quantile(d$temp, c(0.01, 0.99))
mmt <- findmin(cb, m, from = q[1], to = q[2], by = 0.1)
# simulation with an exported matrix of normals so Python can reproduce the draws exactly
nsim <- 200; k <- length(coef)
Z <- matrix(rnorm(k * nsim), nsim)
eig <- eigen(vcov)
coefsim <- coef + eig$vectors %*% diag(sqrt(eig$values), k) %*% t(Z)
minsim <- apply(coefsim, 2, function(ci) { p <- crosspred(cb, coef = ci, vcov = vcov, model.link = "log", from = q[1], to = q[2], by = 0.1); p$predvar[which.min(p$allfit)] })
# attributable measures with cen = mmt
red <- crossreduce(cb, m, cen = mmt)
res <- list(
  coef = coef, vcov = vcov, mmt = mmt, q = unname(q), Z = Z, coefsim = coefsim, minsim = unname(minsim),
  af_back = attrdl(d$temp, cb, d$death, m, cen = mmt),
  an_back = attrdl(d$temp, cb, d$death, m, type = "an", cen = mmt),
  af_forw = attrdl(d$temp, cb, d$death, m, dir = "forw", cen = mmt),
  an_forw = attrdl(d$temp, cb, d$death, m, type = "an", dir = "forw", cen = mmt),
  an_back_daily = attrdl(d$temp, cb, d$death, m, type = "an", cen = mmt, tot = FALSE)[1:200],
  af_forw_daily = attrdl(d$temp, cb, d$death, m, dir = "forw", cen = mmt, tot = FALSE)[1:200],
  an_cold = attrdl(d$temp, cb, d$death, m, type = "an", cen = mmt, range = c(-100, mmt)),
  an_heat = attrdl(d$temp, cb, d$death, m, type = "an", cen = mmt, range = c(mmt, 100)),
  an_extreme_heat = attrdl(d$temp, cb, d$death, m, type = "an", cen = mmt, range = c(quantile(d$temp, .975), 100)),
  af_reduced_forw = attrdl(d$temp, cb, d$death, coef = coef(red), vcov = vcov(red), model.link = "log", dir = "forw", cen = mmt),
  an_back_sim = attrdl(d$temp, cb, d$death, coef = coef, vcov = vcov, model.link = "log", type = "an", cen = mmt, sim = TRUE, nsim = nsim),
  red_coef = unname(coef(red)), red_vcov = unname(vcov(red)),
  # exposure-history matrix input (backward) and matrix of future cases (forward)
  an_back_matrix = attrdl(tsModel::Lag(d$temp, 0:21), cb, d$death, m, type = "an", cen = mmt),
  af_forw_casemat = attrdl(d$temp, cb, as.matrix(tsModel::Lag(d$death, -(0:21))), m, dir = "forw", cen = mmt)
)
# the simulated AN uses R's own rnorm stream inside attrdl; export the seed-independent version by re-running with coefsim
ani <- apply(coefsim, 2, function(ci) attrdl(d$temp, cb, d$death, coef = ci, vcov = vcov, model.link = "log", type = "an", cen = mmt))
res$an_back_sim_fixed <- unname(ani)
write_json(res, file.path(outdir, "attribution.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("mmt", mmt, " af_back", res$af_back, " an_back", res$an_back, "\n")
