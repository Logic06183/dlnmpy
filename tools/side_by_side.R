# Second, independent validation: full analyses run end to end in R (including
# model fitting), on data and specifications NOT used for the unit-test
# fixtures. tools/side_by_side.py repeats each analysis in Python and compares.
#
#   Rscript tools/side_by_side.R && python tools/side_by_side.py

suppressPackageStartupMessages({ library(dlnm); library(splines); library(survival); library(jsonlite) })
set.seed(20260829)
# R glm stops at a looser tolerance (epsilon = 1e-8) than statsmodels; tighten it
# so that any remaining difference is algorithmic rather than convergence.
ctrl <- glm.control(epsilon = 1e-14, maxit = 100)
outdir <- "tests/side_by_side"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
res <- list()
predout <- function(p) list(predvar = p$predvar, cen = p$cen, lag = p$lag, bylag = p$bylag,
  coef = unname(p$coefficients), vcov = unname(p$vcov), matfit = unname(p$matfit),
  matse = unname(p$matse), allfit = unname(p$allfit), allse = unname(p$allse),
  alllow = unname(if (is.null(p$alllow)) log(p$allRRlow) else p$alllow),
  allhigh = unname(if (is.null(p$allhigh)) log(p$allRRhigh) else p$allhigh),
  cumfit = if (is.null(p$cumfit)) NULL else unname(p$cumfit), ci.level = p$ci.level, link = p$model.link)

# A. drug trial: lm, Gaussian, matrix of exposure histories, lag 0-27 -----------
Qdrug <- as.matrix(drug[, rep(7:4, each = 7)])
cbdrug <- crossbasis(Qdrug, lag = 27, argvar = list("lin"), arglag = list(fun = "ns", knots = c(9, 18)))
mdrug <- lm(out ~ cbdrug + sex, drug)
pdrug <- crosspred(cbdrug, mdrug, at = 0:20 * 5, ci.level = 0.90)
res$drug <- list(pred = predout(pdrug), coef_all = unname(coef(mdrug)), names = names(coef(mdrug)),
  sigma = summary(mdrug)$sigma)
write.csv(Qdrug, file.path(outdir, "Qdrug.csv"), row.names = FALSE)

# B. nested case-control: conditional logistic regression, lag 3-40 ------------
Qnest <- t(apply(nested, 1, function(sub) exphist(rep(c(0, 0, 0, sub[5:14]), each = 5), sub["age"], lag = c(3, 40))))
cbnest <- crossbasis(Qnest, lag = c(3, 40), argvar = list("bs", degree = 2, df = 3),
  arglag = list(fun = "ns", knots = c(10, 30), intercept = FALSE))
mnest <- clogit(case ~ cbnest + strata(riskset), nested)
pnest <- crosspred(cbnest, mnest, cen = 0, at = 0:20 * 5)
res$nested <- list(pred = predout(pnest), coef_all = unname(coef(mnest)))
write.csv(Qnest, file.path(outdir, "Qnest.csv"), row.names = FALSE)

# C. Chicago, logistic regression on a binary outcome, thr var + integer lag ---
chi <- chicagoNMMAPS
chi$high <- as.integer(chi$death > median(chi$death))
cbC <- crossbasis(chi$temp, lag = 7, argvar = list(fun = "thr", thr = c(5, 25)), arglag = list(fun = "integer"))
mC <- glm(high ~ cbC + ns(time, 7 * 14) + dow, family = binomial(), chi, control = ctrl)
pC <- crosspred(cbC, mC, at = c(-20, -10, 0, 10, 20, 30, 33), cumul = TRUE)
res$chicago_logit <- list(pred = predout(pC), coef_all = unname(coef(mC)), names = names(coef(mC)))

# D. Chicago, quasi-Poisson, poly var basis, strata lag, matrix 'at', 80% CI ----
cbD <- crossbasis(chi$temp, lag = 10, argvar = list(fun = "poly", degree = 3), arglag = list(fun = "strata", breaks = c(1, 4)))
mD <- glm(death ~ cbD + ns(time, 7 * 14) + dow, family = quasipoisson(), chi, control = ctrl)
pD <- crosspred(cbD, mD, cen = 20, by = 2, ci.level = 0.80)
atmat <- exphist(chi$temp, times = c(200, 2000, 4000), lag = 10)
pDm <- crosspred(cbD, mD, cen = 20, at = atmat)
res$chicago_poly <- list(pred = predout(pD), pred_mat = predout(pDm), atmat = unclass(atmat),
  coef_all = unname(coef(mD)), dispersion = summary(mD)$dispersion,
  red = { r <- crossreduce(cbD, mD, cen = 20, by = 2, ci.level = 0.8); list(coef = unname(coef(r)), vcov = unname(vcov(r)), fit = unname(r$fit), se = unname(r$se)) })

# E. Simulated data with known truth: two-stage design, 4 cities --------------
n <- 3000; ncity <- 4; sim <- list(); truth <- list()
for (k in seq_len(ncity)) {
  tmean <- 15 + 10 * sin(2 * pi * seq_len(n) / 365.25) + rnorm(n, 0, 3) + (k - 2)
  # true exposure-lag-response: U-shape centred at 18, decaying over 10 lags
  ftrue <- function(x) 0.0025 * (x - 18)^2 * ifelse(x > 18, 1.6, 1)
  wtrue <- exp(-seq(0, 10) / 3); wtrue <- wtrue / sum(wtrue)
  Q <- tsModel::Lag(tmean, 0:10)
  eta <- rowSums(sweep(sapply(seq_len(11), function(l) ftrue(Q[, l])), 2, wtrue, "*"))
  mu <- exp(log(40) + 0.1 * sin(2 * pi * seq_len(n) / 365.25 + 1) + eta)
  y <- rpois(n, ifelse(is.na(mu), 40, mu)); y[is.na(mu)] <- NA
  sim[[k]] <- data.frame(city = k, time = seq_len(n), tmean = tmean, y = y)
}
simdf <- do.call(rbind, sim)
write.csv(simdf, file.path(outdir, "sim.csv"), row.names = FALSE)
stage1 <- list()
for (k in seq_len(ncity)) {
  d <- sim[[k]]
  cb <- crossbasis(d$tmean, lag = 10, argvar = list(fun = "ns", knots = c(10, 18, 25)), arglag = list(fun = "ns", df = 3))
  m <- glm(y ~ cb + ns(time, 8 * 3), family = quasipoisson(), d, control = ctrl)
  p <- crosspred(cb, m, cen = 18, at = seq(0, 30, 2))
  r <- crossreduce(cb, m, cen = 18, at = seq(0, 30, 2))
  stage1[[k]] <- list(pred = predout(p), red_coef = unname(coef(r)), red_vcov = unname(vcov(r)),
    knots_var = attr(cb, "argvar")$knots, dispersion = summary(m)$dispersion)
}
res$sim <- list(stage1 = stage1, truth_x = seq(0, 30, 2), truth_all = sapply(seq(0, 30, 2), function(x) 0.0025 * (x - 18)^2 * ifelse(x > 18, 1.6, 1)))

write_json(res, file.path(outdir, "r_results.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
cat("written to", outdir, "\n")
