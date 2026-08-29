# Reference values for dlnmpy.meta from the R package mixmeta (Gasparrini):
# a simulated 12-location two-stage analysis.
suppressPackageStartupMessages({ library(dlnm); library(splines); library(mixmeta); library(jsonlite) })
set.seed(20260831)
outdir <- "tests/fixtures"
n <- 2000; ncity <- 12
sim <- list(); meanT <- numeric(ncity)
for (k in seq_len(ncity)) {
  shift <- rnorm(1, 0, 3); meanT[k] <- 15 + shift
  tmean <- 15 + shift + 9 * sin(2 * pi * seq_len(n) / 365.25) + rnorm(n, 0, 3)
  slope <- 0.0025 * (1 + 0.05 * shift)               # heterogeneity linked to the city predictor
  ftrue <- function(x) slope * (x - 18)^2 * ifelse(x > 18, 1.6, 1)
  wtrue <- exp(-seq(0, 10) / 3); wtrue <- wtrue / sum(wtrue)
  Q <- tsModel::Lag(tmean, 0:10)
  eta <- rowSums(sweep(sapply(seq_len(11), function(l) ftrue(Q[, l])), 2, wtrue, "*"))
  mu <- exp(log(30) + 0.1 * sin(2 * pi * seq_len(n) / 365.25 + 1) + eta)
  y <- rpois(n, ifelse(is.na(mu), 30, mu)); y[is.na(mu)] <- NA
  sim[[k]] <- data.frame(city = k, time = seq_len(n), tmean = tmean, y = y)
}
write.csv(do.call(rbind, sim), file.path(outdir, "meta_sim.csv"), row.names = FALSE)

# stage 1: same knots for every city (on the pooled range) so coefficients are comparable
knots <- c(8, 15, 22); bk <- c(-5, 40)
ymat <- NULL; Slist <- list(); stage1 <- list()
for (k in seq_len(ncity)) {
  d <- sim[[k]]
  cb <- crossbasis(d$tmean, lag = 10, argvar = list(fun = "ns", knots = knots, Boundary.knots = bk), arglag = list(fun = "ns", df = 3))
  m <- glm(y ~ cb + ns(time, 6 * 3), family = quasipoisson(), d, control = glm.control(epsilon = 1e-14))
  r <- crossreduce(cb, m, cen = 18)
  ymat <- rbind(ymat, coef(r)); Slist[[k]] <- vcov(r)
  stage1[[k]] <- list(coef = unname(coef(r)), vcov = unname(vcov(r)))
}
S <- t(sapply(Slist, vechMat))
X <- cbind(1, meanT)

out <- list(y = unname(ymat), S_vech = unname(S), meanT = meanT, knots = knots, bk = bk, stage1 = stage1)
fits <- list()
for (meth in c("reml", "ml", "fixed")) {
  mm <- mixmeta(ymat ~ 1, S, method = meth)
  mr <- mixmeta(ymat ~ meanT, S, method = meth)
  for (nm in c("mm", "mr")) {
    o <- get(nm); q <- qtest(o); s <- summary(o)
    b <- blup(o, se = TRUE, pi = TRUE, vcov = TRUE, format = "list")
    pr <- predict(o, newdata = data.frame(meanT = c(10, 15, 20)), se = TRUE, ci = TRUE, vcov = TRUE, format = "list")
    fits[[paste(meth, nm, sep = "_")]] <- list(
      coef = unname(coef(o)), vcov = unname(vcov(o)), Psi = unname(if (meth == "fixed") matrix(0, ncol(ymat), ncol(ymat)) else o$Psi),
      logLik = as.numeric(logLik(o)), AIC = AIC(o), BIC = BIC(o), df_res = o$df.residual,
      Q = unname(q$Q), Qdf = unname(q$df), Qp = unname(q$pvalue), I2 = unname(s$i2stat),
      blup = lapply(b, function(x) list(blup = unname(x$blup), se = unname(x$se), pi_lb = unname(x$pi.lb), pi_ub = unname(x$pi.ub), vcov = unname(x$vcov))),
      pred = lapply(pr, function(x) list(fit = unname(x$fit), se = unname(x$se), vcov = unname(x$vcov))),
      converged = o$converged)
  }
}
# diagonal and identity structures
fits$reml_diag <- { o <- mixmeta(ymat ~ 1, S, method = "reml", bscov = "diag"); list(coef = unname(coef(o)), Psi = unname(o$Psi), logLik = as.numeric(logLik(o))) }
fits$reml_id <- { o <- mixmeta(ymat ~ 1, S, method = "reml", bscov = "id"); list(coef = unname(coef(o)), Psi = unname(o$Psi), logLik = as.numeric(logLik(o))) }
# univariate case
fits$uni_reml <- { o <- mixmeta(ymat[, 1] ~ 1, S[, 1], method = "reml"); list(coef = unname(coef(o)), vcov = unname(vcov(o)), Psi = unname(o$Psi), logLik = as.numeric(logLik(o)), I2 = unname(summary(o)$i2stat), Q = unname(qtest(o)$Q)) }
# pooled curve from the REML meta-analysis, predicted with a onebasis
mm <- mixmeta(ymat ~ 1, S, method = "reml")
ob <- onebasis(seq(-5, 40, 1), fun = "ns", knots = knots, Boundary.knots = bk)
cp <- crosspred(ob, coef = coef(mm), vcov = vcov(mm), model.link = "log", at = seq(-5, 40, 1), cen = 18)
out$pooled_curve <- list(at = cp$predvar, allRRfit = unname(cp$allRRfit), allRRlow = unname(cp$allRRlow), allRRhigh = unname(cp$allRRhigh))
out$fits <- fits
write_json(out, file.path(outdir, "meta.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("done; REML Psi diag:", diag(mm$Psi), "\n")
