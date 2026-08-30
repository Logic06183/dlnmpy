# Reference values for dlnmpy.penalized from mgcv::gam with paraPen penalties
# (the penalised DLNM framework of Gasparrini et al. 2017).
suppressPackageStartupMessages({ library(dlnm); library(splines); library(mgcv); library(jsonlite) })
outdir <- "tests/fixtures"
d <- chicagoNMMAPS
d$nst <- ns(d$time, 7 * 14)
out <- list()
gamout <- function(m, cbname) {
  cf <- coef(m); ind <- grep(paste0("^", cbname), names(cf))
  list(coef_all = unname(cf), names = names(cf), coef = unname(cf[ind]), Vp = unname(m$Vp[ind, ind]),
    sp = unname(m$sp), scale = m$sig2, edf = unname(m$edf), edf_cb = sum(m$edf[ind]), score = as.numeric(m$gcv.ubre),
    deviance = deviance(m), converged = m$converged)
}
# A. ps/ps cross-basis, Poisson (known scale) and quasi-Poisson (scale estimated)
cb7 <- crossbasis(d$temp, lag = 10, argvar = list(fun = "ps", df = 8), arglag = list(fun = "ps", df = 5))
pen7 <- cbPen(cb7)
mA <- gam(death ~ cb7 + nst + dow, family = poisson(), data = d, paraPen = list(cb7 = pen7), method = "REML")
mAq <- gam(death ~ cb7 + nst + dow, family = quasipoisson(), data = d, paraPen = list(cb7 = pen7), method = "REML")
mAml <- gam(death ~ cb7 + nst + dow, family = poisson(), data = d, paraPen = list(cb7 = pen7), method = "ML")
# criterion at fixed smoothing parameters (checks the REML formula independently of the optimiser)
spfix <- c(10, 100)
pen7fix <- pen7; pen7fix$sp <- spfix
mAfix <- gam(death ~ cb7 + nst + dow, family = poisson(), data = d, paraPen = list(cb7 = pen7fix), method = "REML")
mAqfix <- gam(death ~ cb7 + nst + dow, family = quasipoisson(), data = d, paraPen = list(cb7 = pen7fix), method = "REML")
pA <- crosspred(cb7, mA, cen = 20, by = 1)
out$ps <- list(pen = list(Svar = unname(pen7$Svar), Slag = unname(pen7$Slag)),
  poisson = gamout(mA, "cb7"), quasi = gamout(mAq, "cb7"), ml = gamout(mAml, "cb7"),
  poisson_fixed = c(gamout(mAfix, "cb7"), list(spfix = spfix)), quasi_fixed = c(gamout(mAqfix, "cb7"), list(spfix = spfix)),
  pred = list(predvar = pA$predvar, allfit = unname(pA$allfit), allse = unname(pA$allse), matfit = unname(pA$matfit)))
# B. cr/cr cross-basis, quasi-Poisson
cb9 <- crossbasis(d$temp, lag = 21, argvar = list(fun = "cr", df = 6), arglag = list(fun = "cr", df = 5))
pen9 <- cbPen(cb9)
mB <- gam(death ~ cb9 + nst + dow, family = quasipoisson(), data = d, paraPen = list(cb9 = pen9), method = "REML")
out$cr <- list(quasi = gamout(mB, "cb9"))
# C. additional lag penalty (ridge on lag coefficients) as in the 2017 paper, Poisson
pen7b <- cbPen(cb7, addSlag = rep(1, 5))
mC <- gam(death ~ cb7 + nst + dow, family = poisson(), data = d, paraPen = list(cb7 = pen7b), method = "REML")
out$ps_addslag <- list(pen = list(Svar = unname(pen7b$Svar), Slag = unname(pen7b$Slag), Slag2 = unname(pen7b$Slag2)), poisson = gamout(mC, "cb7"))
write_json(out, file.path(outdir, "penalized.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE, na = "null")
cat("poisson sp", mA$sp, " quasi sp", mAq$sp, "scale", mAq$sig2, " cr sp", mB$sp, "\n")
