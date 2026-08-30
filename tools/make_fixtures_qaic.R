# Reference QAIC values from R for dlnmpy.uncertainty.qaic.
#
# R's logLik() is NA for a quasipoisson glm, which is why the QAIC used
# throughout this literature (Peng, Dominici & Louis 2006; Gasparrini,
# Armstrong & Kenward 2010) evaluates the *Poisson* log-likelihood at the
# fitted values by hand:
#
#     QAIC = -2 * loglik(Poisson) + 2 * phi * k
#
# The grid below is chosen so that the ranking is informative: it spans
# under- and over-fitted lag/df combinations whose dispersions differ, so a
# criterion that mistakenly scales the log-likelihood by phi reorders it.
suppressPackageStartupMessages({ library(dlnm); library(splines); library(jsonlite) })
outdir <- "tests/fixtures"
d <- chicagoNMMAPS
d$nst <- ns(d$time, 7 * 14)

fqaic <- function(m) {
  loglik <- sum(dpois(m$y, m$fitted.values, log = TRUE))
  -2 * loglik + 2 * summary(m)$df[3] * summary(m)$dispersion
}

specs <- list(list(lag = 7, df_var = 3, df_lag = 3), list(lag = 14, df_var = 4, df_lag = 4),
              list(lag = 21, df_var = 4, df_lag = 4), list(lag = 21, df_var = 6, df_lag = 5),
              list(lag = 28, df_var = 5, df_lag = 4), list(lag = 3, df_var = 3, df_lag = 2))

fitone <- function(s, fam) {
  cb <- crossbasis(d$temp, lag = s$lag, argvar = list(fun = "ns", df = s$df_var),
                   arglag = list(fun = "ns", df = s$df_lag))
  # tight convergence so the comparison is not limited by R's default epsilon
  m <- glm(death ~ cb + nst + dow, family = fam, d, control = glm.control(epsilon = 1e-14))
  list(lag = s$lag, df_var = s$df_var, df_lag = s$df_lag,
       qaic = fqaic(m), loglik_pois = sum(dpois(m$y, m$fitted.values, log = TRUE)),
       dispersion = summary(m)$dispersion, rank = summary(m)$df[3], deviance = deviance(m))
}

quasi <- lapply(specs, fitone, fam = quasipoisson())
# with a known scale QAIC reduces to AIC, which R does report
mp <- glm(death ~ crossbasis(d$temp, lag = 21, argvar = list(fun = "ns", df = 4),
                             arglag = list(fun = "ns", df = 4)) + nst + dow,
          family = poisson(), d, control = glm.control(epsilon = 1e-14))
pois <- list(qaic = fqaic(mp), aic = AIC(mp), loglik_pois = sum(dpois(mp$y, mp$fitted.values, log = TRUE)),
             dispersion = summary(mp)$dispersion, rank = summary(mp)$df[3],
             lag = 21, df_var = 4, df_lag = 4)

out <- list(quasipoisson = quasi, poisson = pois,
            # R's ordering, best (lowest QAIC) first, as 0-based indices into `quasipoisson`
            order = order(sapply(quasi, `[[`, "qaic")) - 1L)
write_json(out, file.path(outdir, "qaic.json"), digits = NA, auto_unbox = TRUE, pretty = TRUE)
cat("QAIC fixtures written. R ranking (best first):", out$order, "\n")
cat("best spec: lag", quasi[[out$order[1] + 1]]$lag, "df_var", quasi[[out$order[1] + 1]]$df_var, "\n")
