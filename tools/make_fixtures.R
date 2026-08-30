# Generate reference fixtures from the R package dlnm.
#
# Run from the repository root:
#   Rscript tools/make_fixtures.R
#
# Every number the Python test-suite checks against comes from here. The
# fixtures are language-neutral (CSV + JSON) so a port in any language can be
# validated against the same reference values.

suppressPackageStartupMessages({
  library(dlnm)
  library(splines)
  library(jsonlite)
})
set.seed(13041975)
outdir <- "tests/fixtures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

wjson <- function(obj, name)
  write_json(obj, file.path(outdir, paste0(name, ".json")), digits = NA,
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
wcsv <- function(mat, name)
  write.csv(as.data.frame(unclass(mat)), file.path(outdir, paste0(name, ".csv")),
    row.names = FALSE)

# strip classes/dimnames but keep the useful attributes as a plain list
attrs <- function(b, keep) {
  a <- attributes(b)
  a <- a[names(a) %in% keep]
  lapply(a, function(v) if (is.matrix(v)) unclass(v) else v)
}

# ----------------------------------------------------------------------------
# 0. Data
# ----------------------------------------------------------------------------
write.csv(chicagoNMMAPS, file.path(outdir, "chicagoNMMAPS.csv"), row.names = FALSE)
data(drug); data(nested)
write.csv(drug, file.path(outdir, "drug.csv"), row.names = FALSE)
write.csv(nested, file.path(outdir, "nested.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# 1. Basis functions (onebasis) on a fixed vector, including NA and values
#    outside the boundary knots at prediction time.
# ----------------------------------------------------------------------------
x <- chicagoNMMAPS$temp[1:400]
xpred <- c(-30, -26.6, -10, 0, 5.5, 12, 21, 27, 33, 36, 40)

specs <- list(
  lin            = list(fun = "lin"),
  lin_int        = list(fun = "lin", intercept = TRUE),
  poly3          = list(fun = "poly", degree = 3),
  poly4_int      = list(fun = "poly", degree = 4, intercept = TRUE),
  poly2_scale    = list(fun = "poly", degree = 2, scale = 10),
  strata_df3     = list(fun = "strata", df = 3),
  strata_breaks  = list(fun = "strata", breaks = c(0, 10, 20)),
  strata_ref2    = list(fun = "strata", breaks = c(0, 10, 20), ref = 2),
  strata_int     = list(fun = "strata", breaks = c(0, 10, 20), intercept = TRUE),
  strata_int_ref0= list(fun = "strata", breaks = c(0, 10, 20), intercept = TRUE, ref = 0),
  strata_df1_int = list(fun = "strata", df = 1, intercept = TRUE),
  thr_h          = list(fun = "thr", thr.value = 25),
  thr_l          = list(fun = "thr", thr.value = 5, side = "l"),
  thr_d          = list(fun = "thr", thr.value = c(10, 25)),
  thr_default    = list(fun = "thr"),
  thr_int        = list(fun = "thr", thr.value = 25, intercept = TRUE),
  integer_       = list(fun = "integer"),
  integer_int    = list(fun = "integer", intercept = TRUE),
  ns_df5         = list(fun = "ns", df = 5),
  ns_df5_int     = list(fun = "ns", df = 5, intercept = TRUE),
  ns_knots       = list(fun = "ns", knots = c(-5, 10, 22)),
  ns_knots_bk    = list(fun = "ns", knots = c(-5, 10, 22), Boundary.knots = c(-20, 30)),
  ns_df2         = list(fun = "ns", df = 2),
  ns_df1         = list(fun = "ns", df = 1),
  bs_df5_deg2    = list(fun = "bs", df = 5, degree = 2),
  bs_df6         = list(fun = "bs", df = 6),
  bs_knots_int   = list(fun = "bs", knots = c(0, 15), intercept = TRUE),
  bs_deg1        = list(fun = "bs", degree = 1, df = 3),
  ps_df10        = list(fun = "ps", df = 10),
  ps_df8_int     = list(fun = "ps", df = 8, intercept = TRUE),
  ps_deg2_diff1  = list(fun = "ps", df = 7, degree = 2, diff = 1),
  ps_knots2      = list(fun = "ps", df = 9, knots = c(-20, 30)),
  cr_df5         = list(fun = "cr", df = 5),
  cr_df6_int     = list(fun = "cr", df = 6, intercept = TRUE),
  cr_knots       = list(fun = "cr", knots = c(-10, 0, 10, 20, 28)),
  cr_fx          = list(fun = "cr", df = 4, fx = TRUE)
)
keep <- c("fun", "range", "df", "knots", "Boundary.knots", "degree", "intercept",
  "scale", "breaks", "ref", "thr.value", "side", "values", "fx", "S", "diff")

onebasis_fx <- list()
for (nm in names(specs)) {
  args <- specs[[nm]]
  xx <- if (args$fun == "integer") round(x) else x
  b <- do.call(onebasis, c(list(x = xx), args))
  # prediction-stage reconstruction, exactly as mkXpred does it
  ind <- match(c("fun", names(formals(get(attr(b, "fun"), envir = asNamespace("dlnm"))))), names(attributes(b)), nomatch = 0)
  xp <- if (args$fun == "integer") round(xpred) else xpred
  bp <- do.call(onebasis, c(list(x = xp), attributes(b)[ind]))
  onebasis_fx[[nm]] <- list(spec = args, x = xx, basis = unclass(b)[, , drop = FALSE],
    attributes = attrs(b, keep), xpred = xp, basis_pred = unclass(bp)[, , drop = FALSE])
}
wjson(onebasis_fx, "onebasis")

# NA handling for splines
xna <- x; xna[c(3, 50, 100)] <- NA
b <- onebasis(xna, "ns", df = 4)
wjson(list(x = xna, basis = unclass(b)[, , drop = FALSE], knots = attr(b, "knots"),
  Boundary.knots = attr(b, "Boundary.knots")), "onebasis_na")

# ----------------------------------------------------------------------------
# 2. Knot helpers and pretty()
# ----------------------------------------------------------------------------
wjson(list(
  logknots_30_3        = logknots(30, 3),
  logknots_30_df4      = logknots(30, df = 4),
  logknots_30_df5_bs   = logknots(30, fun = "bs", df = 5),
  logknots_range       = logknots(c(3, 21), 2),
  logknots_vec         = logknots(0:15, df = 3),
  logknots_strata_df3  = logknots(10, fun = "strata", df = 3),
  equalknots_temp_bs   = equalknots(chicagoNMMAPS$temp, fun = "bs", df = 5, degree = 2),
  equalknots_temp_ns   = equalknots(chicagoNMMAPS$temp, df = 4),
  equalknots_nk        = equalknots(chicagoNMMAPS$temp, nk = 3),
  pretty_cases = lapply(list(
      c(-26.6, 33.3), c(0, 1), c(0, 356.2), c(1.3, 1.7), c(-3, 3), c(0, 0),
      c(5, 5), c(-0.0001, 0.0001), c(100, 1e6), c(-26.6, 33.3, 50), c(0, 20, 9)),
    function(v) list(args = v, result = if (length(v) == 3) pretty(v[1:2], n = v[3]) else pretty(v)))
), "helpers")

# mkat behaviour (at/from/to/by) through crosspred on a tiny model
cbx <- crossbasis(chicagoNMMAPS$temp, lag = 3, argvar = list(df = 3), arglag = list(df = 2))
mfit <- glm(death ~ cbx + ns(time, 7 * 14), family = quasipoisson(), chicagoNMMAPS)
mk <- list(
  default    = crosspred(cbx, mfit, cen = 20)$predvar,
  by1        = crosspred(cbx, mfit, by = 1, cen = 20)$predvar,
  by2_5      = crosspred(cbx, mfit, by = 2.5, cen = 20)$predvar,
  from_to    = crosspred(cbx, mfit, from = -10, to = 30, cen = 20)$predvar,
  from_to_by = crosspred(cbx, mfit, from = -10.5, to = 30, by = 3, cen = 20)$predvar,
  at_unsorted= crosspred(cbx, mfit, at = c(30, 0, 0, -20), cen = 20)$predvar,
  auto_cen   = crosspred(cbx, mfit)$cen
)
wjson(mk, "mkat")

# ----------------------------------------------------------------------------
# 3. Cross-bases from the dlnmTS vignette (+ a few extra) with model fits
# ----------------------------------------------------------------------------
cbmeta <- function(cb) {
  a <- attributes(cb)
  list(df = a$df, range = a$range, lag = a$lag,
    argvar = lapply(a$argvar, function(v) if (is.matrix(v)) unclass(v) else v),
    arglag = lapply(a$arglag, function(v) if (is.matrix(v)) unclass(v) else v),
    group = a$group)
}
predmeta <- function(p) {
  out <- list(predvar = p$predvar, cen = p$cen, lag = p$lag, bylag = p$bylag,
    coef = unname(p$coefficients), vcov = unname(p$vcov),
    matfit = unname(p$matfit), matse = unname(p$matse),
    allfit = unname(p$allfit), allse = unname(p$allse),
    ci.level = p$ci.level, model.link = p$model.link)
  if (!is.null(p$cumfit)) { out$cumfit <- unname(p$cumfit); out$cumse <- unname(p$cumse) }
  if (!is.null(p$allRRfit)) {
    out$allRRfit <- unname(p$allRRfit); out$allRRlow <- unname(p$allRRlow)
    out$allRRhigh <- unname(p$allRRhigh); out$matRRfit <- unname(p$matRRfit)
  } else { out$alllow <- unname(p$alllow); out$allhigh <- unname(p$allhigh) }
  out
}
redmeta <- function(r) {
  list(type = r$type, value = r$value, coef = unname(r$coefficients),
    vcov = unname(r$vcov), basis = unname(unclass(r$basis)), predvar = r$predvar,
    cen = r$cen, lag = r$lag, bylag = r$bylag, fit = unname(r$fit), se = unname(r$se),
    RRfit = unname(r$RRfit), RRlow = unname(r$RRlow), RRhigh = unname(r$RRhigh),
    low = unname(r$low), high = unname(r$high))
}
modelmeta <- function(m, cbname) {
  cf <- coef(m); vc <- vcov(m)
  ind <- grep(paste0(cbname, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}"), names(cf))
  if (length(ind) == 0) ind <- which(names(cf) == cbname)  # single-column basis
  list(coef_all = unname(cf), coef_names = names(cf), coef = unname(cf[ind]),
    vcov = unname(vc[ind, ind, drop = FALSE]), dispersion = summary(m)$dispersion,
    link = m$family$link, family = m$family$family)
}

cases <- list()

# Example 1: simple DLM ------------------------------------------------------
cb1.pm <- crossbasis(chicagoNMMAPS$pm10, lag = 15, argvar = list(fun = "lin"),
  arglag = list(fun = "poly", degree = 4))
cb1.temp <- crossbasis(chicagoNMMAPS$temp, lag = 3, argvar = list(df = 5),
  arglag = list(fun = "strata", breaks = 1))
model1 <- glm(death ~ cb1.pm + cb1.temp + ns(time, 7 * 14) + dow,
  family = quasipoisson(), chicagoNMMAPS)
pred1.pm <- crosspred(cb1.pm, model1, at = 0:20, bylag = 0.2, cumul = TRUE)
pred1.temp <- crosspred(cb1.temp, model1, by = 1, cen = 20)
wcsv(cb1.pm, "cb1_pm"); wcsv(cb1.temp, "cb1_temp")
cases$ex1 <- list(
  cb1_pm = cbmeta(cb1.pm), cb1_temp = cbmeta(cb1.temp),
  model_pm = modelmeta(model1, "cb1.pm"), model_temp = modelmeta(model1, "cb1.temp"),
  pred1_pm = predmeta(pred1.pm), pred1_temp = predmeta(pred1.temp))

# Example 2: seasonal analysis with group ------------------------------------
chicagoNMMAPSseas <- subset(chicagoNMMAPS, month %in% 6:9)
cb2.o3 <- crossbasis(chicagoNMMAPSseas$o3, lag = 5,
  argvar = list(fun = "thr", thr = 40.3), arglag = list(fun = "integer"),
  group = chicagoNMMAPSseas$year)
cb2.temp <- crossbasis(chicagoNMMAPSseas$temp, lag = 10,
  argvar = list(fun = "thr", thr = c(15, 25)), arglag = list(fun = "strata", breaks = c(2, 6)),
  group = chicagoNMMAPSseas$year)
model2 <- glm(death ~ cb2.o3 + cb2.temp + ns(doy, 4) + ns(time, 3) + dow,
  family = quasipoisson(), chicagoNMMAPSseas)
pred2.o3 <- crosspred(cb2.o3, model2, at = c(0:65, 40.3, 50.3))
pred2.temp <- crosspred(cb2.temp, model2, by = 1)
wcsv(cb2.o3, "cb2_o3"); wcsv(cb2.temp, "cb2_temp")
cases$ex2 <- list(
  cb2_o3 = cbmeta(cb2.o3), cb2_temp = cbmeta(cb2.temp),
  model_o3 = modelmeta(model2, "cb2.o3"), model_temp = modelmeta(model2, "cb2.temp"),
  pred2_o3 = predmeta(pred2.o3), pred2_temp = predmeta(pred2.temp))

# Example 3: DLNM with bs / ns and logknots ----------------------------------
cb3.pm <- crossbasis(chicagoNMMAPS$pm10, lag = 1, argvar = list(fun = "lin"),
  arglag = list(fun = "strata"))
varknots <- equalknots(chicagoNMMAPS$temp, fun = "bs", df = 5, degree = 2)
lagknots <- logknots(30, 3)
cb3.temp <- crossbasis(chicagoNMMAPS$temp, lag = 30, argvar = list(fun = "bs",
  knots = varknots), arglag = list(knots = lagknots))
model3 <- glm(death ~ cb3.pm + cb3.temp + ns(time, 7 * 14) + dow,
  family = quasipoisson(), chicagoNMMAPS)
pred3.temp <- crosspred(cb3.temp, model3, cen = 21, by = 1)
pred3.temp.sub <- crosspred(cb3.temp, model3, cen = 21, at = c(-20, 0, 27, 33), lag = c(0, 10), bylag = 0.5)
pred3.pm <- crosspred(cb3.pm, model3, at = 0:30)
wcsv(cb3.pm, "cb3_pm"); wcsv(cb3.temp, "cb3_temp")
cases$ex3 <- list(
  cb3_pm = cbmeta(cb3.pm), cb3_temp = cbmeta(cb3.temp), varknots = varknots, lagknots = lagknots,
  model_pm = modelmeta(model3, "cb3.pm"), model_temp = modelmeta(model3, "cb3.temp"),
  pred3_temp = predmeta(pred3.temp), pred3_temp_sub = predmeta(pred3.temp.sub),
  pred3_pm = predmeta(pred3.pm))

# Example 4: reduction -------------------------------------------------------
cb4 <- crossbasis(chicagoNMMAPS$temp, lag = 30,
  argvar = list(fun = "thr", thr = c(10, 25)), arglag = list(knots = lagknots))
model4 <- glm(death ~ cb4 + ns(time, 7 * 14) + dow, family = quasipoisson(), chicagoNMMAPS)
pred4 <- crosspred(cb4, model4, by = 1)
redall <- crossreduce(cb4, model4)
redlag <- crossreduce(cb4, model4, type = "lag", value = 5)
redvar <- crossreduce(cb4, model4, type = "var", value = 33)
redvar_by <- crossreduce(cb4, model4, type = "var", value = 33, bylag = 0.25)
wcsv(cb4, "cb4")
cases$ex4 <- list(cb4 = cbmeta(cb4), model = modelmeta(model4, "cb4"),
  pred4 = predmeta(pred4), redall = redmeta(redall), redlag = redmeta(redlag),
  redvar = redmeta(redvar), redvar_by = redmeta(redvar_by))

# Example 5: ns/ns with centring, reduction with cen, gaussian link ----------
cb5 <- crossbasis(chicagoNMMAPS$temp, lag = 21, argvar = list(fun = "ns", df = 4),
  arglag = list(fun = "ns", df = 4))
model5 <- glm(death ~ cb5 + ns(time, 7 * 14) + dow, family = quasipoisson(), chicagoNMMAPS)
pred5 <- crosspred(cb5, model5, cen = 21, by = 1, cumul = TRUE)
red5all <- crossreduce(cb5, model5, cen = 21, by = 1)
red5var <- crossreduce(cb5, model5, type = "var", value = 33, cen = 21)
red5lag <- crossreduce(cb5, model5, type = "lag", value = 3, cen = 21, by = 1)
model5g <- glm(death ~ cb5 + ns(time, 7 * 14) + dow, family = gaussian(), chicagoNMMAPS)
pred5g <- crosspred(cb5, model5g, cen = 21, by = 1)
wcsv(cb5, "cb5")
cases$ex5 <- list(cb5 = cbmeta(cb5), model = modelmeta(model5, "cb5"),
  model_gauss = modelmeta(model5g, "cb5"),
  pred5 = predmeta(pred5), red5all = redmeta(red5all), red5var = redmeta(red5var),
  red5lag = redmeta(red5lag), pred5g = predmeta(pred5g))

# Example 6: lag-matrix input (exposure histories) and ps basis --------------
Q <- exphist(chicagoNMMAPS$temp[1:200], lag = c(2, 8))
cb6 <- crossbasis(Q, lag = c(2, 8), argvar = list(fun = "ns", df = 3), arglag = list(fun = "poly", degree = 2))
wcsv(Q, "exphist_Q"); wcsv(cb6, "cb6")
Q2 <- exphist(c(1, 2, 3, 4, 5), times = c(2, 5, 7), lag = 3)
Q3 <- exphist(c(1, 2, 3, 4, 5), lag = c(1, 3), fill = -1)
cb7 <- crossbasis(chicagoNMMAPS$temp, lag = 10, argvar = list(fun = "ps", df = 8), arglag = list(fun = "ps", df = 5))
pen7 <- cbPen(cb7)
wcsv(cb7, "cb7")
cb9 <- crossbasis(chicagoNMMAPS$temp, lag = 21, argvar = list(fun = "cr", df = 6), arglag = list(fun = "cr", df = 5))
pen9 <- cbPen(cb9)
wcsv(cb9, "cb9")
cases$ex6 <- list(cb6 = cbmeta(cb6), Q2 = unclass(Q2), Q2_rows = rownames(Q2),
  Q3 = unclass(Q3), cb7 = cbmeta(cb7),
  pen7 = list(Svar = unname(pen7$Svar), Slag = unname(pen7$Slag), rank = unname(pen7$rank)),
  cb9 = cbmeta(cb9), pen9 = list(Svar = unname(pen9$Svar), Slag = unname(pen9$Slag), rank = unname(pen9$rank)))

# Example 7: single-lag (lag=0) crossbasis and onebasis prediction -----------
ob <- onebasis(chicagoNMMAPS$temp, "ns", df = 4)
modelob <- glm(death ~ ob + ns(time, 7 * 14) + dow, family = quasipoisson(), chicagoNMMAPS)
predob <- crosspred(ob, modelob, cen = 21, by = 1)
cb8 <- crossbasis(chicagoNMMAPS$temp, lag = 0, argvar = list(fun = "bs", df = 5))
cases$ex7 <- list(
  ob_attr = attrs(ob, keep), ob_coef = unname(coef(modelob)[grep("^ob", names(coef(modelob)))]),
  ob_vcov = unname(vcov(modelob)[grep("^ob", names(coef(modelob))), grep("^ob", names(coef(modelob)))]),
  predob = predmeta(predob), cb8 = cbmeta(cb8))
wcsv(cb8, "cb8"); wcsv(ob, "ob")

wjson(cases, "cases")
cat("fixtures written to", outdir, "\n")
