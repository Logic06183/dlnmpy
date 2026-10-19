# Number audit of paper/manuscript.qmd

Each number the manuscript states, where it comes from, and whether it matched. "Inline" means the manuscript now computes it at render time, so it cannot drift from its source. "Static, checked" means it stays typed in the text and was checked by hand against the named file. Runs are from a fresh clone of `v0.8.4` (then of `submission-prep`) in a Python 3.11 environment built from `paper/requirements.txt`.

| Where | Number in text | Source | Value in source | Status |
|---|---|---|---|---|
| Abstract, 5.2 | 98 comparisons with R, tolerance 1e-8 | `tools/side_by_side.py` output: 102 rows, 4 of them coverage checks; tolerances in the script | 98, all `1e-10 + 1e-8|r|` or tighter | matched. Abstract now says "comparisons" (each row is a vector or matrix, not one quantity). 5.2 now computes 98 inline |
| Abstract, 5.3 | MMP 89.5, AF 8.94% | `tests/reference/lancet2015.json`; `examples/lancet2015_vs_r.py` output | see log | see log |
| Abstract, 5.3 | MMT 17.1 °C, RR 1.101 at 22 °C, 1.309 at 0 °C | `tests/reference/bmcmrm2013.json` `published/second_stage/{mmt,rr22,rr0}` | 17.1; 1.10119 (1.07872 to 1.12413); 1.30873 (1.24489 to 1.37585) | matched |
| 5.3 | I² 63.8, 16.4, 63.5 | same, `pooled/{all,hot,cold}/qtest/I2` | 63.803, 16.416, 63.499 | matched |
| 5.3 | QAIC 402827.3, 406212.0, 405966.3 | same, `published/first_stage/qaic_total` | 402827.284, 406212.041, 405966.268 | matched |
| 5.3 | log-likelihood 168.244063 (R default) and 168.244068 (optimum) | same, `converged/second_stage/pooled/hot/logLik`, `reltol14/pooled/hot/logLik` | 168.2440626, 168.2440678 | matched |
| 5.3 | Printed 2013 paper values 1.101 (1.078 to 1.124), 1.308 (1.245 to 1.375), I² 63.7/16.0/63.4 | Gasparrini & Armstrong 2013 as printed (your brief) | as stated | matched to your brief; I could not open the paper from here |
| 5.3 | Vicedo RCP8.5 heat 2090 to 2099: 3.95% (0.69 to 6.75); paper 4.0 (0.7 to 6.8) | `tests/reference/vicedo2019.json` `afrel[2090-99, , heat, RCP8.5]` | 3.954 (0.688 to 6.753) | matched; now inline |
| 5.1 | 36 basis specifications | `tests/fixtures/onebasis.json` | 36 | matched; now inline |
| 5.1 | **9 cross-bases** | `tests/fixtures/cases.json` keys `cb*` | **12** | **wrong; fixed** (inline) |
| 5.1 | 12 prediction objects | `cases.json` `pred*` (11) plus `penalized.json` (1) | 12 | matched; now inline |
| 5.1 | 7 reductions | `cases.json` `red*` | 7 | matched; now inline |
| 5.1 | 142 tests | `pytest --collect-only` | see log | now inline |
| 5.5 | 41 defects, 28 silent, 13 not, 16 port, 25 layer, 28 audits, 11 reproductions, 1 user, 1 unknown, 5 without a test | `docs/defects.md` ledger rows | all as stated | matched; now inline |
| Limitations | R runs the Stat Med analysis in "about a minute" | `tests/reference/statmed2014.json` `versions.runtime_seconds` | 30.9 s, for the R script that also refits each model to convergence | **not supported; fixed** (inline from the JSON) |
| Limitations | `fit_clogit` 30 to 90 s per model; `cr` basis about 70 s for 3 million values; Chicago fit 4 s; attribution table 15 s; Stat Med translation 8 min; render 20 min | no committed output records these | n/a | **untraceable; replaced** by times measured during the render (inline) and a profile in the log |
| 5.3, Lancet stage table | coef 1.5e-5 (3.0e-6), Psi 3.5e-6 (3.5e-7), AN 0.22 (1.5e-5), AF 3.4e-5, reduced coef 1.0e-13 | `examples/lancet2015_vs_r.py` group table, this machine | 1.1e-5 (1.1e-6), 3.0e-6 (1.2e-7), 0.21 (1.4e-5), 3.1e-5, 6.0e-14 | **platform-dependent at the flat meta-regression; table now built from the render's own run**. The Mac values in `docs/validation.md` kept, with the Linux run added |
| 5.3 | mvmeta and mixmeta differ by 5.0e-6 | `lancet2015.json`, printed by the script | 5.0e-6 | matched; now inline |
| 5.3 | mvmeta and mixmeta agree to 6e-8 at the optimum | not in `lancet2015.json` (no mvmeta fit at the tight tolerance is stored) | n/a | **[CHECK] in the manuscript** |
| 5.3 | dlnmpy log-likelihood 1.3e-8 higher than R | script output: 52.3899754850 vs 52.3899754715 | 1.35e-8 higher | matched; now inline, sign included |
| Abstract | 161,356 comparisons | sum over the six scripts: 2,545 + 22,771 + 102,224 + 22,690 + 11,060 + 66 | see final render | Lancet 2,545, Leone 22,690, Vicedo 11,060 confirmed so far |
| 5.3 | Stat Med: coefficients and covariance to 1e-9, predictions within 1e-6 of R's SE, peak lags identical; 102,224 comparisons | `statmed2014.py` output | Table II 1.1e-9 (3.6e-11 SE), predictions 4.4e-7 SE, peak lag exact; 102,224, 0 failures | matched |
| 5.3 | Leone: fractions to 1e-13, intervals from R's draws to 4e-6 relative; 22,690 | `leone2014.py` output | converged AF 9.3e-14 relative; eCI R draws 3.5e-6 relative; 22,690 | matched |
| 5.3 | Vicedo: fractions to 1e-12, intervals to 1e-10; 11,060 | `vicedo2019.py` output | AF 8.8e-13; AF eCI 9.2e-11; 11,060 | matched |
| 5.3 | London coefficients 1.5e-8 short (R default glm) | `leone2014.py` output, cb coefficients rel | 1.55e-8 | matched; now inline |
| 5.3 | model 8 of Stat Med 1.6e-6 SE short | `statmed2014.json`, published vs converged / SE | 1.64e-6 | matched; now inline |
| 5.3 | "a reduction summed over 39 lags magnifies [this] to 2e-4 standard errors" | recomputed: `crossreduce` of the stored published and converged model-8 fits | 1.9e-6 SE (reduced coefficients and overall curve) | **wrong; removed** |
| 5.3 | `glm` and `coxph` stop at "about 1e-8 in the deviance" | R defaults (`glm.control` epsilon 1e-8; `coxph.control` eps 1e-9, as the script prints) | 1e-8 and 1e-9 | **imprecise; fixed** |
