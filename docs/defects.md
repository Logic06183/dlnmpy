# Defect ledger, 0.1.0 to 0.8.4

This ledger is built from `CHANGELOG.md` (0.1.0 to 0.8.4), with each entry checked against `git log -p` (tags v0.7.0 to v0.8.3, and the untagged commits before them) and the `tests/` directory. It is current as of 13 October 2026 and was updated on 19 October 2026.

## Counting rule

A defect is a behaviour that gave a wrong number, crashed, hung, or rejected valid input that R accepts, and which was then changed. Pure additions (new functions, exports, features), documentation-only changes, refactors with no behaviour change, and changes the CHANGELOG itself describes as fidelity with no observed effect are excluded. The excluded items are listed in their own section below, with the reason for each.

Silent? is yes when the result was wrong with no error or warning at the point of the call, and no when the call crashed, raised an error, hung or refused the input.

Port / layer is port for a function that exists in R `dlnm` or in the published `attrdl.R`, `findmin.R` or `mixmeta`. It is layer for fitting wrappers, the `dlnm()` workflow, plotting, uncertainty helpers, I/O and DataFrame interop. Plotting is classed as layer, even where it mirrors `plot.crosspred`.

Regression test names a test in `tests/` that targets the behaviour. "none found" means no such test exists in this repository.

## Ledger

| ID | Release fixed | Module / function | Description | Silent? | How found | Port / layer | Regression test | Notes |
|---|---|---|---|---|---|---|---|---|
| D01 | 0.4.1 | `uncertainty.qaic` (and `model_grid`) | Used statsmodels' `results.llf`, which is divided by the dispersion for quasi families, so QAIC was out by about 4700 and the `model_grid` ranking reversed | yes | other (not stated) | layer | `test_uncertainty.py::test_qaic_matches_r`, `::test_model_grid_ranking_matches_r` | Neither the CHANGELOG nor commit 400bf01 says how it was found. The commit says the previous test restated the implementation. |
| D02 | 0.4.2 | `core` basis `to_dataframe()` | Dropped the pandas index, so `data.join(cb.to_dataframe())` on filtered or grouped data was all-NaN with the row count unchanged | yes | reproduction of Gasparrini et al. 2015 *Lancet* | layer | `test_model.py::test_to_dataframe_keeps_the_input_index` | |
| D03 | 0.4.2 | `model.extract_coef_vcov` / `crosspred` (coefficient name matching) | Basis `cb` also matched the columns of a second basis `cb2`, so `crosspred` raised "found 4 coefficients matching basis 'cb'" | no | reproduction of Gasparrini et al. 2015 *Lancet* | port | `test_model.py::test_basis_name_is_not_a_prefix_match` | Classed as port because it is the coefficient lookup inside `crosspred(model=)`. |
| D04 | 0.5.0 | `model.fit_glm` | `offset`, `exposure`, `freq_weights`, `var_weights` were passed to `GLM.fit()` and silently discarded | yes | audit against R | layer | `test_model.py::test_fit_glm_offset_and_weights_reach_the_model` | |
| D05 | 0.5.0 | `crossbasis` (`argvar` `cen`) | Centring value in `argvar` was discarded, so later `crosspred`/`crossreduce` centred on the automatic value (RR 1.343 vs R 1.166) | yes | audit against R | port | `test_predict.py::test_argvar_cen_is_kept_by_crossbasis` | |
| D06 | 0.5.0 | `basis.strata` | Tied quantile breaks gave an empty stratum and a rank-deficient basis; values outside the edges gave zero rows where R gives NA | yes | audit against R | port | `test_basis.py::test_strata_rejects_tied_breaks_and_is_na_outside_the_edges` | Two behaviours in one CHANGELOG bullet, counted once. |
| D07 | 0.5.0 | `plot(..., "slices", cumul=True)` | Lag looked up on `predlag` instead of integer lags, so with `bylag != 1` the wrong lag was drawn (lag 2 drew lag 10) | yes | audit against R | layer | `test_plot.py::test_cumulative_slices_use_integer_lags` | `lag >= 4` raised `IndexError` in the same case; the main case was silent. Mirrors `plot.crosspred`. |
| D08 | 0.5.0 | `model.get_link` (used by `crosspred`) | `cloglog` and `loglog` were reported as `log`, so `crosspred` exponentiated and reported RRs | yes | audit against R | port | `test_model.py::test_get_link_compound_names` | Port of dlnm's internal `getlink`. |
| D09 | 0.5.0 | `mkcen` | `cen=np.True_`/`np.False_` were coerced to 1.0/0.0 | yes | audit against R | port | `test_predict.py::test_numpy_boolean_cen` | |
| D10 | 0.5.0 | `basis.bs` | Refused data whose interior knots all coincide with a boundary knot, which `splines::bs` accepts with a warning | no | audit against R | port | `test_basis.py::test_bs_warns_where_ns_stops_on_tied_boundary_knots` | CHANGELOG lists under "Divergences from R". |
| D11 | 0.5.0 | `basis.thr` (default threshold) | Default threshold dropped NaN; R's `median(x)` without `na.rm` gives NA, so results differed from R on the same input | yes | audit against R | port | `test_basis.py::test_thr_default_threshold_propagates_na_as_in_r` | Borderline: the old number was arguably more useful, but it was not R's. CHANGELOG lists under "Divergences from R" with a behaviour-change flag. |
| D12 | 0.5.0 | `plot.summary_figure` | `var=[one value]` raised `AttributeError` | no | audit against R | layer | `test_plot.py::test_overlay_legend_title_and_single_value` | |
| D13 | 0.6.0 | `PenalizedGLMResults.edf_by` | Prefix match summed the edf of `cb` and `cb_o3` (38.1 vs 9.4) | yes | audit against R | layer | `test_fixes_0_6.py::test_penalty_prefix_does_not_capture_a_longer_basis_name` | Same root cause as D14; one CHANGELOG bullet, split because one is silent and one is not. |
| D14 | 0.6.0 | `penalized.fit_pgam` (penalty placement) | Same prefix match placed penalties on the wrong columns and failed with a shape error | no | audit against R | layer | `test_fixes_0_6.py::test_penalty_prefix_does_not_capture_a_longer_basis_name` | See D13. |
| D15 | 0.6.0 | `uncertainty.qaic` | Raised `AttributeError` on a penalised fit | no | audit against R | layer | `test_fixes_0_6.py::test_qaic_accepts_penalised_fits` | Borderline: R has no `qaic`, and the fix also defines new behaviour (edf as parameter count). CHANGELOG lists it as a fix. |
| D16 | 0.6.0 | `penalized.fit_pgam(offset=)` | Full-length offset was misaligned after rows with missing values were dropped | yes | audit against R | layer | `test_penalized_guards.py::test_d16_full_length_offset_follows_the_rows_dropped_for_missing_values` (added 0.8.5) | Code change confirmed in 7c3ef8e (`penalized.py`); no test touches `fit_pgam` with an offset and missing rows. |
| D17 | 0.6.0 | `meta.MixMeta.qtest` | p-value computed as `1 - cdf`, rounding to exactly 0 for large Q | yes | audit against R | port | none found | Test `test_qtest_pvalue_is_not_rounded_to_zero` was added in 7c3ef8e and removed in ccd3899 (0.8.0) when the meta-analysis moved to mixmetapy. |
| D18 | 0.6.1 | `workflow.dlnm(time=)` | A date column was read as nanoseconds, the seasonal spline asked for 8.5e12 df and the OS killed the process with no traceback | no | user run | layer | `test_workflow.py::test_datetime_time_column_matches_day_index` (also `::test_time_index_converts_datetimes_to_days`, `::test_object_dtype_dates_are_understood`, `::test_absurd_df_per_year_raises_rather_than_exhausting_memory`) | "Running the published 0.6.0 wheel through a full Chicago analysis in a notebook". Who ran it is not stated. |
| D19 | 0.7.0 | `model.fit_glm(exposure=)` | Aliased-column refit passed the stored `log(exposure)` back as `exposure`, logging it twice (intercept 2.17 vs R -9.98) | yes | audit against R | layer | `test_fixes_0_7.py::test_exposure_is_not_logged_twice_on_the_aliased_refit` | |
| D20 | 0.7.0 | `uncertainty.bootstrap_ci` | Centred the draws of an uncentred (`cen=False`) prediction, so the estimate fell outside its own interval | yes | audit against R | layer | `test_fixes_0_7.py::test_bootstrap_ci_keeps_an_uncentred_prediction_uncentred` | |
| D21 | 0.7.0 | `workflow.DLNM.rr_at`, `.summary` | Exponentiated for any link (Gaussian "RR" of 8e9 for an effect of 26) | yes | audit against R | layer | `test_fixes_0_7.py::test_dlnm_reports_an_additive_effect_for_a_gaussian_model` | |
| D22 | 0.7.0 | `attribution.attrdl`, `attr_table` | On a `group=` cross-basis, lagged exposures were formed across group boundaries | yes | audit against R | port | `test_fixes_0_7.py::test_attrdl_lags_within_groups` | Port because `attrdl` is; the `group` argument itself is a Python extension and `attr_table` is layer. |
| D23 | 0.7.0 | `attribution.attr_table` | With a reference beyond the 2.5th/97.5th percentile the cold ranges ran backwards (silently 0) or straddled the heat side | yes | audit against R | layer | `test_fixes_0_7.py::test_attr_table_components_never_straddle_the_reference` | |
| D24 | 0.7.0 | `uncertainty.qaic` | Returned a Poisson-likelihood value for Gaussian or quasi-binomial fits | yes | audit against R | layer | `test_fixes_0_7.py::test_qaic_refuses_non_poisson_families` | |
| D25 | 0.7.0 | `crossreduce(type="lag"/"var")` | A vector `value` silently used its first element; R stops | yes | audit against R | port | `test_fixes_0_7.py::test_crossreduce_refuses_a_vector_value` | |
| D26 | 0.7.0 | `basis.cr`, `cbpen` | `cr` with fewer than three knots gave a one-column basis and zero penalty, and `cbpen` produced an all-NaN penalty with no error | yes | audit against R | port | `test_fixes_0_7.py::test_cr_needs_three_knots_and_poly_a_positive_degree`, `::test_cbpen_refuses_a_zero_penalty` | |
| D27 | 0.7.0 | `basis.poly` | `degree=0` returned a basis with no columns | yes | audit against R | port | `test_fixes_0_7.py::test_cr_needs_three_knots_and_poly_a_positive_degree` | CHANGELOG does not say what R does or whether a downstream call would have failed; classed silent because the call itself raised nothing. |
| D28 | 0.7.0 | `crossbasis(group=)` | Missing group label left rows all-NaN and reported the wrong group count; wrong-length `group` not checked | yes | audit against R | port | `test_fixes_0_7.py::test_group_with_missing_labels_is_refused` | |
| D29 | 0.7.0 | `model.fit_clogit` | Failed with `MissingDataError` on the NaN rows every time-series cross-basis has; `clogit` drops them | no | audit against R | layer | `test_fixes_0_7.py::test_fit_clogit_drops_rows_with_missing_values` | |
| D30 | 0.7.0 | `penalized.fit_pglm(sp=)` | A scalar `sp` with one penalty failed on a 0-d array | no | audit against R | layer | `test_fixes_0_7.py::test_fit_pglm_accepts_a_scalar_sp` | |
| D31 | 0.7.0 | `meta.MixMeta.converged` | Never `False`: tested whether the objective was finite, not the optimiser status | yes | audit against R | port | none found | No test added in e66f671; meta tests later moved to mixmetapy. |
| D32 | 0.7.0 | `PenalizedGLMResults.converged` | Never `False`, same cause as D31 | yes | audit against R | layer | `test_penalized_guards.py::test_d32_converged_is_false_when_the_search_fails` (added 0.8.5) | Refined in f5e5717 (gradient check for Nelder-Mead). `test_penalized.py::test_smoothing_parameter_selection` only asserts `converged` is true. The fix was incomplete: see D36. |
| D33 | 0.8.1 | `model.fit_clogit` | All-NaN coefficients with no error on large-valued designs (statsmodels' Newton overflows); surfaced later as an unrelated `crosspred` error | yes | reproduction of Gasparrini 2014 *Stat Med* | layer | `test_fixes_0_7.py::test_clogit_survives_a_large_valued_design` | Test lives in the 0.7 fixes file. |
| D34 | 0.8.1 | `model.fit_clogit` (covariance) | Score differenced at a fixed 1e-4 step, so SEs and CIs were up to 5e-3 out when coefficients are small | yes | reproduction of Gasparrini 2014 *Stat Med* | layer | `test_fixes_0_7.py::test_clogit_covariance_step_follows_the_coefficient_scale` | CHANGELOG ties it to the same model as D33 by implication, not explicitly. |
| D35 | 0.8.1 | basis argument handling (`onebasis`, `crossbasis`, `ns`, `bs`) | `Bound=` (R partial match of `Boundary.knots`) raised `TypeError` | no | reproduction of Gasparrini & Leone 2014 *BMC MRM* and/or Vicedo-Cabrera et al. 2019 *Epidemiology* | port | `test_fixes_0_7.py::test_boundary_knots_accepts_r_partial_matching` | CHANGELOG lists it under "Changed" as a rough edge. Commit e091d7a names both papers without saying which one hit it. |
| D36 | 0.8.2 | `penalized.fit_pgam` / `fit_pglm` | Returned all-NaN coefficients, covariance and criterion with `converged=True` after sp ran to ~1e19 | yes | reproduction of Gasparrini et al. 2017 *Biometrics* | layer | `test_penalized.py::test_pgam_never_reports_convergence_on_a_nonfinite_fit` | Weak test: half of it asserts source text via `inspect.getsource`, and the behavioural half passes without checking anything if the fit is finite or raises `LinAlgError`. |
| D37 | 0.8.2 | `penalized` starting values (`_start_rho`) | Start taken from an unpenalised fit that diverges when only the penalty makes the model estimable; search aborted with a bare `LinAlgError` | no | reproduction of Gasparrini et al. 2017 *Biometrics* | layer | `test_penalized.py::test_pgam_starts_when_the_unpenalised_fit_diverges` | |
| D38 | 0.8.2 | `penalized` smoothing-parameter search | A trial sp making the penalised information singular aborted the whole search | no | reproduction of Gasparrini et al. 2017 *Biometrics* | layer | `test_penalized_guards.py::test_d38_a_singular_trial_point_does_not_abort_the_search` (added 0.8.5) | Only tolerated, not targeted, by the `except LinAlgError` branch of the D36 test. |
| D39 | 0.8.4 | `attribution.simulate_coef(normals=)` | Did not reproduce R's draws from R's normals, because eigenvector signs differ between LAPACK routines; `eigvecs=` added and signs now fixed deterministically | yes | reproduction of Gasparrini & Leone 2014 *BMC MRM* | layer | `test_fixes_0_8_4.py::test_simulate_coef_sign_does_not_depend_on_eigensolver`, `::test_simulate_coef_eigvecs_reproduces_mvrnorm_construction` | Classed as uncertainty helper; it is the `mvrnorm` step of `attrdl.R` and also feeds `attrdl`, `findmin`, `attr_table` and `bootstrap`. The tests check sign invariance and the mvrnorm construction, not R's actual draws; the R comparison is in `examples/leone2014.py`. Both tests fail on 0.8.3 and pass on 0.8.4. |
| D40 | 0.8.4 | `model.fit_clogit` (step-halving path) | Lost the parameter names, so `crosspred(name=)` failed | no | reproduction of Gasparrini 2014 *Stat Med* | layer | `test_fixes_0_8_4.py::test_fit_clogit_keeps_names_on_step_halving_path` | Regression introduced by the D33 fix in 0.8.1 (`res.params = b` in 3c6c74f). Fails on 0.8.3, passes on 0.8.4. |
| D41 | 0.8.4 | `core.normalise_args` (basis arguments) | `int=` (R partial match of `intercept`) was rejected | no | reproduction of Gasparrini 2014 *Stat Med* and Gasparrini et al. 2017 *Biometrics* | port | `test_fixes_0_8_4.py::test_int_is_partial_match_for_intercept` | Fails on 0.8.3, passes on 0.8.4. |
| D42 | 0.8.5 (unreleased) | `CrossPred._var_index`, `_lag_index` (`slice_var`, plots, `DLNM.rr_at`) | Matched prediction values with `np.isclose` at rtol 1e-5, so on a fine grid of large values the neighbouring row was used | yes | statistical review for the manuscript | layer | `test_fixes_0_8_5.py::test_prediction_lookup_picks_the_exact_row_on_a_fine_grid_of_large_values` | R looks predictions up by name, exactly. |
| D43 | 0.8.5 (unreleased) | `model.get_link` (used by `crosspred`) | statsmodels' `LogC` link was reported as `log`, so `crosspred` exponentiated it into a "relative risk" | yes | statistical review for the manuscript | port | `test_fixes_0_8_5.py::test_logc_link_is_not_reported_as_a_log_link` | Same function as D08. |
| D44 | 0.8.5 (unreleased) | `attribution.attr_table`, `DLNM.attributable` | Accepted a model without a log or logit link and returned meaningless fractions (AF -23,826 for a Gaussian model); `attrdl` refuses | yes | statistical review for the manuscript | layer | `test_fixes_0_8_5.py::test_attr_table_refuses_a_model_without_a_log_or_logit_link` | |
| D45 | 0.8.5 (unreleased) | `model.fit_glm` (dispersion) | Rows with zero prior weight counted in the residual df, so the quasi dispersion and SEs were too small (13% with a quarter of rows at weight zero) | yes | statistical review for the manuscript | layer | `test_fixes_0_8_5.py::test_fit_glm_dispersion_does_not_count_rows_with_zero_prior_weight` | R: `glm.fit`'s `n.ok`. |
| D46 | 0.8.5 (unreleased) | `model.fit_glm` (aliasing) | Aliasing tested at `lm.fit`'s 1e-7, where `glm.fit` uses 1e-11 on the weighted design, so near-collinear columns R keeps were dropped | yes | statistical review for the manuscript | layer | `test_fixes_0_8_5.py::test_fit_glm_aliasing_uses_glm_fit_tolerance` | No stored comparison with R reached the affected range. |
| D47 | 0.8.5 (unreleased) | `model.fit_clogit` (step-halving path) | On separated data returned a non-stationary point with no warning | yes | statistical review for the manuscript | layer | `test_fixes_0_8_5.py::test_fit_clogit_reports_non_convergence_on_separated_data` | `coxph` warns that a coefficient may be infinite. |

## Summary

The ledger lists 47 defects: 41 fixed in releases 0.4.1 to 0.8.4, and 6 fixed for the unreleased 0.8.5. The counts after the table are for the 41 released fixes. The six 0.8.5 defects (D42 to D47) were all silent, five in the layers and one in a ported function, were all found by a statistical review of the package for the manuscript, and all have a targeted regression test.

Their distribution by release is:

| Release | Defects |
|---|---|
| 0.1.0, 0.2.0, 0.3.0, 0.3.1, 0.3.2, 0.4.0 | 0 |
| 0.4.1 | 1 |
| 0.4.2 | 2 |
| 0.5.0 | 9 |
| 0.5.1 | 0 |
| 0.6.0 | 5 |
| 0.6.1 | 1 |
| 0.7.0 | 14 |
| 0.8.0 | 0 |
| 0.8.1 | 3 |
| 0.8.2 | 3 |
| 0.8.3 | 0 |
| 0.8.4 | 3 |
| 0.8.5 (unreleased) | 6 |

Of the 41, 28 were silent and 13 were not.

By location, 16 are in ported functions and 25 in the layers.

The defects were found as follows:

| Source | Defects |
|---|---|
| audit against R (0.5.0, 0.6.0, 0.7.0) | 28 |
| reproduction of a published paper | 11 (Lancet 2015: 2; Stat Med 2014: 4, one shared with Biometrics 2017; Biometrics 2017: 3; Leone 2014: 1; Leone 2014 or Vicedo-Cabrera 2019: 1) |
| user run | 1 |
| other / not stated | 1 |
| unit fixture | 0 |

Of the port defects, 13 came from the audits and 3 from reproductions (D03, D35, D41). Two of those three are R partial-argument matching.

Of the 41 defects, 39 now have a targeted regression test in `tests/`. The tests for D16, D32 and D38 were added for 0.8.5. D17 and D31 are in the meta-analysis code, which moved to mixmetapy in 0.8.0 (D17's test was removed with it). D36's test is weak.

Some calls are borderline, in case a stricter or looser count is wanted. D11, D15, D27 and the convergence-flag entries D31/D32 are included but arguable, and dropping all five gives 36. On the other side, the matplotlib 3.5 crash, the `dlnm(group=)` model change and the `__version__` drift are excluded but arguable (see below), and adding all three gives 44.

## Excluded items

| Release | Item | Why excluded |
|---|---|---|
| 0.1.0, 0.2.0, 0.3.0, 0.3.2, 0.4.0 | Initial port, `attrdl`/`findmin`/`mmt`/`attr_table`, mixmeta port, `simulate_cities`, `cr`, penalised fitting, themes, uncertainty module | Additions |
| 0.3.1 | `contour`/`3d` `ylab`/`zlab` follow R | Axis labels only, no wrong number |
| 0.3.1 | `ns`/`bs` warn when `df` is too small | Adds a warning; values unchanged |
| 0.4.0 | "R-consistent labels" | Labels only |
| 0.5.0 | `overlay_slices` legend titled with the x-axis label | Wrong label text, not a wrong number |
| 0.5.0 | `simulate_coef` exported | Addition |
| 0.5.1 | Spurious numpy 2.0.x `RuntimeWarning`s silenced | CHANGELOG says values never affected |
| 0.6.0 | `overlay_slices(labels=)` too few labels | Input was invalid before and after; only the error message changed |
| 0.6.0 | `bootstrap_ci` exported; workflow, figures, `to_dataframe`, `percentile_knots`/`percentile_of` | Additions |
| 0.6.0 | `mkat` passes `min.n` to `pretty()` | CHANGELOG: fidelity, no observed effect |
| 0.6.0 | `ps(df=4)` refused at construction | R also fails (later, in `crosspred`); not a defect in the port |
| 0.6.1 | Seasonal-spline df guard | Safeguard added alongside D18 |
| 0.6.1 | `__version__` hardcoded and drifted | Package metadata, not a computed result. Arguable. See CHANGELOG issues: the fix is not in the 0.6.1 commit |
| 0.6.1 | `basis.py` identical branches; `ruff` pinned | Refactor with no behaviour change; tooling |
| 0.7.0 | `mixmeta` with a missing outcome | Still refused; only the error message changed (see CHANGELOG issues) |
| 0.7.0 | `DLNM.figure()` missing the MMT band on a fresh fit | Display omission, no wrong number. Arguable |
| 0.7.0 | `DLNM.mmt()` caching | Same values before and after; caching and state |
| 0.7.0 | `CenteringWarning` | Addition |
| 0.7.0 | `dlnm(group=)` now fits per-group intercept and seasonality | CHANGELOG lists it as "Changed". The old pooled model was arguably misspecified. Arguable |
| 0.7.0 | `fit_glm` dispersion with `freq_weights` | Documentation only |
| 0.7.0 | CI, test skips, README, citation | Tooling and docs |
| 0.8.0 | Meta-analysis moved to mixmetapy | Refactor; CHANGELOG says no numerical change |
| 0.8.0 | Matplotlib minimum raised to 3.6 (`summary_figure` crashed on 3.5) | Addressed by a dependency floor, not a code change. Arguable |
| 0.8.0 | Version-independent tests; `.gitlab-ci.yml` | Tests and tooling |
| 0.8.0 (not in CHANGELOG) | ab2e977: sdist missing `conftest.py` and fixtures; clearer error when `mixmetapy` is shadowed | Packaging; error message only |
| 0.8.1 | `simulate_coef(normals=)` | Addition. Its claim to match R was wrong, which is D39 |
| 0.8.2 | README and `docs/penalized.md` no longer claim 1e-4 agreement with mgcv end to end | Documentation |
| 0.8.3 | Citation and Zenodo metadata | Metadata |
| 0.8.5 (unreleased) | `model_grid` warns when the candidate fits used different numbers of rows | A safeguard. The criteria it ranked were computed correctly, only not comparable |

## CHANGELOG issues found while checking

1. 0.8.0 is headed "(unreleased)", but tag `v0.8.0` exists on ab2e977 with `version = "0.8.0"`, and that commit's message refers to "the 0.8.0 sdist". The entry has no date.
2. 0.5.0 says "the pre-0.5.0 implementation would have reversed its published conclusion" (Gasparrini & Armstrong 2013). The QAIC fix shipped in 0.4.1 (400bf01), and commit 99f999c itself calls the reproduction "a direct test of the fix in 0.4.1". It should read pre-0.4.1.
3. The 0.6.1 `__version__` fix is not in 0.6.1. The 0.6.1 release commit (8534dcf) still hardcodes `__version__ = "0.6.0"` against `pyproject.toml` 0.6.1. The fix is in 01a475d, made afterwards, whose message says "the 0.6.1 wheel reported 0.6.0". It first shipped in 0.7.0. There is no `v0.6.1` tag to confirm what was uploaded.
4. 0.6.0 says "the defects below were all in the Python-only layers", but the `qtest` p-value (D17) is in the mixmeta port.
5. 0.7.0 says "everything below is in the Python layers", but six of its fourteen defects are in ported functions: `attrdl` (D22), `crossreduce` (D25), `cr`/`cbpen` (D26), `poly` (D27), `crossbasis(group=)` (D28) and `MixMeta` (D31). The claim that the numerical core needed no change holds only if "core" means the R-validated fixture paths.
6. 0.7.0 says `converged` "now reports the optimiser's own status". For penalised fits this was incomplete: the objective mapped NaN to 1e100, so an all-NaN fit could still report `converged=True`, which 0.8.2 fixed (D36).
7. 0.8.1 says "found by reproducing four more landmark papers" but names three (Stat Med 2014, Leone 2014, Vicedo-Cabrera 2019). The fourth, Biometrics 2017, is the subject of 0.8.2.
8. 0.8.1 calls `Bound=` a "rough edge" under "Changed". It raised `TypeError` on valid R input, which is a defect under this ledger's rule (D35).
9. 0.8.1 says that with `simulate_coef(normals=)` you can "feed in R's draws and the two languages agree without Monte Carlo noise". That was not true in general, because eigenvector signs differ between LAPACK builds (D39).
10. 0.8.1 reports Stat Med 2014 agreeing on 197 quantities, yet the step-halving path it added dropped the parameter names, so `crosspred(name=)` failed on exactly the model that needed it (D40). The 0.8.1 reproduction presumably worked around this. That cannot be checked from the repository.

For the 0.8.4 entry, note that fixing eigenvector signs in `_coefsim_from_normals` changes seeded Monte Carlo output from `attrdl(sim=True)`, `findmin` (and so `mmt`), `attr_table` and `uncertainty.bootstrap` relative to 0.8.3 for the same seed, wherever numpy returned an eigenvector whose largest entry was negative. The results are statistically equivalent but not bit-identical, and the CHANGELOG should say so.
