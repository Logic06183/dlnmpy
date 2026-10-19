# Submission-prep log

Working machine: Claude Code cloud container (Ubuntu 24.04, 4 cores), not the Mac. Branch `submission-prep` from `8c5edd7` (main).

## 2026-10-19

### Setup
- Installed Quarto 1.6.43 (GitHub release), Python 3.11.15 (uv), TeX Live from Ubuntu apt. No R in the container: CRAN is blocked by the network policy. Every R comparison therefore runs against the stored results in `tests/reference/` and `tests/side_by_side/`, as a journal rerun would.

### Task 1: journal guidelines
- The journal's domains do not resolve from here and the proxy refuses them. Wrote `.handoff/journal-checklist.md` from search-engine snippets, with conservative defaults where the rule is unknown. Flagged: ISSN pending, not yet in DOAJ, founding-article waiver not visible in any snippet.

### Task 5 (started early, because it is slow): cold clone of v0.8.4
- `git clone` from GitHub into a scratch directory, `uv venv -p 3.11`, `pip install -r paper/requirements.txt`.
- Finding: `paper/requirements.txt` pins `dlnmpy==0.8.4` **from PyPI**, so the paper runs the wheel, not the clone. For the release tag those should be the same code, but a reviewer changing the source would not see the change. Fixed below.
- Finding: the manuscript says `paper/.python-version` exists. It does not.

### Tasks 2 to 4 (in progress)
- Abstract cut from 383 to 249 words. Same claims. "98 quantities" is now "98 comparisons", because each row of `tools/side_by_side.py` compares a vector or matrix. The script prints 102 rows, of which 4 are coverage checks against the simulated truth.
- Introduction said "the two published reproductions in Section 5 are rerun". There are six. Fixed.
- Section 5.1 said "9 cross-bases". `tests/fixtures/cases.json` holds 12. Fixed, and the fixture, prediction, reduction and test counts are now computed from the repository at render time.
- Defect counts (41, 28/13, 16/25, 28/11/1/1, 5 untested) match `docs/defects.md` row by row. Now computed from the ledger at render time, with an assertion that the categories add up.
- Vicedo 2019, RCP8.5 heat: 3.95 (0.69 to 6.75) matches both `tests/reference/vicedo2019.json` (3.954, 0.688 to 6.753) and the script's printed output. Now parsed from the output at render time. The text said "heat-attributable increase"; it is the change in the attributable fraction from 2010 to 2019, in percentage points, and now says so.
- Gasparrini & Armstrong 2013 values (MMT 17.1, RRs, I², QAIC totals, the two log-likelihoods) all match `tests/reference/bmcmrm2013.json`.
- "8 of 20 eigenvectors with opposite signs" (Leone 2014): reproduced with `numpy.linalg.eigh` against R's stored eigenvectors (8 of 20; `scipy.linalg.eigh` gives 4, which is the point of the fix).
- **Stale numbers printed in the paper.** `examples/lancet_2015.py` printed hard-coded "[published code]" intervals 8.42-9.43, 8.13-9.12 and 0.27-0.34. `docs/validation.md` already says these came from an R run with no recorded seed and cannot be reproduced; the stored R result is 8.40-9.48, 8.08-9.19, 0.27-0.35. The manuscript prints this output and calls the bracketed figures "that script's results". Both headline scripts (`lancet_2015.py`, `bmcmrm_2013.py`) now read their bracketed R values from `tests/reference/`. No package code changed.
- Limitations: the timings (clogit 30 to 90 s per model, `cr` 70 s, Chicago 4 s and 15 s, Stat Med 8 min against "about a minute" in R, render 20 min) are not recorded by any committed script. The R figure is contradicted by `tests/reference/statmed2014.json`, which records 31 s for the R script including the refits. Replaced with times measured during the render.
- README and CHANGELOG reported who suggested the mixmetapy split (private correspondence). Removed.
- `docs/validation.md` said the old QAIC "reversed the paper's conclusion". Reworded: the paper does not use QAIC; the authors' script ranks the models and the old criterion ranked them differently. Same correction annotated on the 0.5.0 CHANGELOG entry.
- Manuscript said "There is no hosted continuous integration". The repository has GitHub Actions workflows, but every run since at least 21 September fails in about 2 seconds with no steps executed (account-level, most likely billing or Actions minutes). Including the release workflow, so PyPI uploads were manual. Text corrected; listed for you.
- Added software references (R Core Team, mvmeta, mixmeta, mgcv and Wood 2017, survival and Therneau and Grambsch 2000, MASS, tsModel, Quarto, geopandas, pandas Zenodo) and a computational-details section whose version table is read from the rendering environment and from the stored R results.
- DOIs: `doi.org` is unreachable from this container, so no DOI could be resolved live. Wrote `tools/check_dois.py`; every DOI was checked by hand for form against my knowledge of the papers. See REPORT.
- DOIs confirmed through PubMed (registered, and matched to the right paper): the 14 journal articles in biomedical journals (Gasparrini 2010, 2012, 2013, 2014 x2, 2015, 2017, 2021; Sera 2019; Tobías 2017; Vicedo-Cabrera 2019; Peng 2011; Harris 2020; Virtanen 2020). Not checkable from here: JSS (dlnm 2011), JRSS-B (Wood 2011), the two SciPy-proceedings DOIs, Matplotlib (CiSE), xarray (JORS), the CRAN package DOIs, the Zenodo DOIs and the three books. `python tools/check_dois.py` on any machine with normal internet checks all of them in a few seconds.
- `examples/biomet2017.py` docstring said dlnmpy does not accept `int=`; it has since 0.8.4. Fixed.

### Statistical review (added goal: until a top statistical reviewer would approve)
- Round 1, three independent reviewers (core DLNM; attribution, MMT and QAIC; model fitting and penalised REML), every claim backed by code they ran. Verified correct by brute force or simulation: cross-basis, crosspred/crossreduce SEs and covariance, centring, ci_level, spline extrapolation, penalties, quantile/pretty ports, attrdl backward/forward, MMT search, eigen square root, QAIC (Gasparrini's fqaic), exact conditional likelihood and its covariance, REML/ML scores, edf, Vp, Fletcher scale.
- Found and fixed on `submission-prep`, each with a regression test that fails without the fix, a CHANGELOG entry (0.8.5, unreleased) and a ledger row: D42 prediction lookup by np.isclose (wrong row on fine grids of large values); D43 LogC link read as log; D44 attr_table accepted non-log links; D45 fit_glm zero-weight rows in the residual df; D46 fit_glm aliasing tolerance (lm.fit's, not glm.fit's); D47 fit_clogit silent non-convergence on separated data; D49 fit_pgam zero-weight rows in the scale. Safeguards: model_grid warns on unequal n; penalised search warns when it fails.
- D48, attr_table double-counting days on the 2.5th/97.5th percentile in both extreme and moderate components (moderate heat about 60% too high for Chicago), changes results already reported, so it is on branch `fix/attr-table-partition` for your decision. The manuscript states the current convention honestly.
- Manuscript corrections from the review: the forward perspective's assumption (checked against Gasparrini and Leone 2014 full text through PubMed); 7 df/year is the vignette's, the Lancet used 8; Chicago interval limits rounded and their Monte Carlo SD computed over five seeds; clogit algorithm described as it is; log-determinant limitation.
- Round 2, a fresh referee: "minor revision", with a list of what would make it an accept. Done: D49; offline test tier on the stored R results (`tests/test_reference_offline.py`, 7 tests, 27 s) and a slow tier running the six reproductions under pytest; abstract states the tested scope; function-by-layer coverage table computed from the source; counting rule of the 161,356 comparisons stated (five scripts count elements, biomet2017 counts 66 quantities over 142,275 elements); the untraceable "6e-8" claim removed; independence caveat extended to the R driver scripts; inference caveats paragraph and docs/theory.md section; side-by-side "coverage" rows renamed. Remaining minors documented (qaic ignores weights, as fqaic; warm-start noise in the penalised search; MMT grid).
- Also found while auditing: the "2e-4 SE" magnification claim for the Stat Med reduction does not reproduce (1.9e-6 SE); removed. The Lancet stage table is platform-dependent at the flat meta-regression; now computed from the render. README showed withdrawn Lancet intervals as R's; fixed. docs/penalized.md overstated agreement with mgcv relative to what its tests check; fixed.

### Reproducibility
- v0.8.4 cold baseline (fresh clone, PyPI wheel): all six reproductions pass, 2,545 + 22,771 + 102,224 + 22,690 + 11,060 + 66 = 161,356 comparisons, 0 failures. Times on 4 shared cores under heavy contention: Lancet 23 min, BMC 2013 50 min, Stat Med 28 min, Leone 12 min, Vicedo 12 min, biomet 60 min. The test run was cut off by a session time limit.
- `paper/requirements.txt` now pins all 123 packages of the rendering environment; dlnmpy installs editable from the clone and the manuscript refuses to run any other copy.
