# Statistical review of dlnmpy and its manuscript

Reviewer: Claude (Anthropic), acting as statistical reviewer at Craig Parker's request. I am an AI, so this report is evidence for the co-authors' own judgement and not a substitute for it.

Branch `submission-prep`, dlnmpy 0.8.5.dev0, manuscript `paper/manuscript.qmd`. 20 October 2026.

## Recommendation

**Accept.** The package computes what it says it computes, its agreement with R is demonstrated on independent published code, its intervals behave as nominal intervals should where I tested them, and the manuscript states what the evidence supports and no more. Two items remain for the authors (below). Neither is statistical.

## What I examined, and what I found

**Method and port.** Four rounds of independent review, every claim backed by code that was run, covered the cross-basis and its tensor-product construction, prediction and reduction (standard errors from `M V M'`, centring, confidence levels), the spline bases and their extrapolation, the penalty matrices, R's `quantile`/`pretty`, attributable risk in both perspectives, the minimum mortality temperature and its simulated interval, QAIC (Gasparrini's `fqaic`), the GLM fitter (aliasing, dispersion, weights, offsets), conditional logistic regression (exact likelihood and its covariance), penalised REML/ML (scores against a brute-force restricted likelihood, edf, Vp, the Fletcher scale), the two-stage bridge to mixmetapy, exposure histories, the simulation helpers and the one-call workflow. These were confirmed correct by brute force, simulation or comparison with stored R results.

**Defects.** The review found nine silent defects (D42 to D50 in `docs/defects.md`). All are edge cases that none of the six published analyses reached, which is why the comparisons with R did not find them. The two that mattered most for applied users: a data frame with a repeated index made `dlnm()` multiply its rows (standard errors 1.24 times too small), and rows with zero prior weight shrank standard errors by about 14% in `fit_glm` and `fit_pgam`. Eight are fixed on this branch, each with a regression test that fails without the fix. The ninth (D48, `attr_table` counting days on the percentile cut-offs in two components) changes results 0.8.4 already reports, so it waits on the authors' decision on branch `fix/attr-table-partition`. The manuscript states the current behaviour.

**Agreement with R.** The six published analyses rerun from their authors' own R code agree in all 161,356 comparisons at the stated tolerances, rerun twice from fresh clones. The manuscript says how the comparisons are counted, which configurations they cover (daily time series and one cohort, log links, unit weights, penalised fits at mgcv's smoothing parameters) and which they do not.

**Calibration of the intervals (my own check).** I simulated 400 data sets of 3,000 days from a known temperature-lag-response surface and refitted each with `dlnm()`. The 95% intervals for the overall cumulative relative risk covered the true value 95.8%, 95.5% and 95.3% of the time at -10, 0 and 28 °C (Monte Carlo SE about 1 point): `.handoff/review_coverage.py`, output in `.handoff/review_coverage_output.txt`. An earlier referee's simulation found 97% coverage for the MMT interval and 93% for the attributable fraction, acceptable for simulation-based intervals that hold the reference and dispersion fixed, which the manuscript now says.

**The manuscript.** Every number is computed at render time or traced to a committed file (`.handoff/number-audit.md`), and the abstract is checked against the computed values. The inferential caveats are stated: the simulated intervals condition on the coefficients only, the MMT is a fixed reference, Ψ is treated as known for BLUPs, penalised fits use Vp, and Monte Carlo error is reported. The independence of the validation is described honestly: the authors' published analysis code is independent, and the R driver scripts and tolerances came from the same process as the port.

## What remains for the authors

1. Decide on D48 (`fix/attr-table-partition`). Statistically the fix is right: the extreme and moderate components should partition the days. Note that the original authors' scripts behave as the current code does.
2. The [CHECK] items (author details, grant number, manuscript tag, release 0.8.5 and its DOI). None affects the statistics.

## For each co-author, before approving

A short route to satisfy yourself, about an hour:

- Read Sections 5 and 8 of `paper/manuscript.pdf` (the checks and the limitations): are the claims ones you are willing to sign?
- Run `pytest` from a clone (four minutes) and, if you can spare the hour, `DLNMPY_SLOW=1 pytest -m slow`, which reruns the six published analyses against R.
- Fit one analysis you know from your own work with `dl.dlnm()` and compare it with what you got in R.
- Tell Craig your CRediT roles and anything you checked, so it can go in the Author contributions.
