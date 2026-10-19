# Reproducing the paper

This file explains how to rerun everything in `paper/manuscript.qmd` from a fresh clone: the test suite, the six published analyses checked against R, and the manuscript itself. R is not needed. The R results the comparisons use are stored in `tests/reference/`, `tests/side_by_side/` and `tests/fixtures/`.

## What you need

- Python 3.11 (the version in `paper/.python-version`). Other versions from 3.9 work for the package, but the pins in `paper/requirements.txt` were resolved for 3.11.
- Quarto 1.6 or later, to render the manuscript.
- For the PDF, a LaTeX installation with XeLaTeX or pdfLaTeX. `quarto install tinytex` is enough.
- Network access to `raw.githubusercontent.com` on the first run. The data for the six published analyses are downloaded from pinned commits of their authors' repositories and checked against SHA-256 hashes. They are not redistributed here.
- About 1 GB of disk and, for the slowest analysis, about 2 GB of memory.

Nothing else is assumed. No credentials, no `~/.pypirc`, no local data and no R installation.

## Steps

```sh
git clone https://github.com/Logic06183/dlnmpy.git
cd dlnmpy
python3.11 -m venv .venv
. .venv/bin/activate
pip install -r paper/requirements.txt        # every Python package, pinned
pip install --no-deps -e .                    # dlnmpy from this clone

pytest                                        # the unit tests against R fixtures

cd examples                                   # the six published analyses
for f in lancet2015_vs_r bmcmrm2013_vs_r statmed2014 leone2014 vicedo2019 biomet2017; do
    python $f.py | tail -1                    # each exits non-zero on any failure
done
cd ..
python tools/side_by_side.py | tail -1        # the five side-by-side analyses

quarto render paper/manuscript.qmd            # HTML, PDF and DOCX; reruns all of the above
quarto render paper/manuscript.qmd --to pdf   # or one format at a time
```

Each format reruns every computation, so rendering all three takes three times as long as one.

On a machine with several cores, setting `OMP_NUM_THREADS=1` (and `OPENBLAS_NUM_THREADS=1`) avoids oversubscription when several scripts run at once. It does not change any result beyond the tolerances the scripts check.

## Where files are written

| What | Where | Notes |
|---|---|---|
| England and Wales, London, uranium miners data | `examples/*.csv` | downloaded on first use, ignored by git |
| Climate projection data | `examples/vicedo2019_data/` | downloaded on first use, ignored by git |
| Data for the penalised DLNM paper | `~/.cache/dlnmpy/biomet2017-531c5c8/` | set `DLNMPY_CACHE` to put it elsewhere |
| Rendered manuscript | `paper/manuscript.{html,pdf,docx}` | |

Delete any of these to force a fresh download. A file whose SHA-256 does not match is refused, never used.

## Times

Measured on the machine named in the computational details of the rendered manuscript. The slowest step is `examples/statmed2014.py` (Cox models on the uranium miners cohort), then `examples/biomet2017.py` (penalised DLNMs). The rendered manuscript prints the times measured during its own render.

## Regenerating the R side

Only needed to check the stored R results themselves. Each `tools/reproduce_*.R` pins the authors' repository, runs their code, and writes `tests/reference/<name>.json` with the R and package versions that produced it. The stored results were produced with R 4.3.3, dlnm 2.4.7, mvmeta 1.0.3, mixmeta 1.2.2, mgcv 1.9-1, survival 3.5-8, MASS 7.3-60.0.1 and tsModel 0.6-2. `tools/make_fixtures*.R` regenerate the unit fixtures and `tools/side_by_side.R` the side-by-side results.

```sh
Rscript tools/reproduce_lancet2015.R          # and likewise for the other five
```

## Checking the references

`python tools/check_dois.py` asks the DOI resolver for every DOI in `paper/references.bib` and exits non-zero if any fails.
