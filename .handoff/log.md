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
