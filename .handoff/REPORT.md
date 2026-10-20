# dlnmpy submission prep: report

Branch `submission-prep` (pushed). Details in `.handoff/log.md`, `.handoff/number-audit.md` and `.handoff/journal-checklist.md`.

## What changed

- **Manuscript.** The abstract went from 383 to 250 words. Every count, agreement level, runtime and the Lancet stage table are now computed at render time, and a final cell stops the render if the typed abstract disagrees with them. Code availability, computational details (with a version table read from the environment), author contributions, data availability and a function-by-check coverage table were added. The prose was cleaned of AI tells. Corrected errors: "9 cross-bases" (12); "two" reproductions rerun (six); the forward attribution assumption (checked against Gasparrini and Leone 2014); 7 df/year attributed to the Lancet analysis (it used 8); clogit's algorithm; untraceable or wrong runtimes; a "2e-4 SE" magnification claim (1.9e-6); stale "[published code]" intervals printed by `examples/lancet_2015.py`; and the counting rule behind 161,356.
- **Statistical review.** Four rounds of independent AI referees, each verifying claims by running code. Final verdicts: round 3 **accept**, round 4 **accept with minor edits** (those edits are done). They found nine silent defects, D42 to D50. The most consequential: `dlnm()` with a repeated data-frame index (stacked per-city frames) multiplied the rows, giving SEs 1.24 times too small; zero prior weights shrank SEs by about 14% in `fit_glm` and `fit_pgam`; `fit_glm` used `lm`'s aliasing tolerance; and `attr_table` double-counted days on the percentile cut-offs. Eight are fixed on `submission-prep`, each with a regression test that fails without it, a CHANGELOG entry and a ledger row; the ninth (D48) awaits you. Also added: an offline test tier against the stored R results, a slow tier for the six reproductions, workflow safeguards (time-step and repeated-index warnings), inference caveats and Monte Carlo error reporting. The referees were AI reviewers; a human statistician's review (a co-author's, or the journal's) is still the real test.
- **Reproducibility.** `REPRODUCE.md`; `paper/requirements.txt` pins all 123 packages; the paper runs dlnmpy from the clone and refuses any other copy; `paper/.python-version` added. Report of private correspondence removed from README and CHANGELOG.
- **Cold runs (fresh clones from GitHub, Linux x86_64, 4 CPUs, Python 3.11, nothing cached).** On `v0.8.4`: all six reproductions pass (161,356 comparisons, 0 failures). On `submission-prep`: 161 tests pass and 7 slow ones skip by design. `paper/manuscript.pdf` (23 pages) and `paper/manuscript.docx` were rendered from a fresh clone of commit `dd87878`, in 10 and 9 minutes; each render reran all six translations (0 failures) and the side-by-side analyses. I checked the PDF page by page: tables, equations, figure, code wrapping and references.

## [CHECK] items in the manuscript

1. Grant number U54TW012083 (Funding).
2. The Git tag of the submitted manuscript (Code availability).
3. Release 0.8.5 and its Zenodo version DOI, in two places (Reproducibility, Code availability). Until then the text says 0.8.5.dev0.

## Release decision (pending, yours)

0.8.5 is justified: eight silent defects are fixed. Nothing in the six published analyses, side-by-side or the paper changes. Branch `release-0.8.5` holds the version bumps only (not tagged, not uploaded, no Zenodo version). Fill in the date and conda sha256.

Separate decision: **`fix/attr-table-partition`** (D48). `attr_table` counted days exactly on the 2.5th/97.5th percentile in both extreme and moderate components. Chicago has 10 and 21 such days, so moderate heat is about 60% too high. The authors' own scripts do the same. The fix changes results 0.8.4 reports, so I stopped. To accept it: merge that branch into `submission-prep` (it is up to date with it), then the release commit. The manuscript on `submission-prep` states the current convention honestly.

## Only you can do

- Confirm grant number U54TW012083.
- Co-authors: Gia, Ben, Sam, Sibusiso Makhanya, Nicholas Brink (only if he approves the manuscript), Admire and Matthew are now in the author list, in that order, after you. Each has a visible [CHECK] for surname, affiliation and ORCID, and Author contributions has a [CHECK] for everyone's CRediT roles and any independent checks they made. If "Ben" is Ben Armstrong, the Acknowledgements (which thank him for dlnm and say none of the dlnm authors reviewed dlnmpy) and the statement that dlnmpy is not affiliated with the dlnm authors must change. Tell me and I will adjust them.
- Consider retitling to "... checked against R and against the original authors' published code": with several authors, "the authors'" in the title can be read as yours.
- Open the journal's author guidelines (unreachable from my machine) and tick `.handoff/journal-checklist.md`, in particular the abstract limit, required sections, reference style (add a `csl:` if they supply one) and the founding-article fee waiver in writing.
- Decide on 0.8.5 and the partition fix, then tag, release on GitHub and PyPI, mint the Zenodo version, and replace the three [CHECK]s.
- Run `python tools/check_dois.py` on a normal connection. The 14 biomedical DOIs were confirmed through PubMed. The rest (JSS, JRSS-B, SciPy proceedings, CiSE, JORS, CRAN, Zenodo, books) I could not resolve. Check the mixmetapy version DOI 10.5281/zenodo.22876733 in particular.
- GitHub Actions: every run fails in about 2 seconds with no steps executed, which points to a billing or Actions-minutes setting on the account. The manuscript says CI does not currently run.
- Merge `submission-prep` to main when you are happy.
- On your Mac: `rm -rf ~/Projects/dlnmpy/.release ~/Projects/dlnmpy/.git/stale-maintenance-lock-claude` and delete any `tmp_pack*` files (none of these are in the clone).
- Sign in to submit.fairpressjournals.com and upload the PDF with the repository link.

## Not done

Optional task 7 (speeding up `fit_clogit` and the `cr` basis) was not attempted. Every change to the numerical core reopens the full R validation, and the review fixes took priority.
