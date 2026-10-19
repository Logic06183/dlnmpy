# Journal checklist: Journal of Reproducible Statistics (FAIR Press Journals)

Compiled 19 October 2026.

## Could the guidelines be read?

No. From this machine, `fairpressjournals.com`, `www.fairpressjournals.com`, `jrepstat.fairpressjournals.com` and `submit.fairpressjournals.com` did not resolve in DNS (WebFetch: `ENOTFOUND`), and the egress proxy refused them. The Wayback Machine was also unreachable. A web search does index two of the publisher's pages:

- <https://jrepstat.fairpressjournals.com/> ("Journal of Reproducible Statistics — FAIR Press Journals")
- <https://fairpressjournals.com/> ("FAIR Press Journals — reproducibility-first academic publishing")

From the search snippets only, which I could not open:

| Item | What the snippets say | Confirmed? |
|---|---|---|
| Publisher | FAIR Press Journals GmbH, Swiss | snippet only |
| Licence | Authors keep copyright, CC-BY | snippet only |
| Verification | Results are rerun and verified before publication | snippet only (matches your brief) |
| Fees | Low fee, need-blind waivers, submission free; a prize paid to authors of the strongest papers | snippet only |
| Indexing | ISSN pending; working towards DOAJ listing | snippet only |
| Founding-article waiver, 31 Oct 2026 | Not found in any snippet | **not confirmed** |
| Software article type | Not found | **not confirmed** |

**Before submitting, open the author guidelines yourself and tick the rows below.** The failure is almost certainly my network and not the journal: the same machine also cannot resolve `doi.org` or `zenodo.org`, while allow-listed hosts such as PyPI resolve. The pending ISSN and absent DOAJ listing are normal for a new journal, but it is worth having the founding-article fee waiver in writing before you submit.

## Requirements, with what the manuscript now does

Where the journal's rule is unknown, the target is a conservative default that meets the usual rules for software papers (JSS, JOSS, JORS), so it should need little change once you have the real numbers.

| # | Requirement | Journal rule | Target used | Manuscript state | Done |
|---|---|---|---|---|---|
| 1 | Article type | "software" (your brief) | `article-type: software` in YAML | set | yes |
| 2 | Abstract length | unknown | at most 250 words | see log; was 383 | see REPORT |
| 3 | Keywords | unknown | 4 to 6 | 6 | yes |
| 4 | Main-text length | unknown | no change made | about 7,000 words plus code | check |
| 5 | Required sections | unknown | Introduction, Methods/framework, Software design, Illustrations, Validation, Limitations, Reproducibility / computational details, Declarations, References | all present | check |
| 6 | Code availability | rerun before publication | repository URL, tag `v0.8.4`, Zenodo version DOI 10.5281/zenodo.23362540, concept DOI 10.5281/zenodo.22876735, licence | added | yes |
| 7 | Data availability | unknown | sources, pinned commits, SHA-256, not redistributed | present | yes |
| 8 | Computational details / versions | rerun before publication | Python, Quarto, every package version, R and R package versions used for reference results, hardware, runtimes | added | yes |
| 9 | Single entry point to reproduce | rerun before publication | `REPRODUCE.md` at repo root, one command | added | yes |
| 10 | Seeds | rerun before publication | every seed stated | present | yes |
| 11 | Funding statement | unknown | stated | grant number awaiting your confirmation | **you** |
| 12 | Competing interests | unknown | stated | "None" | yes |
| 13 | AI-use statement | unknown | stated | present | yes |
| 14 | Author contributions (CRediT) | unknown | stated for sole author; add co-authors | added for sole author | **you** if co-authors join |
| 15 | ORCID | unknown | present | present | yes |
| 16 | Ethics | unknown | stated | present | yes |
| 17 | Reference style | unknown | Quarto default (Chicago author-date via citeproc) with DOIs | DOIs on every reference that has one | check against journal CSL |
| 18 | DOIs must resolve | checked automatically (your brief) | every DOI checked | see log | see REPORT |
| 19 | File format | PDF upload plus repo link (your brief) | `paper/manuscript.pdf` | rendered | yes |
| 20 | Licence of software | likely open | GPL-2.0-or-later, stated | present | yes |

If the journal supplies a CSL file or LaTeX template, add `csl:` to the YAML header and rerender. Nothing else in the source would need to change.
