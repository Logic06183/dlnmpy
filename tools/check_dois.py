"""Check that every DOI in paper/references.bib resolves.

    python tools/check_dois.py

Asks the DOI handle service for each DOI (no publisher pages are fetched) and
exits non-zero if any DOI is unknown or the service cannot be reached. Uses the
standard library only.
"""
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

BIB = Path(__file__).resolve().parent.parent / "paper" / "references.bib"


def dois(text):
    for key, body in re.findall(r"@\w+\{([^,]+),(.*?)\n\}", text + "\n", flags=re.S):
        m = re.search(r"doi\s*=\s*\{([^}]+)\}", body)
        if m:
            yield key.strip(), m.group(1).strip()


def resolves(doi):
    url = "https://doi.org/api/handles/" + urllib.request.quote(doi, safe="/")
    try:
        with urllib.request.urlopen(url, timeout=30) as r:          # noqa: S310
            return json.load(r).get("responseCode") == 1, ""
    except urllib.error.HTTPError as e:
        return False, f"HTTP {e.code}"
    except OSError as e:
        return False, str(e)


def main():
    bad = 0
    for key, doi in dois(BIB.read_text(encoding="utf-8")):
        ok, why = resolves(doi)
        bad += not ok
        print(f"{'ok  ' if ok else 'FAIL'} {key:28s} {doi} {why}")
    print(f"{bad} DOIs did not resolve" if bad else "all DOIs resolve")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
