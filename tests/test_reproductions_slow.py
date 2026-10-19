"""The six published reproductions and the side-by-side analyses, end to end.

Each script compares every quantity with the stored R results and exits
non-zero on any failure. They download the authors' data on first use and
take about an hour in all, so they run only when asked:

    DLNMPY_SLOW=1 pytest -m slow
"""

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ["lancet2015_vs_r.py", "bmcmrm2013_vs_r.py", "statmed2014.py", "leone2014.py",
           "vicedo2019.py", "biomet2017.py"]

pytestmark = [pytest.mark.slow,
              pytest.mark.skipif(os.environ.get("DLNMPY_SLOW") != "1",
                                 reason="set DLNMPY_SLOW=1 to run the full reproductions")]


@pytest.mark.parametrize("script", SCRIPTS)
def test_published_reproduction(script):
    out = subprocess.run([sys.executable, script], cwd=ROOT / "examples", capture_output=True, text=True)
    assert out.returncode == 0, out.stdout[-3000:] + out.stderr[-3000:]
    last = out.stdout.strip().splitlines()[-1]
    assert re.search(r"comparisons with R, .* 0 failures$", last), last


def test_side_by_side():
    out = subprocess.run([sys.executable, "tools/side_by_side.py"], cwd=ROOT, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr[-3000:]
    assert out.stdout.strip().splitlines()[-1].endswith(" 0 failures")
