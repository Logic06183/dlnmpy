"""Runs the end-to-end side-by-side validation (tools/side_by_side.py) as a
test: every analysis is fitted in Python with statsmodels and compared with
the R results stored in tests/side_by_side/r_results.json."""

import runpy
from pathlib import Path

import pytest

pytest.importorskip("statsmodels")


def test_side_by_side_has_no_failures(capsys):
    script = Path(__file__).resolve().parents[1] / "tools" / "side_by_side.py"
    with pytest.raises(SystemExit) as exc:
        runpy.run_path(str(script), run_name="__main__")
    out = capsys.readouterr().out
    assert exc.value.code == 0, out
    assert "0 failures" in out
