"""Run ../figures/rd_charts.py unchanged, but save each chart as PDF as well as PNG, for the paper.

rd_charts.py writes PNG only, into ../figures/out/ (not committed). This wrapper patches matplotlib's savefig so that
every PNG it writes also gets a PDF next to it, then runs the script with the arguments given.

  uv run --with matplotlib python rare-disease/paper/charts_pdf.py [--study indel|cancer]
"""
import runpy
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.figure  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "figures" / "rd_charts.py"
_save = matplotlib.figure.Figure.savefig


def savefig(self, fname, *args, **kwargs):
    _save(self, fname, *args, **kwargs)
    path = Path(fname)
    if path.suffix == ".png":
        _save(self, path.with_suffix(".pdf"), *args, **kwargs)


matplotlib.figure.Figure.savefig = savefig
sys.argv = [str(SCRIPT)] + sys.argv[1:]
runpy.run_path(str(SCRIPT), run_name="__main__")
