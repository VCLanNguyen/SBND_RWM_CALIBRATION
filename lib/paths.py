"""Input / output paths (str, absolute) -- every file path goes through here.

- Inputs: one TimingCalibAna ntuple (module v2) per run, any readable path (local or Fermilab, e.g. /exp/sbnd/data/...);
  set by `set_paths` (calibration notebook, Setup)
- Outputs: `OUT/run1/`, `OUT/run2/` = calibration CSV `sbnd_fine_calibration.csv` + lean-entry cache
- CSV output: `write_csv(df, path)` -- one timestamp format everywhere (pandas `to_csv` picks it per chunk of rows)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = str(Path(__file__).resolve().parents[1])

RAW_TREE = {1: None, 2: None}    # TimingCalibAna ntuple per run (`set_paths`)

OUT = f'{ROOT}/output'           # `set_paths(out=)`
RUN1_SBND = f'{OUT}/run1'
RUN2_SBND = f'{OUT}/run2'


def set_paths(ntuple, out=None):
    """`ntuple` = {run: path of the run's ntuple}; `out` = output folder (default `<repo>/output`); creates `out/run1`, `out/run2`."""
    global OUT, RUN1_SBND, RUN2_SBND
    for run, p in ntuple.items():
        if not p or not os.path.isfile(p):
            raise FileNotFoundError(f'Run {run} ntuple not found: {p!r} -- set NTUPLE in the notebook Setup cell')
        RAW_TREE[run] = os.path.abspath(p)
    OUT = os.path.abspath(out) if out else OUT
    RUN1_SBND, RUN2_SBND = f'{OUT}/run1', f'{OUT}/run2'
    for d in (RUN1_SBND, RUN2_SBND):
        os.makedirs(d, exist_ok=True)


def write_csv(df, path):
    """`df.to_csv(path, index=False)` with every datetime column (datetime64 or object of Timestamps) as 'YYYY-MM-DD HH:MM:SS.fffffffff',
    empty for NaT; other columns unchanged."""
    out = df.copy()
    for c in out.columns:
        s = out[c]
        if s.dtype == object and len(s.dropna()) and s.dropna().map(lambda v: isinstance(v, pd.Timestamp)).all():
            s = pd.to_datetime(s)
        if pd.api.types.is_datetime64_any_dtype(s):
            txt = np.char.replace(np.datetime_as_string(s.values.astype('datetime64[ns]'), unit='ns'), 'T', ' ')
            out[c] = np.where(s.isna().values, '', txt)
    out.to_csv(path, index=False)
