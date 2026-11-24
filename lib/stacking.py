"""Stacked bunch per segment (figure notebooks): calibration CSV loader, entry → bucket assignment, stacked-bunch fit.

- `load_calibration`: the calibration CSV (`sbnd_fine_calibration.csv`) with parsed times / lists + skewness `rl_asym`
- `apply_to_tree`: every entry of the ntuple → its chunk and nearest saved bucket (core + edge), residual to that bucket
- `stacked_bunch_fit`: residuals of many chunks stacked → symmetric + asymmetric Gaussian
"""

import ast

import numpy as np
import pandas as pd

from calib_fit_helpers import (PRODUCTION_INT_COLS, PRODUCTION_LIST_COLS, _fit_asym_gaussian, _fit_gaussian_windowed,
                               _gaussian_chi2ndof, asymmetric_gaussian)

STATUSES = ('in_spill', 'before_spill', 'after_spill', 'bad_time', 'no_chunk')
BAD_TIME = -9000   # opflash_t0_corrected_ns ≤ this = no corrected time (−9999)


def load_calibration(path, include_bad=False):
    """Calibration CSV → DataFrame by start; times / list columns parsed, floats read exactly (round trip).
    - `include_bad=False`: fitted chunks only (empty `skip_reason`)
    - Adds `rl_asym` = σR − σL and `rl_asym_err` (σL, σR errors in quadrature)
    """
    df = pd.read_csv(path, float_precision='round_trip')
    df['start'] = pd.to_datetime(df['start'])
    df['end'] = pd.to_datetime(df['end'])
    df['skip_reason'] = df['skip_reason'].fillna('')
    df['logs'] = df['logs'].fillna('')
    df['segment'] = df['segment'].astype('Int64')
    if not include_bad:
        df = df[df['skip_reason'] == ''].copy()
    for col in PRODUCTION_LIST_COLS:
        df[col] = df[col].apply(ast.literal_eval)
    for col in PRODUCTION_INT_COLS:
        df[col] = df[col].astype('int64' if df[col].notna().all() else 'Int64')
    df['rl_asym'] = df['sigma_r'] - df['sigma_l']
    df['rl_asym_err'] = np.hypot(df['sigma_l_err'], df['sigma_r_err'])
    return df.sort_values('start').reset_index(drop=True)


def all_buckets(row):
    """Core + edge buckets of one calibration CSV row → (slots, positions), position order, fixed frame (`window_offset`)."""
    core = np.arange(len(row['peak_positions']))
    slots = np.concatenate([core, row['edge_slots']]).astype(int) + int(row['window_offset'])
    pos = np.concatenate([row['peak_positions'], row['edge_positions']])
    o = np.argsort(pos)
    return slots[o], pos[o]


def assign_hits(t, row):
    """Entry times of one chunk → dict of arrays: status, slot, pos, resid.
    - Spill window = first bucket − P/2 ≤ t < last bucket + P/2; inside: nearest bucket (ties to the earlier), resid = t − pos
    - Outside: 'before_spill' / 'after_spill'; slot / pos / resid NaN
    """
    t = np.asarray(t, dtype=float)
    slots, pos = all_buckets(row)
    lo, hi = pos[0] - row['period'] / 2, pos[-1] + row['period'] / 2
    j = np.clip(np.searchsorted(pos, t), 1, len(pos) - 1)
    k = np.where(np.abs(t - pos[j - 1]) <= np.abs(t - pos[j]), j - 1, j)
    inside = (t >= lo) & (t < hi)
    status = np.where(t < lo, 'before_spill', 'after_spill').astype(object)
    status[inside] = 'in_spill'
    nan = np.full(len(t), np.nan)
    return {'status': status, 'slot': np.where(inside, slots[k], nan), 'pos': np.where(inside, pos[k], nan),
            'resid': np.where(inside, t - pos[k], nan)}


def apply_to_tree(raw_df, db):
    """Every entry of `raw_df` (`load_run_raw` + 'opflash_t0_corrected_ns') → copy with `calib_chunk` (row of `db`, −1 = none),
    `calib_status` (`STATUSES`), `calib_slot`, `calib_pos`, `calib_resid`.
    - `db`: fitted chunks (`load_calibration`); an entry belongs to a chunk if its run is in the row's runs and start ≤ ts ≤ end
    - Asserts each entry lands in at most one chunk
    """
    out = raw_df.copy()
    n = len(out)
    chunk = np.full(n, -1)
    cols = {k: np.full(n, np.nan) for k in ('slot', 'pos', 'resid')}
    status = np.full(n, 'no_chunk', dtype=object)
    run, ts, t = out['run'].values, out['ts'].values, out['bunch_corrected'].values
    bad = out['opflash_t0_corrected_ns'].values <= BAD_TIME
    for i, row in db.iterrows():
        m = np.isin(run, row['run']) & (ts >= row['start'].to_datetime64()) & (ts <= row['end'].to_datetime64())
        assert not (chunk[m] >= 0).any(), f'entry in two chunks (row {i})'
        chunk[m] = i
        status[m & bad] = 'bad_time'
        g = m & ~bad
        a = assign_hits(t[g], row)
        status[g] = a['status']
        for k in cols:
            cols[k][g] = a[k]
    out['calib_chunk'] = chunk
    out['calib_status'] = pd.Categorical(status, categories=STATUSES)
    for k, v in cols.items():
        out[f'calib_{k}'] = v
    return out


def _asym_gaussian_chi2ndof(bins_mid, raw_counts, bin_width, popt, mask):
    """Neyman chi2/ndof of `asymmetric_gaussian` vs raw counts in `mask`; model scaled by the full histogram sum."""
    model_counts = asymmetric_gaussian(bins_mid[mask], *popt) * raw_counts.sum() * bin_width
    chi2 = np.sum(((raw_counts[mask] - model_counts) ** 2) / np.clip(raw_counts[mask], 1, None))
    return chi2 / max(mask.sum() - 5, 1)


def stacked_bunch_fit(resid, period, window_sigma=2.5, nbins=77):
    """Stacked residuals of one segment → sym + asym Gaussian in a ±`window_sigma` window (2.5σ: 3σ pulls the fit into the tails).
    - Returns dict: `mid`, `bw`, `n` (counts), `mask` (fit window), `xx` (curve x), `k` (density → entries / bin),
      `ps` / `es` (sym popt / errors), `chi_s`, `pa` / `ea` (asym popt / errors, σL = pa[3], σR = pa[4]), `chi_a`
    """
    b = np.linspace(-period / 2, period / 2, nbins)
    mid, bw = (b[1:] + b[:-1]) / 2, b[1] - b[0]
    h, _ = np.histogram(resid, b, density=True)
    n, _ = np.histogram(resid, b)
    ps, cs, mask = _fit_gaussian_windowed(mid, h, n, [0, h.max(), 0, 2], window_sigma)
    pa, ca = _fit_asym_gaussian(mid, h, n, mask, ps)
    return dict(mid=mid, bw=bw, n=n, mask=mask, xx=np.linspace(mid[mask][0], mid[mask][-1], 400), k=len(resid) * bw,
                ps=ps, es=np.sqrt(np.diag(cs)), chi_s=_gaussian_chi2ndof(mid, n, bw, ps, mask),
                pa=pa, ea=np.sqrt(np.diag(ca)), chi_a=_asym_gaussian_chi2ndof(mid, n, bw, pa, mask))
