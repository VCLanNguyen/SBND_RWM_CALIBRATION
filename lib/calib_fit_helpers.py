"""Fitting and data-loading helpers of the SBND RWM timing calibration.

- Ntuple loader (`load_run_raw` / `load_raw_df`), time chunking (`build_fine_chunks`), bucket finding, period / intercept / σ fits
- Calibration CSV columns + `calibration_table`
- Imported by `calibration/calibration.ipynb`: `from calib_fit_helpers import *` (`lib/` on `sys.path`)
"""

import os
import uproot
import numpy as np
import pandas as pd

from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from scipy.stats import poisson

import paths

# fit config
SIGMA_FIT_WINDOW = 3
PERIOD_WEIGHTED = True
POISSON_WEIGHT = True

RWM_MAX_UINT64 = np.iinfo(np.uint64).max   # invalid-timestamp default on sbnobj develop (gen2: 0)

# edge buckets: partly filled buckets just outside the 81-bucket window
EDGE_MAX_SLOTS = 2          # per side
EDGE_MAX_PVALUE = 1e-4      # SBND: Poisson excess over the background; ≈ 600 edge slots per run → ≈ 0.06 false
EDGE_SIG_HALF = 4           # SBND significance window: ±4 bins (1-ns bins) of the period-line slot
EDGE_KEYS = ('edge_slots', 'edge_positions', 'edge_pos_err')
# weak buckets: the slot just past the outermost saved bucket on each side (real when summed, not significant per chunk)
WEAK_KEYS = ('weak_slots', 'weak_positions', 'weak_pvalue')
SKEW_KEYS = ('sigma_l', 'sigma_l_err', 'sigma_r', 'sigma_r_err')   # asym fit to the stacked bunch (`fit_bunch_skew`)


RAW_COLUMNS = ('run', 'subrun', 'event', 'slice_id', 'opflash_t0_ns', 'nu_tof_light_ns', 'nu_tof_charge_ns',
               'opflash_t0_corrected_ns', 'frame_tdc_rwm_ns', 'frame_hlt_beam_gate_ns', 'tdc_rwm_ns',
               'hlt_beam_gate_ns')   # pre-v2 tree branches
SLICE_KEY = ['run', 'subrun', 'event', 'slice_id']

CACHE_COLUMNS = RAW_COLUMNS + ('particle_prop_time_ns', 'nu_score', 'bfm_score', 'nu_vtx_x', 'nu_vtx_y', 'nu_vtx_z')


def raw_cache_path(run):
    """Lean-entry cache of the run's ntuple: `<OUT>/run<N>/<ntuple name>_lean.pkl` (output folder, so the input may be read-only)."""
    name = os.path.basename(paths.RAW_TREE[run]).rsplit('.', 1)[0]
    return f"{(paths.RUN1_SBND, paths.RUN2_SBND)[run - 1]}/{name}_lean.pkl"


def load_run_raw(run, columns=RAW_COLUMNS, cache=True, **kw):
    """`load_raw_df` of the run's ntuple (`paths.RAW_TREE`).
    - `cache`: `load_raw_df` rows of `CACHE_COLUMNS` pickled once (`raw_cache_path`), read from there (~1 s vs ~15 s per run)
      - rebuilt when older than the ntuple or missing a requested column; columns outside `CACHE_COLUMNS` or
        extra `kw` → direct ntuple read
    """
    cols = list(dict.fromkeys(columns))
    if not cache or kw or not set(cols) <= set(CACHE_COLUMNS):
        return load_raw_df(paths.RAW_TREE[run], columns=cols, **kw)
    p = raw_cache_path(run)
    df = pd.read_pickle(p) if os.path.exists(p) and os.path.getmtime(p) > os.path.getmtime(paths.RAW_TREE[run]) else None
    if df is None or not set(cols) <= set(df.columns):
        df = load_raw_df(paths.RAW_TREE[run], columns=CACHE_COLUMNS)
        df.to_pickle(p)
    return df[cols + ['bunch_corrected', 'ts']]


def load_raw_df(path, columns=RAW_COLUMNS, step_size='200 MB'):
    """Load `timingCalibTree` into a per-entry df.
    - Reads only `columns` (the module-v2 tree has 71 branches), filtered per read step to keep memory low
    - Keeps status 0 only (`corr_flash_status`: entry with a corrected flash); trees without the branch are not cut
    - Drops rows with an invalid RWM: `tdc_rwm_ns` = 0 (gen2 default) or max uint64 (develop default)
    - De-duplicates on (run, subrun, event, slice_id), first row kept
    - Adds `bunch_corrected` (front-face corrected) and `ts` (from `tdc_rwm_ns`)
    """
    tree = uproot.open(path)['timingcalib/timingCalibTree']
    has_status = 'corr_flash_status' in tree.keys()
    extra = ['corr_flash_status'] if has_status else []
    read = list(dict.fromkeys(list(columns) + extra))
    parts = []
    for df in tree.iterate(read, library='pd', step_size=step_size):
        keep = (df['tdc_rwm_ns'] != 0) & (df['tdc_rwm_ns'] != RWM_MAX_UINT64)
        if has_status:
            keep &= df['corr_flash_status'] == 0
        parts.append(df[keep])
    df = pd.concat(parts)
    df = df.drop(columns=[c for c in extra if c not in columns])
    df = df.drop_duplicates(SLICE_KEY)

    df['bunch_corrected'] = df['opflash_t0_corrected_ns'] - df['nu_tof_charge_ns']
    df['ts'] = pd.to_datetime(df['tdc_rwm_ns'], unit='ns')
    return df


def linear(x, m, b):
    return m*x + b

def gaussian(x, H, a, x0, sigma):
    return H + a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def asymmetric_gaussian(x, H, a, x0, sigma_left, sigma_right):
    """Gaussian with `sigma_left` below x0, `sigma_right` above."""
    sigma = np.where(x < x0, sigma_left, sigma_right)
    return H + a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def find_dynamic_buckets(hist, centers, expected_n=81, min_distance=9, window=8):
    """Find `expected_n` evenly spaced buckets in `hist`.
    - Height-cut search over loosening steps; median spacing merges double detections
    - Exactly `expected_n` found is kept only if evenly spaced (`_regular_spacing`), else the comb decides
    - 84 RF slots, 81 filled: weak satellites sit just outside, so one height cut fails at the margin
    - Returns the `expected_n` contiguous comb slots with the most counts (`_select_best_window`)
    - No gap-fill: `_require_good_buckets` raises instead
    """
    noise = np.std(hist)

    best_under = None
    over = None
    buckets = None

    for step in [0.15, 0.2, 0.3, 0.5, 0.8, 1, 1.3, 1.5, 2, 2.5, 3]:
        found, _ = find_peaks(hist, height=step*noise, distance=min_distance)
        n = len(found)
        if n == expected_n and buckets is None:
            buckets = found
        if n < expected_n:
            if best_under is None or n > best_under[0]:
                best_under = (n, found)
        else:
            # smallest over-detection: loose cuts add noise bumps (124 vs 82 real case)
            if over is None or n < over[0]:
                over = (n, found)

    if buckets is None:
        if over is not None:
            buckets = over[1]
        elif best_under is not None and best_under[0] >= 2:
            buckets = best_under[1]
        else:
            return best_under[1] if best_under is not None else np.array([], dtype=int)

    if len(buckets) < 2:
        return buckets

    # merge double detections of one peak
    spacing_guess = np.median(np.diff(centers[buckets]))
    merged = [buckets[0]]
    for p in buckets[1:]:
        if centers[p] - centers[merged[-1]] < 0.6 * spacing_guess:
            if hist[p] > hist[merged[-1]]:
                merged[-1] = p
        else:
            merged.append(p)
    buckets = np.array(merged)

    # an exact count can still be wrong: a missed bucket + a side bump (Run 2 run 20040, gaps 37 / 23 + 16 ns)
    if len(buckets) == expected_n and _regular_spacing(centers[buckets]):
        return buckets

    comb = _select_best_window(hist, centers, buckets, expected_n, window)
    if comb is not None:
        return comb

    # fallback: height-based trim
    if len(buckets) > expected_n:
        keep = np.argsort(hist[buckets])[::-1][:expected_n]
        return np.array(sorted(buckets[keep]))
    return buckets


def _regular_spacing(positions, lo=0.6, hi=1.4):
    """True if every gap between neighbouring `positions` is within [`lo`, `hi`] x the median gap."""
    gaps = np.diff(positions)
    med = np.median(gaps)
    return bool(np.all((gaps >= lo * med) & (gaps <= hi * med)))


def _select_best_window(hist, centers, buckets, expected_n, window, margin=4):
    """Comb at the median spacing, `margin` slots past each end; the `expected_n` contiguous slots with the most counts, or None."""
    spacing = np.median(np.diff(centers[buckets]))
    if not np.isfinite(spacing) or spacing <= 0:
        return None

    grid = centers[buckets[0]] - margin * spacing + spacing * np.arange(
        len(buckets) + 2 * margin)

    slots = []
    for pos in grid:
        nbhd = np.where(np.abs(centers - pos) < 0.4 * spacing)[0]
        if len(nbhd) == 0:
            continue
        best = nbhd[np.argmax(hist[nbhd])]
        count = hist[max(0, best - window):best + window].sum()
        slots.append((best, count))

    if len(slots) < expected_n:
        return None

    idxs = np.array([s[0] for s in slots])
    counts = np.array([s[1] for s in slots], dtype=float)
    window_sums = np.array([counts[i:i + expected_n].sum()
                             for i in range(len(counts) - expected_n + 1)])
    best_start = np.argmax(window_sums)
    return idxs[best_start:best_start + expected_n]


def fit_edge_buckets(hist, centers, positions, line, H, window=8, poisson_weight=False,
                     max_slots=EDGE_MAX_SLOTS, max_pvalue=EDGE_MAX_PVALUE, sig_half=EDGE_SIG_HALF,
                     max_chi2ndof=20, max_pos_err=5, max_sigma=5, max_shift=0.3):
    """Partly filled buckets just outside the core window → dict of `EDGE_KEYS` lists, in slot order (empty if none).
    - Slots −1, −2 … / n, n + 1 … (fit-index frame, n = len(`positions`)) at the period line `line` = (period, intercept)
    - Walk outward, stop at the first slot that fails: counts within ±`sig_half` bins a Poisson excess over 2·`sig_half`·`H`
      (p < `max_pvalue`), fixed-`H` fit (start x0 = line prediction, so a sub-peak can't capture it) passes the core
      guards + within ±`max_shift` period of the line
    - Significance, not a fixed occupancy / count cut: scales with chunk size (thin buckets in small chunks)
    - Core fit untouched; a failing edge never skips the chunk
    """
    period, intercept = line
    found = []
    for side in (np.arange(-1, -max_slots - 1, -1), np.arange(len(positions), len(positions) + max_slots)):
        for slot in side:
            pred = intercept + period * slot
            p = int(np.argmin(np.abs(centers - pred)))
            if p - window < 0 or p + window > len(hist):
                break
            pval = float(poisson.sf(hist[p - sig_half:p + sig_half].sum() - 1, 2 * sig_half * H))
            if pval >= max_pvalue:
                break
            try:
                _, _, popts, pcovs, chi2s = fit_peaks(hist, centers, [p], window=window, poisson_weight=poisson_weight,
                                                      fixed_H=H, x0_start=[pred])
            except (RuntimeError, ValueError):
                break
            popt, err, sig = popts[0], np.sqrt(pcovs[0][2, 2]), np.abs(popts[0][3:]).max()
            if not (np.isfinite(err) and err <= max_pos_err and chi2s[0] <= max_chi2ndof and sig <= max_sigma
                    and abs(popt[2] - pred) <= max_shift * period):
                break
            found.append((int(slot), float(popt[2]), float(err)))
    found.sort()
    return {k: [f[i] for f in found] for i, k in enumerate(EDGE_KEYS)}


def weak_buckets(hist, centers, positions, line, H, edge_slots, window=8, sig_half=EDGE_SIG_HALF):
    """Weak buckets = the slot just past the outermost saved bucket on each side → dict of `WEAK_KEYS` lists.
    - Position = period-line prediction; p = Poisson excess within ±`sig_half` bins over 2·`sig_half`·`H`
    - Recorded only (not fitted, not stacked); slots outside the histogram are left out
    """
    period, intercept = line
    n = len(positions)
    lo = min([s for s in edge_slots if s < 0], default=0) - 1
    hi = max([s for s in edge_slots if s >= n], default=n - 1) + 1
    weak = []
    for slot in (lo, hi):
        pred = intercept + period * slot
        p = int(np.argmin(np.abs(centers - pred)))
        if p - window < 0 or p + window > len(hist):
            continue
        n4 = int(hist[p - sig_half:p + sig_half].sum())
        weak.append((int(slot), float(pred), float(poisson.sf(n4 - 1, 2 * sig_half * H))))
    return {k: [w[i] for w in weak] for i, k in enumerate(WEAK_KEYS)}


def _require_good_buckets(hist, centers, expected_n=81, min_bucket_counts=15, window=8):
    """`expected_n` detected buckets, each >= `min_bucket_counts` entries, else ValueError.
    - Thin buckets fit falsely tight and wrong (3 counts: chi2/ndof ~0, position off by 112 ns)
    """
    buckets = find_dynamic_buckets(hist, centers, expected_n=expected_n)
    if len(buckets) != expected_n:
        raise ValueError(f'found {len(buckets)}/{expected_n} genuine buckets -- need more stats')
    counts = np.array([hist[max(0, p - window):p + window].sum() for p in buckets])
    if counts.min() < min_bucket_counts:
        raise ValueError(f'thinnest bucket has {counts.min()} raw counts '
                          f'(< {min_bucket_counts}) -- need more stats')
    return buckets


def fit_peaks(hist, centers, buckets, window=8, window_sigma=3, poisson_weight=False, fixed_H=None, x0_start=None):
    """Per-bucket Gaussian fits within ±`window` bins of each bucket.
    - Both passes: fit over ±`window`, refit within `window_sigma`·sigma so bucket tails can't pull; full-window fit kept if
      < 8 bins or the refit fails (narrow buckets, ~0.8 % pass 1 / 0.4 % pass 2)
    - `fixed_H`: H fixed (a, x0, sigma free); removes the H/amplitude degeneracy (wide/flat fits)
    - `poisson_weight`: sqrt(N) weights
    - `x0_start`: start x0 per bucket (None = tallest bin)
    - Returns bucket x, y, popts, pcovs, raw-count chi2/ndof, in `buckets` order
    """
    bucket_x, bucket_y, popts, pcovs, chi2ndofs = [], [], [], [], []
    for i, p in enumerate(buckets):
        x = centers[p-window:p+window]
        y = hist[p-window:p+window]
        x0 = x[np.argmax(y)] if x0_start is None else x0_start[i]
        if fixed_H is None:
            guess = [0, max(y), x0, 2]
            popt, pcov, mask = _fit_gaussian_windowed(x, y, y, guess, window_sigma, poisson_weight=poisson_weight)
        else:
            popt, pcov, mask = _fit_fixed_H_windowed(x, y, fixed_H, window_sigma, poisson_weight=poisson_weight, x0_start=x0)
        model = gaussian(x[mask], *popt)
        chi2 = np.sum(((y[mask] - model) ** 2) / np.clip(y[mask], 1, None))
        ndof = max(mask.sum() - (3 if fixed_H is not None else 4), 1)
        bucket_x.append(x)
        bucket_y.append(y)
        popts.append(popt)
        pcovs.append(pcov)
        chi2ndofs.append(chi2 / ndof)
    return bucket_x, bucket_y, popts, pcovs, np.array(chi2ndofs)


def _fit_fixed_H_windowed(x, y, fixed_H, window_sigma, poisson_weight=False, min_window_bins=8, x0_start=None):
    """`_fit_gaussian_windowed` with H fixed: fit a, x0, sigma, then refit within `window_sigma`·sigma of the center.
    - Full-window fit if < `min_window_bins` or the refit fails; returns 4-parameter (popt, pcov, mask), pcov row/col 0 = 0
    """
    def gaussian_fixed(x, a, x0, sigma):
        return fixed_H + a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    w = np.sqrt(np.clip(y, 1, None))
    kw = dict(sigma=w, absolute_sigma=True) if poisson_weight else {}
    x0 = x[np.argmax(y)] if x0_start is None else x0_start
    popt3, pcov3 = curve_fit(gaussian_fixed, x, y, p0=[max(y) - fixed_H, x0, 2], **kw)
    mask = np.ones(len(x), dtype=bool)
    if window_sigma:
        m = np.abs(x - popt3[1]) <= window_sigma * abs(popt3[2])
        if m.sum() >= min_window_bins:
            try:
                kw_m = dict(sigma=w[m], absolute_sigma=True) if poisson_weight else {}
                popt3, pcov3 = curve_fit(gaussian_fixed, x[m], y[m], p0=[max(y[m]) - fixed_H, popt3[1], abs(popt3[2])],
                                         **kw_m)
                mask = m
            except Exception:
                pass
    pcov = np.zeros((4, 4))
    pcov[1:, 1:] = pcov3
    return np.array([fixed_H, *popt3]), pcov, mask


def _merge_fit_passes(popts1, pcovs1, chi2ndofs1, popts2, pcovs2, chi2ndofs2, max_pos_err=5):
    """Per bucket, take pass 2 (fixed H) unless its position error > `max_pos_err` and pass 1's isn't.
    - Rare: local background differs from the chunk floor (29 counts: 0.50 → 406 ns)
    - Threshold, not min: pass 2 is sometimes marginally larger in normal chunks
    """
    popts, pcovs, chi2ndofs = [], [], []
    for p1, c1, x1, p2, c2, x2 in zip(popts1, pcovs1, chi2ndofs1, popts2, pcovs2, chi2ndofs2):
        err1, err2 = np.sqrt(c1[2, 2]), np.sqrt(c2[2, 2])
        use1 = (not np.isfinite(err2) or err2 > max_pos_err) and np.isfinite(err1) and err1 <= max_pos_err
        popts.append(p1 if use1 else p2)
        pcovs.append(c1 if use1 else c2)
        chi2ndofs.append(x1 if use1 else x2)
    return popts, pcovs, np.array(chi2ndofs)


def _require_good_fits(chi2ndofs, pcovs, popts, max_chi2ndof=20, max_pos_err=5, max_sigma=5):
    """Three independent guards; ValueError names the worst bucket.
    - chi2/ndof > 20: poor fit with enough counts (max 16.25 over 22113 fits)
    - position error > 5 ns: degenerate covariance (27779 ns seen; healthy 0.15-0.5)
    - sigma > 5 ns: wide/flat fit (healthy max 4.56 vs degenerate 5.27-34.6)
    """
    worst = int(np.argmax(chi2ndofs))
    if chi2ndofs[worst] > max_chi2ndof:
        raise ValueError(f'bucket {worst} fit chi2/ndof={chi2ndofs[worst]:.1f} '
                          f'(> {max_chi2ndof}) -- poor fit despite adequate stats')

    pos_errs = np.array([np.sqrt(np.diag(pcov)[2]) for pcov in pcovs])
    finite = np.isfinite(pos_errs)
    worst_pos = int(np.argmax(np.where(finite, pos_errs, np.inf)))
    if not finite[worst_pos] or pos_errs[worst_pos] > max_pos_err:
        raise ValueError(f'bucket {worst_pos} position error={pos_errs[worst_pos]:.3g}ns '
                          f'(non-finite or > {max_pos_err}) -- degenerate covariance despite adequate stats')

    sigmas = np.array([np.abs(p[3:]).max() for p in popts])
    worst_sigma = int(np.argmax(sigmas))
    if sigmas[worst_sigma] > max_sigma:
        raise ValueError(f'bucket {worst_sigma} fitted sigma={sigmas[worst_sigma]:.3g}ns '
                          f'(> {max_sigma}) -- degenerate wide/flat fit despite fine chi2/ndof and position error')


def peak_by_peak_gaussian_fit(data, bins, weighted_period=False, poisson_weight=False):
    """Period line from 81 per-bucket fits, 2 passes.
    - Pass 1: free H, seeds positions only
    - H: background of the stacked bunch (`_stacked_background_H`)
    - Pass 2: fixed H; merged per bucket (`_merge_fit_passes`), then `_require_good_fits`
    - Returns period line (popt = (period, value at fit index 0), pcov), period-fit chi2/ndof, positions + errors, H
    """
    centers = (bins[:-1] + bins[1:]) / 2
    hist, _, = np.histogram(data, bins = bins )

    buckets = _require_good_buckets(hist, centers)

    _, _, popts1, pcovs1, chi2ndofs1 = fit_peaks(hist, centers, buckets, poisson_weight=poisson_weight)
    H = _stacked_background_H(data, np.array([p[2] for p in popts1]))

    _, _, popts2, pcovs2, chi2ndofs2 = fit_peaks(hist, centers, buckets, poisson_weight=poisson_weight, fixed_H=H)
    popts, pcovs, chi2ndofs = _merge_fit_passes(popts1, pcovs1, chi2ndofs1, popts2, pcovs2, chi2ndofs2)
    _require_good_fits(chi2ndofs, pcovs, popts)

    xp = [p[2] for p in popts]
    yp = np.arange(0,len(xp), 1)
    peak_pos_err = np.array([np.sqrt(np.diag(p)[2]) for p in pcovs])
    if weighted_period:
        linear_popt, linear_pcov = curve_fit(linear, yp, xp, sigma=peak_pos_err, absolute_sigma=True)
    else:
        linear_popt, linear_pcov = curve_fit(linear, yp, xp)

    period_resid = np.array(xp) - linear(yp, *linear_popt)
    period_chi2 = np.sum((period_resid / peak_pos_err) ** 2)
    period_chi2_ndof = period_chi2 / (len(xp) - 2)

    return linear_popt, linear_pcov, period_chi2_ndof, np.array(xp), peak_pos_err, H


def _fit_gaussian_windowed(bins_mid, hist, raw_counts, guess, window_sigma, min_window_bins=8, poisson_weight=False):
    """Gaussian fit, then refit within `window_sigma`·sigma of the center; returns (popt, pcov, mask).
    - Keeps the bucket tails / histogram edge out of the fit; full fit if < `min_window_bins`
    - Poisson weights from `raw_counts`, scaled to the density `hist`
    """
    if poisson_weight:
        nz = raw_counts > 0
        scale = np.median(hist[nz] / raw_counts[nz]) if nz.any() else 1.0
        poisson_sigma = np.sqrt(np.clip(raw_counts, 1, None)) * scale
        fit_kwargs = dict(sigma=poisson_sigma, absolute_sigma=True)
    else:
        fit_kwargs = {}
    popt0, pcov0 = curve_fit(gaussian, bins_mid, hist, p0=guess, **fit_kwargs)
    full_mask = np.ones(len(bins_mid), dtype=bool)
    if not window_sigma:
        return popt0, pcov0, full_mask

    center0, sigma0 = popt0[2], abs(popt0[3])
    mask = np.abs(bins_mid - center0) <= window_sigma * sigma0
    if mask.sum() < min_window_bins:
        return popt0, pcov0, full_mask

    try:
        w_kwargs = dict(sigma=poisson_sigma[mask], absolute_sigma=True) if poisson_weight else {}
        popt_w, pcov_w = curve_fit(gaussian, bins_mid[mask], hist[mask],
                                    p0=[guess[0], max(hist[mask]), center0, sigma0], **w_kwargs)
        return popt_w, pcov_w, mask
    except Exception:
        return popt0, pcov0, full_mask


def _gaussian_chi2ndof(bins_mid, raw_counts, bin_width, popt, mask):
    """Neyman chi2/ndof vs raw counts in `mask`; model scaled by the full histogram sum."""
    model_counts = gaussian(bins_mid[mask], *popt) * raw_counts.sum() * bin_width
    chi2 = np.sum(((raw_counts[mask] - model_counts) ** 2) / np.clip(raw_counts[mask], 1, None))
    ndof = max(mask.sum() - 4, 1)
    return chi2 / ndof


def _fit_asym_gaussian(bins_mid, hist, raw_counts, mask, sym_popt, poisson_weight=False):
    """`asymmetric_gaussian` on the sym fit's mask and units, seeded from `sym_popt`."""
    guess = [sym_popt[0], sym_popt[1], sym_popt[2], abs(sym_popt[3]), abs(sym_popt[3])]
    if poisson_weight:
        nz = raw_counts > 0
        scale = np.median(hist[nz] / raw_counts[nz]) if nz.any() else 1.0
        poisson_sigma = np.sqrt(np.clip(raw_counts[mask], 1, None)) * scale
        fit_kwargs = dict(sigma=poisson_sigma, absolute_sigma=True)
    else:
        fit_kwargs = {}
    popt, pcov = curve_fit(asymmetric_gaussian, bins_mid[mask], hist[mask], p0=guess, **fit_kwargs)
    return popt, pcov


def _stacked_residual_histogram(data, positions, window_sigma=3, nbins=39, min_gap_frac=0.5, spill_period=None):
    """Residual of each entry to its nearest bucket, histogrammed.
    - Entries kept: `spill_period=None` → first … last bucket position (outer halves of the outermost buckets dropped; step-4 H);
      `spill_period` = P → first − P/2 ≤ t < last + P/2, ties to the earlier bucket
    - ValueError: non-finite / crossed positions, or a gap < `min_gap_frac` × median spacing (healthy min 0.875)
    - ValueError: binning drops an in-range residual
    - Returns (obins_mid, ohist, oraw_counts, obin_width, xp, resid)
    """
    xp = np.sort(np.asarray(positions))

    if not np.all(np.isfinite(xp)):
        raise ValueError('stacked bunch: non-finite bucket position(s) -- degenerate fit')

    gaps = np.diff(xp)
    period_est = np.median(gaps)
    worst = int(np.argmin(gaps))
    if gaps[worst] <= 0:
        raise ValueError(f'stacked bunch: bucket positions {worst}/{worst+1} not strictly '
                          f'increasing (gap={gaps[worst]:.3g}ns) -- duplicate/crossed fit')
    if gaps[worst] < min_gap_frac * period_est:
        raise ValueError(f'stacked bunch: gap between buckets {worst}/{worst+1} is '
                          f'{gaps[worst]:.3g}ns (< {min_gap_frac}x the typical '
                          f'{period_est:.3g}ns spacing) -- degenerate/wobbled bucket '
                          f'fit would corrupt the nearest-neighbor assignment boundary')

    data_arr = np.asarray(data)
    idx = np.clip(np.searchsorted(xp, data_arr), 1, len(xp) - 1)
    left, right = xp[idx - 1], xp[idx]
    if spill_period is None:
        nearest = np.where(np.abs(data_arr - left) < np.abs(data_arr - right), left, right)
        in_range = (data_arr >= xp[0]) & (data_arr <= xp[-1])
    else:
        nearest = np.where(np.abs(data_arr - left) <= np.abs(data_arr - right), left, right)
        in_range = (data_arr >= xp[0] - spill_period / 2) & (data_arr < xp[-1] + spill_period / 2)
    resid_all = data_arr - nearest
    resid = resid_all[in_range]

    half_window = max(np.abs(resid).max(), period_est / 2)
    obins = np.linspace(-half_window, half_window, nbins)
    obins_mid = np.convolve(obins, [0.5, 0.5], "valid")
    ohist, _ = np.histogram(resid, bins=obins, density=True)
    oraw_counts, _ = np.histogram(resid, bins=obins, density=False)
    obin_width = obins[1] - obins[0]

    if oraw_counts.sum() != resid.size:
        raise ValueError(f'stacked bunch: histogram captured {int(oraw_counts.sum())}/'
                          f'{resid.size} in-range residuals -- binning dropped '
                          f'{resid.size - int(oraw_counts.sum())} point(s)')

    return obins_mid, ohist, oraw_counts, obin_width, xp, resid


def _stacked_background_H(data, xp, window_sigma=3, nbins=39, target_bin_width=1.0):
    """Background floor H per bucket per `target_bin_width` bin, for the fixed-H refits.
    - Linear fit to the stacked-bunch bins outside the `window_sigma` core (flat cosmic rate)
    """
    obins_mid, ohist, oraw_counts, obin_width, xp, resid = \
        _stacked_residual_histogram(data, xp, window_sigma, nbins)

    # locates the core only; window_sigma=None still needs an off-peak region
    oguess = [0, max(ohist), 0, 2]
    _, _, omask = _fit_gaussian_windowed(obins_mid, ohist, oraw_counts, oguess, window_sigma or 3)

    # raw counts, so H is in fixed_H units
    mask_out = ~omask
    x_bg, y_bg = obins_mid[mask_out], oraw_counts[mask_out]
    lpopt, _ = curve_fit(linear, x_bg, y_bg,
                          sigma=np.sqrt(np.clip(y_bg, 1, None)), absolute_sigma=True)
    H_pooled = linear(0, *lpopt)

    # pooled -> per bucket: 81 filled buckets, not len(xp)
    return H_pooled * (target_bin_width / obin_width) / 81


def fit_bunch_sigma(data, positions, window_sigma=3, nbins=39, spill_period=None):
    """Bunch width σ: every entry's residual to its nearest bucket (`positions`), stacked, symmetric Gaussian fit.
    - `window_sigma`: core refit within ±`window_sigma`·σ (`_fit_gaussian_windowed`)
    - `spill_period`: entry window of `_stacked_residual_histogram` (None = first … last bucket position)
    - Returns (sigma, sigma_err, chi2ndof)
    """
    xp = np.sort(np.asarray(positions))
    obins_mid, ohist, oraw_counts, obin_width, xp, resid = \
        _stacked_residual_histogram(data, xp, window_sigma, nbins, spill_period=spill_period)

    oguess = [0, max(ohist), 0, 2]
    opopt, opcov, omask = _fit_gaussian_windowed(obins_mid, ohist, oraw_counts, oguess, window_sigma)
    sigma = abs(opopt[3])
    sigma_err = np.sqrt(np.diag(opcov)[3])
    sigma_chi2ndof = _gaussian_chi2ndof(obins_mid, oraw_counts, obin_width, opopt, omask)
    return sigma, sigma_err, sigma_chi2ndof


def fit_bunch_skew(data, positions, window_sigma=3, nbins=39, poisson_weight=False, spill_period=None):
    """Bunch skewness: the stacked bunch of `fit_bunch_sigma` fitted with `asymmetric_gaussian` → `sigma_l` / `sigma_r` (+ `_err`).
    - Poisson-weighted symmetric fit sets the core window and the start values; the asym fit uses the same window
    - NaN if a fit fails (the chunk is never skipped for it)
    """
    try:
        obins_mid, ohist, oraw_counts, obin_width, xp, resid = \
            _stacked_residual_histogram(data, positions, window_sigma, nbins, spill_period=spill_period)
        oguess = [0, max(ohist), 0, 2]
        opopt_sym, _, omask = _fit_gaussian_windowed(obins_mid, ohist, oraw_counts, oguess, window_sigma,
                                                     poisson_weight=poisson_weight)
        popt, pcov = _fit_asym_gaussian(obins_mid, ohist, oraw_counts, omask, opopt_sym, poisson_weight=poisson_weight)
    except Exception:
        return dict.fromkeys(SKEW_KEYS, np.nan)
    errs = np.sqrt(np.diag(pcov))
    return {'sigma_l': abs(popt[3]), 'sigma_l_err': errs[3], 'sigma_r': abs(popt[4]), 'sigma_r_err': errs[4]}


SKIP_COLS = ['start', 'end', 'entries', 'reason', 'runs']


def _skip_entry(chunk, reason):
    """One `skipped_reasons` entry (`SKIP_COLS` order)."""
    return (chunk['ts'].min(), chunk['ts'].max(), len(chunk), reason, sorted(int(r) for r in chunk['run'].unique()))


def build_fine_chunks(fine_base, df, target_entries, min_entries, weighted_period=True, poisson_weight=True,
                      sigma_fit_window=3, recovery_min_entries=None):
    """Split each `fine_base` row into ~`target_entries` chunks; fit period / intercept / σ / skewness per chunk.
    - A tail < `min_entries` merges into the previous chunk; never cuts inside one timestamp
    - `recovery_min_entries`: chunks between it and `min_entries` are kept if they pass the guards
    - Edge buckets (`fit_edge_buckets`) fitted after the period line; core 81, period and intercept unchanged;
      `window_offset` / `intercept` (`add_window_offset`) on top
    - σ (`fit_bunch_sigma`) and skewness (`fit_bunch_skew`): stack of every saved bucket (core + edges), entries in
      first − P/2 … last + P/2 (spill window); weak buckets (`weak_buckets`) recorded, not used
    - Returns (fine_df, skipped, skipped_reasons); `calibration_table` turns them into the calibration CSV
    - `skipped_reasons`: (start, end, entries, reason, runs) per skipped chunk (`SKIP_COLS`)
    """
    def fit_chunk(chunk):
        """One chunk's entries → (row dict, None) or (None, skip reason)."""
        bins = np.linspace(400, 2200, 1801)
        try:
            linearpopt, linearpcov, period_chi2ndof, xp_arr, peak_pos_err, background_H = \
                peak_by_peak_gaussian_fit(chunk['bunch_corrected'], bins, weighted_period=weighted_period,
                                          poisson_weight=poisson_weight)
        except Exception as e:
            return None, f'period fit: {e}'

        period = linearpopt[0]
        period_err = np.sqrt(np.diag(linearpcov)[0])

        # edge buckets: need only the period line + H; saved buckets = core + edges, all in the σ / skewness stack below
        edge_hist, edge_centers = np.histogram(chunk['bunch_corrected'], bins)[0], (bins[:-1] + bins[1:]) / 2
        edges = fit_edge_buckets(edge_hist, edge_centers, xp_arr, linearpopt, background_H, poisson_weight=poisson_weight)
        weak = weak_buckets(edge_hist, edge_centers, xp_arr, linearpopt, background_H, edges['edge_slots'])
        saved_pos = np.sort(np.concatenate([xp_arr, edges['edge_positions']]))

        try:
            sigma, sigma_err, sigma_chi2ndof = fit_bunch_sigma(chunk['bunch_corrected'], saved_pos,
                                                               window_sigma=sigma_fit_window, spill_period=period)
        except Exception as e:
            return None, f'sigma fit: {e}'
        skew = fit_bunch_skew(chunk['bunch_corrected'], saved_pos, window_sigma=sigma_fit_window,
                              poisson_weight=poisson_weight, spill_period=period)

        return {
            'run': sorted(chunk['run'].unique().tolist()),
            'start': chunk['ts'].min(),
            'end': chunk['ts'].max(),
            'entries': len(chunk),
            'fit_entries': len(chunk),
            'period': period, 'period_err': period_err,
            **line_terms(linearpopt, linearpcov),
            'sigma': sigma, 'sigma_err': sigma_err,
            'period_chi2ndof': period_chi2ndof, 'sigma_chi2ndof': sigma_chi2ndof,
            **skew,
            'peak_positions': xp_arr.tolist(), 'peak_pos_err': peak_pos_err.tolist(),
            'n_buckets': len(xp_arr) + len(edges['edge_slots']),
            **edges, **weak,
        }, None

    def too_few(n):
        return n < min_entries and (recovery_min_entries is None or n < recovery_min_entries)

    fine_rows = []
    skipped_reasons = []

    for _, row in fine_base.iterrows():
        thisdf = chunk_entries(df, row['run'], row['start'], row['end']).sort_values('ts', kind='stable')

        n_entries = len(thisdf)
        if n_entries == 0:
            continue

        edges = list(range(0, n_entries, target_entries)) + [n_entries]
        if len(edges) > 2 and (edges[-1] - edges[-2]) < min_entries:
            edges = edges[:-2] + [edges[-1]]
        # never cut inside one timestamp, so chunk_entries(start, end) recovers the chunk
        ts_sorted = thisdf['ts'].values
        edges = [edges[0]] + [e + int(np.sum(ts_sorted[e:] == ts_sorted[e - 1])) for e in edges[1:-1]] + [edges[-1]]
        edges = sorted(set(edges))
        chunks = [thisdf.iloc[i0:i1] for i0, i1 in zip(edges[:-1], edges[1:])]

        for chunk in chunks:
            if too_few(len(chunk)):
                skipped_reasons.append(_skip_entry(chunk, f'too few entries ({len(chunk)} < {min_entries})'))
                continue
            # recovery: the fit guards decide
            r, why = fit_chunk(chunk)
            if r is None:
                skipped_reasons.append(_skip_entry(chunk, why))
                continue
            fine_rows.append(r)

    fine_df = pd.DataFrame(fine_rows)
    if len(fine_df):
        fine_df = fine_df.sort_values('start', kind='stable').reset_index(drop=True)
        add_window_offset(fine_df)   # fixed slot frame + intercept over this call's chunks
    skipped_reasons.sort(key=lambda s: s[0])
    return fine_df, len(skipped_reasons), skipped_reasons


# --- Calibration CSV: list / integer columns (skipped chunks get empty lists / cells) ---
CALIB_LIST_COLS = ('run', 'peak_positions', 'peak_pos_err', 'edge_slots', 'edge_positions', 'edge_pos_err', *WEAK_KEYS)
CALIB_INT_COLS = ('fit_entries', 'n_buckets', 'window_offset')


# --- Calibration CSV columns (explicit whitelist: fixed schema) ---
# - one list for Run 1 + Run 2
CALIB_COLS = ['run', 'start', 'end', 'entries', 'fit_entries', 'period', 'period_err',
                   'intercept', 'intercept_err', 'sigma', 'sigma_err', 'period_chi2ndof',
                   'sigma_chi2ndof', *SKEW_KEYS, 'peak_positions', 'peak_pos_err',
                   'n_buckets', 'window_offset', *EDGE_KEYS, *WEAK_KEYS]


LABEL_COLS = ('segment', 'segment_status', 'constants', 'bunch_condition', 'logs')   # `labels=` below


def calibration_table(fine_df, skipped_reasons, cols, labels):
    """Fitted chunks + skipped chunks (no constants), by start → the calibration CSV rows.
    - `cols`: calibration CSV columns of the fitted chunks; adds `skip_reason` ('' unless skipped)
    - Skipped rows: `run`, `start`, `end`, `entries`, `skip_reason`; fit columns empty, list columns `[]`
    - `labels(start, end, skip_reason)` → dict of `LABEL_COLS` (`event_timestamps_run{1,2}.chunk_labels`)
    """
    fit = fine_df[cols].copy()
    fit['skip_reason'] = ''
    sk = pd.DataFrame(skipped_reasons, columns=SKIP_COLS)
    bad = pd.DataFrame({'run': sk['runs'], 'start': sk['start'], 'end': sk['end'], 'entries': sk['entries'],
                        'skip_reason': sk['reason']})
    for c in CALIB_LIST_COLS:
        if c in cols and c != 'run':
            bad[c] = [[] for _ in range(len(bad))]
    for c in CALIB_INT_COLS:
        if c in cols:
            fit[c] = fit[c].astype('Int64')   # keeps integers in the CSV next to the empty bad cells
    out = pd.concat([fit, bad], ignore_index=True)
    lab = pd.DataFrame([labels(s, e, r) for s, e, r in zip(out['start'], out['end'], out['skip_reason'])], index=out.index)
    out = pd.concat([out, lab.astype({'segment': 'Int64'})], axis=1)[[*cols, *LABEL_COLS, 'skip_reason']]
    return out.sort_values('start').reset_index(drop=True)


def chunk_entries(raw_df, runs, start, end):
    """Entries of one row (run group, fine chunk, coarse row): run in `runs`, start <= ts <= end.
    - start/end are the row's own first/last entry, so both ends are inclusive
    - The one entry rule used everywhere; returns exactly the entries a chunk was fit on (`fit_entries`)
    """
    return raw_df[raw_df['run'].isin(runs) & (raw_df['ts'] >= start) & (raw_df['ts'] <= end)]


def find_spotty_beam_gaps(raw_df, min_gap_minutes=10):
    """Mid-run entries gaps >= `min_gap_minutes` (same run both sides: beam off, DAQ on)."""
    d = raw_df.sort_values('ts').reset_index(drop=True)
    gap = d['ts'].diff()
    run_same = d['run'] == d['run'].shift(1)
    mask = (gap >= pd.Timedelta(minutes=min_gap_minutes)) & run_same
    rows = [{
        'run': int(d['run'].iloc[i]), 'gap_start': d['ts'].iloc[i - 1], 'gap_end': d['ts'].iloc[i],
        'duration_min': gap.iloc[i].total_seconds() / 60,
    } for i in d.index[mask]]
    return pd.DataFrame(rows, columns=['run', 'gap_start', 'gap_end', 'duration_min']
                        ).sort_values('gap_start').reset_index(drop=True)


# --- fixed slot frame across chunks: window offset + intercept ---

def line_terms(popt, pcov):
    """Period-line terms of one fit (popt = (period, value at fit index 0)) → `line0`, `line0_err`, `line_cov` (cov period, line0)."""
    return {'line0': popt[1], 'line0_err': np.sqrt(pcov[1, 1]), 'line_cov': pcov[0, 1]}


def unwrap_bucket0(bucket0, period):
    """Rows in time order → (phase-continuous bucket-0 value, window offset in buckets); bucket0 = unwrapped + offset · period.
    - The 81-window can start one bucket early / late (partly filled edge buckets) → the bucket-0 value jumps by one period
    - Row to row, a jump > period / 2 = window shift; smaller = real phase step (e.g. Run 1 LLRF switch 04-23, −8.5 ns)
    - Offsets are relative to the first row's window
    """
    fp, P = np.asarray(bucket0, float), np.asarray(period, float)
    off = np.zeros(len(fp), dtype=int)
    un = fp.copy()
    for i in range(1, len(fp)):
        off[i] = int(np.rint((fp[i] - un[i - 1]) / P[i]))
        un[i] = fp[i] - off[i] * P[i]
    return un, off


def add_window_offset(df, time_col='start', good=None):
    """`window_offset` + `intercept` / `intercept_err` columns, frame fixed over the rows (`unwrap_bucket0` of `line0`, `time_col` order).
    - `window_offset`: slot of fit index 0; physical slot = fit index (or edge slot) + `window_offset`; same number = same bucket in every row
    - `intercept`: period line at slot 0 of the fixed frame = `line0` − `window_offset` · period; error from the line covariance
      (= spill position from all 81 buckets: no one-period window jumps, not pulled by a half-filled first bucket)
    - `good`: rows to unwrap over (default all); others get NA
    """
    m = np.ones(len(df), dtype=bool) if good is None else np.asarray(good, dtype=bool)
    d = df[m].sort_values(time_col)
    _, off = unwrap_bucket0(d['line0'], d['period'])
    df['window_offset'] = pd.Series(off, index=d.index).reindex(df.index).astype('Int64')
    o = df['window_offset'].astype(float)
    df['intercept'] = df['line0'] - o * df['period']
    df['intercept_err'] = np.sqrt(df['line0_err'] ** 2 + o ** 2 * df['period_err'] ** 2 - 2 * o * df['line_cov'])
    return df
