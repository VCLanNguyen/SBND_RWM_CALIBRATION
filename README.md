# SBND RWM Timing Calibration

Per-time-chunk calibration of SBND BNB bucket timing against the RWM: bucket positions, period, bunch width σ and skewness σR − σL. Applying it puts every entry on its bucket: `calib_time_ns` = entry time − bucket position.

Method and results: SBND timing tech note (SBN-doc-XXXXX, TODO), Sec. "RWM Calibration".

## Layout

```
calibration/
  calibration.ipynb         ntuple → calibration CSV (one per run)
  overview.ipynb            σ, σR − σL, period, intercept per time chunk vs time (tech note Fig. 22)
  segmentation.ipynb        per-segment σ / σR − σL and stacked bunch (tech note Figs. 23-24, Tab. 2)
  apply_calibration.ipynb   how to apply the calibration CSV to an uncalibrated ntuple, step by step
lib/
  calib_fit_helpers.py      ntuple loader, chunking, bucket finding, period / σ / skewness fits, CSV columns
  run_table.py              run table, run groups and splits
  event_timestamps_run1.py  Run 1 registry: hard boundaries, bad runs, segments, chunk labels (+ figure overlays)
  event_timestamps_run2.py  Run 2 registry (same layout)
  registry.py               window / label logic shared by the two registries
  stacking.py               figures: CSV loader, entry → bucket assignment, stacked-bunch fit
  plot_helpers.py           figures: style, segment colours, per-chunk series
  paths.py                  input / output paths (`set_paths`), CSV writer
```

## Inputs

| Input | What | Set in |
|---|---|---|
| Run 1 ntuple | `TimingCalibAna` module v2, tree `timingcalib/timingCalibTree` (`run1_timing_calib_ntuple.root`, 1.7 GB) | `NTUPLE` in the Setup cell |
| Run 2 ntuple | same, Run 2 (`run2_timing_calib_ntuple.root`, 1.3 GB) | `NTUPLE` in the Setup cell |

- Any readable path works (local or Fermilab); the input files are never written to

## Data on Fermilab

`/pnfs/sbnd/persistent/users/lnguyen/RWM_Calibration/`:

| File | What | md5 |
|---|---|---|
| `run1_timing_calib_ntuple.root` | Run 1 input ntuple (1.7 GB) | `ee97140b68395c09360ca5ece60dff20` |
| `run2_timing_calib_ntuple.root` | Run 2 input ntuple (1.3 GB) | `85b9ce86b14592f0d08889a36c015225` |
| `run1_calibration.csv` | Run 1 calibration CSV = `calibration.ipynb` output `run1/sbnd_fine_calibration.csv` from these ntuples | `276e6ea68ec1bda1ca81e8fde97d9532` |
| `run2_calibration.csv` | Run 2 calibration CSV (same, Run 2) | `cc156df89a5d95069efcc09c667eadae` |

- The calibration CSVs are the ones used in the tech note; to apply them directly, point `CALIB_CSV` in `apply_calibration.ipynb` at them
- How the ntuples were made (module, release, input datasets): `calibration/calibration.ipynb`, Setup › "Where the Input Ntuples Come From"

## Environment

Python 3.10: `pip install -r requirements.txt`

## Run

First set the paths at the top of each notebook's Setup cell ("Set these paths"):

| Notebook | Paths |
|---|---|
| `calibration`, `segmentation`, `apply_calibration` | `NTUPLE` = your copy of the two ntuples; `CALIB_DIR` = output folder |
| `overview` | `CALIB_DIR` |

- `NTUPLE` is a placeholder (`/path/to/…`) and must be set; `CALIB_DIR` defaults to `<repo>/output` and must be the same in all four
- Then:

```bash
cd calibration
for nb in calibration overview segmentation apply_calibration; do
    jupyter nbconvert --to notebook --execute --inplace $nb.ipynb
done
```

| Notebook | Time | Peak RAM | Writes |
|---|---|---|---|
| `calibration` | ~2.5 min (first run: builds the entry caches) | 2.2 GB | `CALIB_DIR/run{1,2}/sbnd_fine_calibration.csv` + `*_lean.pkl` caches |
| `overview` | ~10 s | 0.3 GB | `CALIB_DIR/figures/calib_results_overview.png` |
| `segmentation` | ~15 s | 1.2 GB | `CALIB_DIR/figures/calib_results_segments_run{1,2}.png` |
| `apply_calibration` | ~15 s | 1.6 GB | nothing (table in memory) |

Run order: `calibration` first; the other three read its CSVs (`segmentation` and `apply_calibration` also read the ntuples).

## Output: `sbnd_fine_calibration.csv` (one per run)

One row per time chunk, in time order; uncalibrated time chunks are kept with empty constants.

| Column | Meaning |
|---|---|
| `run`, `start`, `end`, `entries` | DAQ runs, first / last entry time (UTC, ns), entries of the time chunk |
| `period`, `intercept` (+ `_err`) | bucket spacing [ns]; position of bucket 0 of the fixed bucket frame [ns] |
| `sigma`, `sigma_l`, `sigma_r` (+ `_err`) | bunch width σ; σL / σR of the asymmetric fit (skewness σR − σL) [ns] |
| `peak_positions`, `peak_pos_err` | positions of the 81 core buckets [ns] |
| `edge_slots`, `edge_positions`, `edge_pos_err` | partly filled edge buckets (≤ 2 per side): slot and position |
| `window_offset` | bucket number = fit index (or edge slot) + `window_offset` (same number = same RF bucket in every time chunk) |
| `weak_slots`, `weak_positions`, `weak_pvalue` | next slot past the outermost bucket on each side: recorded, not used |
| `period_chi2ndof`, `sigma_chi2ndof`, `n_buckets`, `fit_entries` | fit quality, buckets used, entries fitted |
| `segment`, `segment_status`, `constants`, `bunch_condition`, `logs`, `skip_reason` | stacking segment (1, 2, …), its status (Stable / Unstable), Yes / LowStats / FitFailed, beam condition and supporting log, why skipped |

How to apply it to an uncalibrated ntuple: `apply_calibration.ipynb`.

## Changing the calibration

- Time-chunk boundaries, bad runs, segments, chunk labels: `lib/event_timestamps_run{1,2}.py` (`FINE_SPLITS`, `BAD_RUNS`, `SEGMENTS`, `CHUNK_NOTES`), then rerun `calibration.ipynb`
- Chunk size and fit settings: Setup cell of `calibration.ipynb` and the constants at the top of `lib/calib_fit_helpers.py`

## TODO

- **Incomplete production:** The data processed for this tech note had failed jobs, so about 2 % of entries are missing.
- **Production Run 2A:** Data from 31-Oct-2025 to 24-Jan-2026 (runs 19683–19854) are not processed yet. They need to be calibrated and segmented once processed.
- **Next production release:** This calibration is for release `v10_14_02_10`. Remake it for the next production release.
