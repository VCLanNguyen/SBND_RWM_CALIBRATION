# SBND RWM Timing Calibration

Per-chunk calibration of SBND BNB bucket timing against the RWM: bucket centres, period, σ and σR − σL. Applying it puts every neutrino slice on its bucket: `calib_time_ns` = slice time − bucket centre.

Method and results: SBND timing technote (SBN-doc-XXXXX, TODO), Sec. "RWM Calibration".

## Layout

```
calibration/
  calibration.ipynb         ntuple → calibration table (one CSV per run)
  overview.ipynb            σ, σR − σL, period, intercept per chunk vs time (technote Fig. 22)
  segmentation.ipynb        per-segment σ / σR − σL and stacked bunch (technote Figs. 23-24, Tab. calib_segments)
  apply_calibration.ipynb   how to apply the table to the slices of an ntuple, step by step
lib/
  calib_fit_helpers.py      ntuple loader, chunking, bucket finding, period / σ / skewness fits, CSV columns
  run_table.py              run table, run groups and splits, DAQ run-log reader
  event_timestamps_run1.py  Run 1 registry: hard boundaries, bad runs, segments, chunk labels (+ figure overlays)
  event_timestamps_run2.py  Run 2 registry (same layout)
  registry.py               window / label logic shared by the two registries
  stacking.py               figures: CSV loader, slice → bucket assignment, stacked-bunch fit
  plot_helpers.py           figures: style, segment colours, per-chunk series
  paths.py                  input / output paths (`set_paths`), CSV writer
```

## Inputs

| Input | What | Variable |
|---|---|---|
| Run 1 ntuple | `TimingCalibAna` module v2, tree `timingcalib/timingCalibTree` (`run1_timing_calib_ntuple.root`, 1.7 GB) | `$SBND_RWM_NTUPLE_RUN1` |
| Run 2 ntuple | same, Run 2 (`run2_timing_calib_ntuple.root`, 1.3 GB) | `$SBND_RWM_NTUPLE_RUN2` |
| DAQ run log | `DAQ Runs - DAQ Runs (Run 2).csv`; figures only (time span of Run 2 segment A) | `$SBND_RWM_DAQ_RUN2` |

- Any readable path works (local or Fermilab); the input files are never written to

## Data on Fermilab

`/pnfs/sbnd/persistent/users/lnguyen/RWM_Calibration/`:

| File | What | md5 |
|---|---|---|
| `run1_timing_calib_ntuple.root` | Run 1 input ntuple (1.7 GB) | `ee97140b68395c09360ca5ece60dff20` |
| `run2_timing_calib_ntuple.root` | Run 2 input ntuple (1.3 GB) | `85b9ce86b14592f0d08889a36c015225` |
| `run1_calibration.csv` | Run 1 calibration table = `calibration.ipynb` output `run1/sbnd_fine_calibration.csv` from these ntuples | `276e6ea68ec1bda1ca81e8fde97d9532` |
| `run2_calibration.csv` | Run 2 calibration table (same, Run 2) | `cc156df89a5d95069efcc09c667eadae` |

- The calibration tables are the ones used in the technote; to apply them directly, point `CALIB_CSV` in `apply_calibration.ipynb` at them
- How the ntuples were made (module, release, input datasets): `calibration/calibration.ipynb`, Setup › "Where the Input Ntuples Come From"
- Each notebook's Setup cell reads these variables; the paths can also be typed there

## Environment

Python 3.10 with numpy 2.2, pandas 2.3, scipy 1.15, uproot 5.7, matplotlib 3.10, jinja2 3.1 (pandas table styling), jupyter / nbconvert.

## Run

```bash
export SBND_RWM_NTUPLE_RUN1=/pnfs/sbnd/persistent/users/lnguyen/RWM_Calibration/run1_timing_calib_ntuple.root
export SBND_RWM_NTUPLE_RUN2=/pnfs/sbnd/persistent/users/lnguyen/RWM_Calibration/run2_timing_calib_ntuple.root
export SBND_RWM_DAQ_RUN2="/path/to/DAQ Runs - DAQ Runs (Run 2).csv"
export SBND_RWM_OUT=$PWD/output          # default: <repo>/output

cd calibration
for nb in calibration overview segmentation apply_calibration; do
    jupyter nbconvert --to notebook --execute --inplace $nb.ipynb
done
```

| Notebook | Time | Peak RAM | Writes |
|---|---|---|---|
| `calibration` | ~2.5 min (first run: builds the slice caches) | 2.2 GB | `$SBND_RWM_OUT/run{1,2}/sbnd_fine_calibration.csv` + `*_lean.pkl` caches |
| `overview` | ~10 s | 0.3 GB | `$SBND_RWM_OUT/figures/calib_results_overview.png` |
| `segmentation` | ~15 s | 1.2 GB | `$SBND_RWM_OUT/figures/calib_results_segments_run{1,2}.png` |
| `apply_calibration` | ~15 s | 1.6 GB | nothing (table in memory) |

Run order: `calibration` first; the other three read its CSVs.

## Output: `sbnd_fine_calibration.csv` (one per run)

One row per time chunk, in time order; skipped chunks are kept with empty constants.

| Column | Meaning |
|---|---|
| `run`, `start`, `end`, `entries` | DAQ runs, first / last slice time (UTC, ns), slices of the chunk |
| `period`, `intercept` (+ `_err`) | bucket spacing [ns]; centre of bucket 0 of the fixed bucket frame [ns] |
| `sigma`, `sigma_l`, `sigma_r` (+ `_err`) | bunch width; σL / σR of an asymmetric fit (skewness = σR − σL) [ns] |
| `peak_positions`, `peak_pos_err` | centres of the 81 core buckets [ns] |
| `edge_slots`, `edge_positions`, `edge_pos_err` | partly filled edge buckets (≤ 2 per side): slot and centre |
| `window_offset` | bucket number = fit index (or edge slot) + `window_offset` (same number = same RF bucket in every chunk) |
| `weak_slots`, `weak_positions`, `weak_pvalue` | next slot past the outermost bucket on each side: recorded, not used |
| `period_chi2ndof`, `sigma_chi2ndof`, `n_buckets`, `fit_entries` | fit quality, buckets used, slices fitted |
| `segment`, `segment_status`, `constants`, `bunch_condition`, `logs`, `skip_reason` | stacking segment (1, 2, …), its status (Stable / Unstable), Yes / LowStats / FitFailed, beam condition and supporting log, why skipped |

How to use it on slices: `apply_calibration.ipynb`.

## Changing the calibration

- Hard chunk boundaries, bad runs, segments, chunk labels: `lib/event_timestamps_run{1,2}.py` (`FINE_SPLITS`, `BAD_RUNS`, `SEGMENTS`, `CHUNK_NOTES`), then rerun `calibration.ipynb`
- Chunk size and fit settings: Setup cell of `calibration.ipynb` and the constants at the top of `lib/calib_fit_helpers.py`

## TODO

- **Incomplete production:** The data processed for this technote had failed jobs, so about 2 % of slices are missing.
- **Run 2 segment A:** Runs 19683–19854 at the start of Run 2 are not processed yet. They need to be calibrated and segmented once processed.
- **Timing reconstruction changes:** Remake the calibration after any change to the timing reconstruction modules (technote Sec. "Timing Reconstruction").
