"""Run table and row-splitting helpers.

- `build_run_table(df)`: one row per run (start/end/duration/entries) from the raw ntuple
- Bad runs: the registries (`event_timestamps_run{1,2}.BAD_RUNS`)
- `split_row_by_run`: split a run-group row at a timestamp (registry `FINE_SPLITS`)
- Row convention (every table here, and every chunk built from them):
  - `start`/`end` = the row's own first/last entry (exact ns `ts`)
  - the row's entries = run in `run` list and start <= ts <= end (`calib_fit_helpers.chunk_hits`)
"""

import pandas as pd


def build_run_table(df):
    """Per-run start/end (first/last entry `ts`, exact ns), duration, entries."""
    g = df.groupby('run')['ts']
    out = pd.DataFrame({'start': g.min(), 'end': g.max(), 'entries': g.size()}).reset_index()
    out['duration'] = out['end'] - out['start']
    return out[['run', 'start', 'end', 'duration', 'entries']]


# --- Run grouping ---
def grab_combined_run(df, combined_run):
    thisdf = pd.DataFrame()
    for run in combined_run:
        maskeddf = df[df['run']==run]
        thisdf = pd.concat([thisdf, maskeddf])

    start = thisdf['start'].min()
    end = thisdf['end'].max()
    duration = end - start
    entries = sum(thisdf['entries'])

    newrow = pd.DataFrame({'run':[combined_run], 'start':[start], 'end':[end], 'duration':[duration], 'entries':[entries]})
    return newrow


# --- Row splitting ---
def split_row_by_run(df, ts_df, target_run, split_time, ts_col='ts'):
    """Split the row of `df` whose 'run' list contains `target_run` at `split_time`.

    - `target_run`: a single run number or list/tuple -- matches any row whose
      'run' list overlaps with it; if multiple rows match, picks the one whose
      [start, end] actually contains `split_time`
    - Halves: entries with ts < split_time / ts >= split_time (row's runs only); each half's
      start/end snap to its own first/last entry, entries counted from `ts_df`
    - No-op if `split_time` already falls between two rows of the run, or leaves one half empty
      (incl. < 1 s before the row start: a split truncated to seconds on a run's first entry)
    - `split_time`: anything `pd.Timestamp()` accepts
    """
    split_time = pd.Timestamp(split_time)
    target_runs = set(target_run) if isinstance(target_run, (list, tuple, set)) else {target_run}

    mask = df['run'].apply(lambda runs: bool(target_runs & set(runs)))
    if mask.sum() == 0:
        print(f"run {target_run} not found, skipping")
        return df

    candidates = df.index[mask]
    if len(candidates) > 1:
        fits = [i for i in candidates if df.loc[i, 'start'] <= split_time <= df.loc[i, 'end']]
        if len(fits) == 0:
            if df.loc[candidates, 'start'].min() < split_time < df.loc[candidates, 'end'].max():
                return df
            raise ValueError(
                f"run {target_run} found in {len(candidates)} rows, "
                f"but split_time {split_time} falls inside none of them"
            )
        if len(fits) > 1:
            raise ValueError(
                f"run {target_run} found in {len(candidates)} rows, "
                f"and split_time {split_time} falls inside more than one ({fits})"
            )
        idx = fits[0]
    else:
        idx = candidates[0]

    row = df.loc[idx]

    if pd.Timedelta(0) < row['start'] - split_time < pd.Timedelta(seconds=1):
        return df   # split truncated to seconds, row starts at its first entry: already a boundary
    if not (row['start'] <= split_time <= row['end']):
        raise ValueError(f"split_time {split_time} is outside row's [start, end] range")

    ts = ts_df[ts_col]
    ts = ts[ts_df['run'].isin(row['run']) & (ts >= row['start']) & (ts <= row['end'])]
    ts_a, ts_b = ts[ts < split_time], ts[ts >= split_time]
    if len(ts_a) == 0 or len(ts_b) == 0:
        return df

    row_a = row.copy()
    row_a['end'] = ts_a.max()
    row_a['entries'] = len(ts_a)
    row_a['duration'] = row_a['end'] - row_a['start']

    row_b = row.copy()
    row_b['start'] = ts_b.min()
    row_b['entries'] = len(ts_b)
    row_b['duration'] = row_b['end'] - row_b['start']

    return pd.concat(
        [df.iloc[:idx], pd.DataFrame([row_a, row_b]), df.iloc[idx + 1:]],
        ignore_index=True,
    )

# --- Fine-calibration base rows (calibration notebook) ---
def consecutive_groups(runs):
    """Consecutive run numbers -> groups; a missing number breaks the group."""
    runs = sorted(runs)
    groups = [[runs[0]]]
    for r in runs[1:]:
        if r == groups[-1][-1] + 1:
            groups[-1].append(r)
        else:
            groups.append([r])
    return groups


def build_fine_base(df, dfrun, groups, splits, spotty_min_gap=45, spotty_skip_s=120, counts=False):
    """Fine-calibration base rows: run groups → split at every (run, time) of `splits`, in order → split at every
    mid-run spotty-beam gap >= `spotty_min_gap` min (gap end; skipped within `spotty_skip_s` of a named split).
    - `df`: raw entries (spotty gaps + split entries); `dfrun`: `build_run_table(df)`; `groups`: run-number groups
    - `counts=True` → (base, spotty splits applied, spotty gaps found)
    """
    from calib_fit_helpers import find_spotty_beam_gaps
    base = pd.concat([grab_combined_run(dfrun, g) for g in groups], ignore_index=True)
    base = base.sort_values('start').reset_index(drop=True)
    for run, t in splits:
        base = split_row_by_run(base, df, run, t, ts_col='ts')
    named = [pd.Timestamp(t) for _, t in splits]
    gaps = find_spotty_beam_gaps(df, min_gap_minutes=spotty_min_gap)
    n_applied = 0
    for _, g in gaps.iterrows():
        if any(abs((g['gap_end'] - t).total_seconds()) < spotty_skip_s for t in named):
            continue
        base = split_row_by_run(base, df, g['run'], g['gap_end'], ts_col='ts')
        n_applied += 1
    return (base, n_applied, len(gaps)) if counts else base
