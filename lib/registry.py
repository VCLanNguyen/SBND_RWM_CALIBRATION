"""Window / label logic shared by the run registries (`event_timestamps_run1.py`, `event_timestamps_run2.py`).

- The registries hold the data (windows, `CHUNK_NOTES`, `SEGMENTS`) and bind these functions to it
- `SEGMENTS` rows = (number, name, windows, status, beam state, evidence); number None = own beam state, not segmented;
  a chunk takes the innermost row with a window that contains it (`segment_row`)
- A chunk [start, end] is inside a window (a, b) if start >= a and end <= b + `TRANSITION_SLACK` (chunk edges = its own first / last entry)
"""

from datetime import timedelta

import numpy as np

TRANSITION_SLACK = timedelta(minutes=1)


def inside(start, end, windows, slack=TRANSITION_SLACK):
    """True if a chunk [start, end] lies inside one (a, b) of `windows`."""
    return any(start >= a and end <= b + slack for a, b in windows)


def segment_row(segments, start, end, slack=TRANSITION_SLACK):
    """Row of a registry `SEGMENTS` table for a chunk [start, end]; None if no window contains it (chunk straddles a boundary).
    - Window (a, b) contains the chunk if start >= a and end <= b + `slack`; a / b None = open (run start / end)
    - Boundaries = first entry after the change (as `FINE_SPLITS`), so every boundary sits on a chunk edge
    - Nested windows (narrow-bunch periods, own beam states): the innermost (shortest) window wins
    """
    best, best_len = None, None
    for row in segments:
        for a, b in row[2]:
            if (a is None or start >= a) and (b is None or end <= b + slack):
                n = b - a if a is not None and b is not None else timedelta.max
                if best is None or n < best_len:
                    best, best_len = row, n
    return best


def stack_segment(segments, start, end):
    """Segment name of a chunk [start, end]; None if not segmented (own beam state) or straddling a boundary."""
    row = segment_row(segments, start, end)
    return row[1] if row is not None and row[0] is not None else None


def is_unsegmented(segments, start, end):
    """True if a chunk [start, end] lies in an own beam state (`SEGMENTS` row with number None)."""
    row = segment_row(segments, start, end)
    return row is not None and row[0] is None


def chunk_labels(notes, segments, start, end, skip_reason=''):
    """Calibration CSV label columns of one chunk [start, end] (`calib_fit_helpers.production_table(labels=)`).
    - `notes` = registry `CHUNK_NOTES` ((a, b, bunch_condition, logs), ...): a chunk inside a window takes its label
    - segment: number of its `segments` row (`segment_row`), None = no constants or not segmented
    - segment_status: Stable / Unstable from the `segments` row, '' when no segment
    - constants: Yes / LowStats ('too few' skip) / FitFailed
    - bunch_condition: from the note, else 'Stable' (fitted) / '-' (skipped, not investigated); logs: from the note, else ''
    """
    note = next((n for n in notes if start >= n[0] and end <= n[1] + TRANSITION_SLACK), None)
    fitted = not skip_reason
    row = segment_row(segments, start, end) if fitted else None
    seg = row if row is not None and row[0] is not None else None
    return {'segment': seg[0] if seg else None,
            'segment_status': seg[3] if seg else '',
            'constants': 'Yes' if fitted else 'LowStats' if skip_reason.startswith('too few') else 'FitFailed',
            'bunch_condition': note[2] if note else ('Stable' if fitted else '-'),
            'logs': note[3] if note else ''}


def mask(df, test):
    """Boolean array over df rows (needs 'start' / 'end'): `test(start, end)` per row."""
    return np.array([test(s, e) for s, e in zip(df['start'], df['end'])], dtype=bool)


def assign_segments(db, ev):
    """Segment columns on a table of fitted chunks `db`, in place; `ev` = run registry module (figures).
    - `seg`: `ev.stack_segment` name, None if not segmented; `unsegmented`: own beam state (`ev.NO_SEGMENT_WINDOWS`)
    - Asserts no other chunk straddles a segment boundary
    """
    db['seg'] = [ev.stack_segment(s, e) for s, e in zip(db['start'], db['end'])]
    db['unsegmented'] = mask(db, ev.is_unsegmented)
    assert (db['seg'].notna() | db['unsegmented']).all(), 'chunk straddles a segment boundary'
    return db
