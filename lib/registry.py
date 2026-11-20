"""Window / label logic shared by the run registries (`event_timestamps_run1.py`, `event_timestamps_run2.py`).

- The registries hold the data (windows, `CHUNK_NOTES`, `SEGMENT_NUMBER`, `stack_segment`) and bind these functions to it
- A chunk [start, end] is inside a window (a, b) if start >= a and end <= b + `TRANSITION_SLACK` (chunk edges = its own first / last hit)
"""

from datetime import timedelta

import numpy as np

TRANSITION_SLACK = timedelta(minutes=1)


def inside(start, end, windows, slack=TRANSITION_SLACK):
    """True if a chunk [start, end] lies inside one (a, b) of `windows`."""
    return any(start >= a and end <= b + slack for a, b in windows)


def chunk_labels(notes, segment_number, stack_segment, start, end, skip_reason=''):
    """Production label columns of one chunk [start, end] (`calib_fit_helpers.production_table(labels=)`).
    - `notes` = registry `CHUNK_NOTES` ((a, b, bunch_condition, logs), ...): a chunk inside a window takes its label
    - segment: `segment_number[stack_segment(start, end)]`, None = no constants or unsegmented
    - constants: Yes / LowStats ('too few' skip) / FitFailed
    - bunch_condition: from the note, else 'Stable' (fitted) / '-' (skipped, not investigated); logs: from the note, else ''
    """
    note = next((n for n in notes if start >= n[0] and end <= n[1] + TRANSITION_SLACK), None)
    fitted = not skip_reason
    seg = stack_segment(start, end) if fitted else None
    return {'segment': segment_number.get(seg),
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
