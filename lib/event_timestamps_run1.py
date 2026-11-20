"""Run 1 registry: beam-condition timestamps, hard chunk boundaries, bad runs, stacking segments and chunk labels.

- Used by `calibration/calibration.ipynb` (`ev1`): `FINE_SPLITS`, `BAD_RUNS`, `chunk_labels`; Run 2: `event_timestamps_run2.py`
- Timestamps in chronological order, numbered t05, t06, ...; letter suffixes (t07a) keep the order; `[SPLIT]` = listed in `FINE_SPLITS`
- Split time = first entry after the change, truncated to whole seconds (`ts < split` goes before)
- All times UTC; ELOG / DAQ-log / ECL beam-checklist times (US Central) converted
- Sources: SBND ntuple (bunch position / width per chunk), AD RWM traces, AD ELOG, SBND ECL beam checklist, DAQ run log
"""

from datetime import datetime, timedelta

import registry as _reg

t05_downtime_mar_start = datetime(2025, 3, 7, 5, 7, 22)       # beam downtime, KRF3 klystron replaced
t05_downtime_mar_first_entry = datetime(2025, 3, 8, 4, 35, 24)   # first entry after the downtime (run 18301, low-intensity linac beam); before_a / before_b segment boundary
t06_downtime_mar20_start = datetime(2025, 3, 20, 9, 58, 59)   # scheduled BNB shutdown (run 18334 end)
t06_downtime_mar20_end = datetime(2025, 3, 21, 1, 41, 56)     # run 18340 start = first entry after the shutdown; before_b / before_c segment boundary

t06a_beam_return_min = datetime(2025, 3, 21, 1, 41, 56)       # beam back after the shutdown, "expect higher intensity" (checklist); not segmented
t06a_beam_return_max = datetime(2025, 3, 22, 15, 15, 15)      # first hit after the 14:10 → 15:15 beam gap

t07_rotation_ramp_end = datetime(2025, 4, 3, 20, 0, 0)    # [SPLIT] run 18388 -- bunch rotation commissioning: width step (AD σ 1.83-2.11 → 2.06-2.34); before → during
t07a_rotation_settled = datetime(2025, 4, 7, 21, 1, 34)   # [SPLIT] run 18394 -- rotation in operation, width 2.32 → 1.93 ns; during → after

t14_llrf1 = datetime(2025, 4, 23, 16, 40, 4)              # [SPLIT] run 18443 -- LLRF1 analog → digital, inside the 15:10-16:40 beam gap
t14a_llrf1_settled = datetime(2025, 4, 23, 19, 0, 0)      # [SPLIT] run 18443 -- LLRF1 settled (width back to ~2.4 ns, phase stable)

t15_llrf2_prelude_min = datetime(2025, 4, 29, 14, 0, 0)   # [SPLIT] run 18445 -- pre-switch disturbance: rate halves, width ~3.0 ns, phase +0.8 ns
t15_llrf2_prelude_max = datetime(2025, 4, 29, 15, 10, 0)  # [SPLIT] run 18445 -- disturbance over, back to the pre-switch state
t15a_spottybeam1_min = datetime(2025, 4, 29, 17, 0, 12)   # run 18445, ~50 min beam gap starts
t15a_spottybeam1_max = datetime(2025, 4, 29, 17, 50, 27)  # [SPLIT] run 18445 -- LLRF2 intermediate state starts (phase −8.2 → −4.9 ns)
t16_llrf2 = datetime(2025, 4, 29, 21, 3, 48)              # [SPLIT] run 18445 -- LLRF2 back to analog, at the end of a ~44 min beam gap

t18_llrf3_switch = datetime(2025, 5, 6, 14, 40, 29)       # [SPLIT] run 18466 -- LLRF3 analog → digital: half-bucket jump (+9.2 ns)
t18_llrf3_flip_back = datetime(2025, 5, 6, 18, 54, 0)     # [SPLIT] run 18466 -- back to the pre-switch lock point (−1.5 ns)
t18_llrf3_settled = datetime(2025, 5, 6, 20, 26, 3)       # [SPLIT] run 18466 -- digital again (= following chunks)

t19a_may22_phase_flip_min = datetime(2025, 5, 22, 14, 35, 27)  # [SPLIT] run 18488 -- undeclared half-bucket (~9.5 ns) phase flip
t19a_may22_phase_flip_max = datetime(2025, 5, 22, 17, 27, 13)  # [SPLIT] run 18488 -- back and settled
t19b_may23_narrow_start = datetime(2025, 5, 23, 16, 59, 47)    # run 18488 -- first hit after a ~2 h beam gap = narrow state starts; not segmented
t19b_may23_narrow_end = datetime(2025, 5, 23, 20, 58, 0)       # [SPLIT] run 18488 -- narrow stretch ends (σ 1.49 vs 1.77 after)

t21_llrf4 = datetime(2025, 6, 2, 14, 29, 57)              # [SPLIT] run 18504 -- LLRF4 → digital, inside the 14:14-14:29 beam gap
t21_llrf4_settled = datetime(2025, 6, 2, 17, 20, 0)       # [SPLIT] run 18504 -- post-switch offset (+0.7 ns) / width (~3 ns) settled
t21a_jun09_offset_step = datetime(2025, 6, 9, 15, 51, 1)  # [SPLIT] run 18507 -- −0.7 ns offset step, cause unknown

t22_bnb_study_min = datetime(2025, 6, 16, 14, 5, 1)       # [SPLIT] run 18518 start -- BNB bunch-width study (narrower); not segmented
t22_bnb_study_max = datetime(2025, 6, 16, 21, 6, 17)      # run 18518 end = run-group edge

t23a_waterpump_onset = datetime(2025, 6, 26, 12, 44, 11)  # [SPLIT] run 18543 -- −1.7 ns offset onset at the Booster study start (ECL)
t24_booster_start = datetime(2025, 6, 26, 12, 31, 0)      # Booster study: beam stopped ~12:30, back 17:48 (ECL)
t24_booster_end = datetime(2025, 6, 26, 17, 48, 0)
t25a_waterpump_step2 = datetime(2025, 6, 27, 15, 52, 22)  # [SPLIT] run 18551 -- +0.55 ns offset step, inside an ECL "intermittent BNB, tests" window
t26a_jun30_widening_start = datetime(2025, 6, 30, 16, 0, 0)   # [SPLIT] run 18553 -- wider beam from 16:00 (no beam gap)
t28a_waterpump_end = datetime(2025, 6, 30, 19, 47, 20)    # [SPLIT] run 18553 -- +0.6 ns step after BRF Station 8 work (ECL), ~1 h ramp
t28b_waterpump_ramp_end = datetime(2025, 6, 30, 21, 0, 0)  # [SPLIT] run 18553 -- offset ramp done (−1.7 → −1.2 ns)

t29_flatinjection_min = datetime(2025, 7, 2, 14, 14, 49)  # [SPLIT] run 18576 -- flat-injection study, filling +1 bucket
t29_flatinjection_max = datetime(2025, 7, 3, 0, 35, 47)   # [SPLIT] run 18577 -- study ends


# --- Hard chunk boundaries: (run, split time), chronological ---
# - Every beam transition (ramp, settling, disturbance) is bracketed on both sides, so it sits in its own chunk
FINE_SPLITS = (
    (18388, t07_rotation_ramp_end),
    (18394, t07a_rotation_settled),
    (18443, t14_llrf1),
    (18443, t14a_llrf1_settled),
    (18445, t15_llrf2_prelude_min),
    (18445, t15_llrf2_prelude_max),
    (18445, t15a_spottybeam1_max),
    (18445, t16_llrf2),
    (18466, t18_llrf3_switch),
    (18466, t18_llrf3_flip_back),
    (18466, t18_llrf3_settled),
    (18488, t19a_may22_phase_flip_min),
    (18488, t19a_may22_phase_flip_max),
    (18488, t19b_may23_narrow_end),
    (18504, t21_llrf4),
    (18504, t21_llrf4_settled),
    (18507, t21a_jun09_offset_step),
    (18518, t22_bnb_study_min),
    (18543, t23a_waterpump_onset),
    (18551, t25a_waterpump_step2),
    (18553, t26a_jun30_widening_start),
    (18553, t28a_waterpump_end),
    (18553, t28b_waterpump_ramp_end),
    (18576, t29_flatinjection_min),
    (18577, t29_flatinjection_max),
)

# --- Bad runs (dropped before chunking) ---
BAD_RUNS = (
    18266,   # 02-24 linac trip, spotty run
    18452,
    18457,
    18549,
)

# --- Stacking segments: before_a / before_b / before_c / during / after (TN Sec. 5.4, App. B) ---
# - One row per segment, like `FINE_SPLITS`: (number, name, windows, status, beam state, evidence); number None = own beam state, not segmented
# - Boundary = first entry after the change (on a chunk edge); a chunk takes the innermost row whose window contains it (`registry.segment_row`)
# - 'before' cut at two beam downtimes (chunk-σ stability scan over the ≥ 1 h SBND gaps)
#   - before_b: KRF3 klystron recovery, reduced and stepping intensity, σ follows it → Unstable (no b / response / slice-correction product)
#   - before_c starts at the end of the `beam_return` window (2 chunks, not segmented): first entry 22-Mar-2025 15:15
# - Own beam state with only a few chunks → not segmented
#   - bnb_study: Jun 16 BNB bunch-width study, run 18518 = one chunk, σ 1.59 vs after-rotation mean 1.77 (std 0.05)
#   - may23_narrow: 05-23 16:59 → 20:57 = one chunk, narrow throughout (σ 1.53 vs neighbours 1.79 / 1.78), flat inside
SEGMENTS = (
    (1, 'before_a', ((None, t05_downtime_mar_first_entry),), 'Stable', 'Before rotation',
     'Run start to the 07-Mar-2025 beam downtime (KRF3 klystron replaced, checklist)'),
    (2, 'before_b', ((t05_downtime_mar_first_entry, t06_downtime_mar20_end),), 'Unstable', 'Before rotation',
     'Between the 07-Mar-2025 and 20-Mar-2025 beam downtimes; KRF3 klystron recovery, bunch width follows the stepping intensity'),
    (3, 'before_c', ((t06a_beam_return_max, t07_rotation_ramp_end),), 'Stable', 'Before rotation',
     '20-Mar-2025 scheduled BNB shutdown to the start of the bunch rotation commissioning'),
    (4, 'during', ((t07_rotation_ramp_end, t07a_rotation_settled),), 'Stable', 'Rotation commissioning',
     'Bunch rotation commissioning: width step at the start, width drop at the end (raw RWM traces, AD ELOG)'),
    (5, 'after', ((t07a_rotation_settled, None),), 'Stable', 'After commissioning',
     'End of the bunch rotation commissioning to the run end'),
    (None, 'beam_return', ((t06a_beam_return_min, t06a_beam_return_max),), '', 'Short beam state',
     'Beam return after the 20-Mar-2025 shutdown, "expect higher intensity" (checklist); 2 chunks'),
    (None, 'may23_narrow', ((t19b_may23_narrow_start, t19b_may23_narrow_end),), '', 'Short beam state',
     'Narrow bunches after a ~1 h BNB downtime (checklist); 1 chunk'),
    (None, 'bnb_study', ((t22_bnb_study_min, t22_bnb_study_max),), '', 'Short beam state',
     'BNB bunch-width study (DAQ log, AD ELOG); 1 chunk'),
)
SEGMENT_NUMBER = {r[1]: r[0] for r in SEGMENTS if r[0] is not None}
STACK_SEGMENTS = tuple(sorted(SEGMENT_NUMBER, key=SEGMENT_NUMBER.get))
UNSTABLE_SEGMENTS = tuple(r[1] for r in SEGMENTS if r[3] == 'Unstable')   # drawn black, not stacked
NO_SEGMENT_WINDOWS = {r[1]: r[2][0] for r in SEGMENTS if r[0] is None}
NO_SEGMENT_LABEL = 'Own beam state, not segmented'

def is_unsegmented(start, end):
    """True if a chunk [start, end] lies in an own beam state (`SEGMENTS` row with number None)."""
    return _reg.is_unsegmented(SEGMENTS, start, end)

def stack_segment(start, end):
    """Stacking segment of a chunk [start, end] (`registry.segment_row`); None if not segmented or straddling a boundary."""
    return _reg.stack_segment(SEGMENTS, start, end)


# --- Production label columns ---
# - segment: 1-5 = before_a / before_b / before_c / during / after (`SEGMENTS`), None = no constants or not segmented
# - segment_status: Stable / Unstable (`SEGMENTS`), '' when no segment
# - constants: Yes / LowStats / FitFailed (from `skip_reason`)
# - bunch_condition: from `CHUNK_NOTES` where investigated AND a log supports it; fitted chunks otherwise 'Stable', skipped '-'
# - logs: one-line summary of the supporting log (checklist = ECL beam checklist, AD ELOG, DAQ log); '' if none
# - a chunk takes the note of the window it lies inside (+ `registry.TRANSITION_SLACK`); windows don't overlap
CHUNK_NOTES = (
    (datetime(2025, 3, 8, 4, 35, 24), datetime(2025, 3, 8, 6, 32, 50), '-', 'Linac issue, low-intensity beam overnight (checklist)'),
    (datetime(2025, 3, 10, 18, 11, 30), datetime(2025, 3, 10, 19, 33, 30), '-', 'KRF3 klystron work, BNB down ~30 min (checklist)'),
    (t14_llrf1, t14a_llrf1_settled, 'Width Changing', 'LLRF1 analog to digital, "new low-level RF" (checklist)'),
    (t15_llrf2_prelude_min, t15_llrf2_prelude_max, 'Width + Step Changing', 'Reduced intensity "due to an issue" before LLRF2 (checklist)'),
    (t15_llrf2_prelude_max, t15a_spottybeam1_min, '-', 'Beam back 15:31, off 17:08 (checklist)'),
    (t15a_spottybeam1_max, t16_llrf2, 'Step Changing', 'LLRF2 intermediate state (LLRF test #2)'),
    (t18_llrf3_flip_back, t18_llrf3_settled, 'Step Changing', 'LLRF3 test, back to the analog lock point (AD ELOG)'),
    (t19a_may22_phase_flip_min, t19a_may22_phase_flip_max, 'Width + Step Changing', 'Half-bucket phase flip during accelerator work (checklist)'),
    (t19b_may23_narrow_start, t19b_may23_narrow_end, 'Stable', 'Narrow state after ~1 h BNB downtime, beam in and out (checklist)'),
    (t21_llrf4, t21_llrf4_settled, 'Width + Step Changing', 'LLRF4 to digital (LLRF test #4)'),
    (t22_bnb_study_min, t22_bnb_study_max, 'Stable', 'BNB bunch-width study 15:40-19:00 (DAQ log, AD ELOG)'),
    (datetime(2025, 6, 23, 10, 27, 58), datetime(2025, 6, 23, 10, 53, 18), '-', 'Source trouble, beam off / on (checklist)'),
    (datetime(2025, 6, 23, 13, 26, 39), datetime(2025, 6, 23, 14, 21, 2), '-', 'Source trouble, beam back (checklist)'),
    (t24_booster_start, t24_booster_end, '-', 'Booster study, beam stopped ~12:30 (checklist)'),
    (t26a_jun30_widening_start, t28a_waterpump_end, 'Width + Step Changing', '"Some instabilities in beam power" (checklist 18:38)'),
    (t28a_waterpump_end, t28b_waterpump_ramp_end, 'Step Changing', 'BRF Station 8 work, then Booster study (checklist)'),
    (t29_flatinjection_min, t29_flatinjection_max, 'Stable', 'Flat-injection study, filling +1 bucket (AD ELOG)'),
)


def chunk_labels(start, end, skip_reason=''):
    """Production label columns of one chunk [start, end] (`registry.chunk_labels`)."""
    return _reg.chunk_labels(CHUNK_NOTES, SEGMENTS, start, end, skip_reason)


# ===== Figures (overview / segmentation notebooks): transition chunks + beam-condition overlays =====
# - Not used by the calibration itself

t02_downtime_feb_start = datetime(2025, 2, 18, 13, 50, 0)     # beam downtime -- run 18255 last slice (DAQ "beam down" 13:54)
t02_downtime_feb_end = datetime(2025, 2, 18, 23, 9, 45)       # run 18258, beam back
t05_downtime_mar_end = datetime(2025, 3, 8, 7, 48, 49)        # run 18304, beam back
t09_downtime_apr08_start = datetime(2025, 4, 8, 8, 59, 39)    # beam downtime -- APSAIA power-cycle + network (runs 18407 / 18409)
t09_downtime_apr08_end = datetime(2025, 4, 8, 17, 50, 18)
t10_downtime_apr12_start = datetime(2025, 4, 12, 13, 28, 53)  # beam downtime -- no-beam config (run 18417)
t10_downtime_apr12_end = datetime(2025, 4, 13, 0, 47, 3)
t11_downtime_apr15_start = datetime(2025, 4, 15, 10, 59, 28)  # beam downtime -- beam shutdown + DAQ tests (runs 18423-18430)
t11_downtime_apr15_end = datetime(2025, 4, 15, 20, 55, 52)
t12_downtime_apr18_start = datetime(2025, 4, 18, 10, 59, 57)  # beam downtime -- DAQ test runs (18432-18438)
t12_downtime_apr18_end = datetime(2025, 4, 18, 20, 11, 22)
t13_downtime2_start = datetime(2025, 4, 21, 11, 1, 39)        # beam downtime -- run 18440 end
t13_downtime2_end = datetime(2025, 4, 21, 19, 12, 44)         # run 18442 start
t17_vacuum_outage_min = datetime(2025, 4, 30, 10, 7, 6)       # vacuum outage -- run 18445 end
t17_vacuum_outage_max = datetime(2025, 5, 1, 19, 20, 59)      # run 18458 start
t19_downtime_may14_start = datetime(2025, 5, 14, 9, 59, 59)   # beam downtime -- test runs 18480 / 18481
t19_downtime_may14_end = datetime(2025, 5, 14, 18, 18, 53)
t20_downtime_may28_start = datetime(2025, 5, 28, 12, 56, 47)  # beam downtime -- off-beam trigger tests (runs 18491-18502)
t20_downtime_may28_end = datetime(2025, 5, 28, 21, 44, 1)
t23_downtime_june_start = datetime(2025, 6, 23, 15, 59, 47)   # beam downtime -- run 18529 end
t23_downtime_june_end = datetime(2025, 6, 25, 12, 23, 25)     # run 18541 start
t25_downtime_jun26_start = datetime(2025, 6, 26, 13, 43, 59)  # beam downtime -- run 18543 end
t25_downtime_jun26_end = datetime(2025, 6, 26, 18, 1, 39)     # run 18550 start
t27_booster_start = datetime(2025, 6, 30, 20, 30, 0)          # Booster study, ECL "down for an hour" -- no SBND effect of its own
t27_booster_end = datetime(2025, 6, 30, 21, 0, 0)

# --- Transition windows: beam changing between two stable states (chunks drawn grey, not segmented) ---
TRANSITION_WINDOWS = (
    (t14_llrf1, t14a_llrf1_settled),
    (t15_llrf2_prelude_min, t15_llrf2_prelude_max),
    (t15a_spottybeam1_max, t16_llrf2),
    (t18_llrf3_flip_back, t18_llrf3_settled),   # flip-back piece only; 14:40 → 18:31 = stable digital state
    (t19a_may22_phase_flip_min, t19a_may22_phase_flip_max),
    (t21_llrf4, t21_llrf4_settled),
    (t26a_jun30_widening_start, t28b_waterpump_ramp_end),
)
STUDY_WINDOWS = (
    (t19b_may23_narrow_start, t19b_may23_narrow_end),   # narrow-bunch study
    (t22_bnb_study_min, t22_bnb_study_max),
    (t24_booster_start, t24_booster_end),
    (t27_booster_start, t27_booster_end),
    (t29_flatinjection_min, t29_flatinjection_max),
)
DOWNTIME_WINDOWS = (
    (t02_downtime_feb_start, t02_downtime_feb_end),
    (t05_downtime_mar_start, t05_downtime_mar_end),
    (t06_downtime_mar20_start, t06_downtime_mar20_end),
    (t09_downtime_apr08_start, t09_downtime_apr08_end),
    (t10_downtime_apr12_start, t10_downtime_apr12_end),
    (t11_downtime_apr15_start, t11_downtime_apr15_end),
    (t12_downtime_apr18_start, t12_downtime_apr18_end),
    (t13_downtime2_start, t13_downtime2_end),
    (t17_vacuum_outage_min, t17_vacuum_outage_max),
    (t19_downtime_may14_start, t19_downtime_may14_end),
    (t20_downtime_may28_start, t20_downtime_may28_end),
    (t23_downtime_june_start, t23_downtime_june_end),
    (t25_downtime_jun26_start, t25_downtime_jun26_end),
)
DOWNTIME_MIN_DURATION = timedelta(hours=24)   # shaded only if longer


def is_transition(start, end):
    """True if a chunk [start, end] lies inside one `TRANSITION_WINDOWS` entry."""
    return _reg.inside(start, end, TRANSITION_WINDOWS)


def transition_mask(df):
    """Boolean array over df rows (needs 'start' / 'end')."""
    return _reg.mask(df, is_transition)


def plot_studies(ax, label=False):
    """Orange band over each `STUDY_WINDOWS` entry; `label` adds one 'Accelerator Studies' legend entry."""
    for k, (a, b) in enumerate(STUDY_WINDOWS):
        ax.axvspan(a, b, color='tab:orange', alpha=0.3, lw=0, label='Accelerator Studies' if label and k == 0 else None)


def plot_llrf(ax):
    """Red dashed line at each LLRF switch."""
    ax.axvline(t14_llrf1, linestyle='--', lw=1, color='tab:red', label=str(t14_llrf1) + ' Switched to digital LLRF')
    ax.axvline(t16_llrf2, linestyle='--', lw=1, color='tab:red', label=str(t16_llrf2) + ' Switched back to analog LLRF')
    ax.axvline(t18_llrf3_switch, linestyle='--', lw=1, color='tab:red', label=str(t18_llrf3_switch) + ' Switched to digital LLRF (test)')
    ax.axvline(t21_llrf4, linestyle='--', lw=1, color='tab:red', label=str(t21_llrf4) + ' Switched to digital LLRF')


def plot_waterpump(ax):
    """Red dotted line at each late-June offset step."""
    for t in (t23a_waterpump_onset, t25a_waterpump_step2, t28a_waterpump_end):
        ax.axvline(t, linestyle=':', lw=1, color='tab:red')


def plot_bunch_rotation(ax):
    """Purple band over the rotation commissioning (= 'during' segment), solid line at its end."""
    ax.axvspan(t07_rotation_ramp_end, t07a_rotation_settled, color='tab:purple', alpha=0.3)
    ax.axvline(t07a_rotation_settled, linestyle='-', lw=1, color='tab:purple')


def plot_events(ax):
    """All Run 1 beam-condition overlays: rotation, LLRF, offset steps, studies."""
    plot_bunch_rotation(ax)
    plot_llrf(ax)
    plot_waterpump(ax)
    plot_studies(ax)


def plot_extended_downtime(ax):
    """Grey band over each `DOWNTIME_WINDOWS` entry longer than `DOWNTIME_MIN_DURATION`."""
    for t_min, t_max in DOWNTIME_WINDOWS:
        if t_max - t_min > DOWNTIME_MIN_DURATION:
            ax.axvspan(t_min, t_max, color='grey', alpha=0.3)
