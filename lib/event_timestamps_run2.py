"""Run 2 registry: beam-condition timestamps, hard chunk boundaries, bad runs, stacking segments and chunk labels (same API as Run 1).

- Used by `calibration/calibration.ipynb` (`ev2`): `FINE_SPLITS`, `BAD_RUNS`, `chunk_labels`
- Run 2 in the ntuple: runs 19857-20219 (24-Jan → 12-Jun-2026); runs 19683-19854 at the start of Run 2 are not processed yet
- Timestamps in chronological order, numbered t04, t05, ...; letter suffixes keep the order; `[SPLIT]` = listed in `FINE_SPLITS`
- Split time = SBND hit after the change, truncated to whole seconds (`ts < split` goes before)
- All times UTC; AD times are 10-min grid points; ELOG / DAQ-log / ECL beam-checklist times (US Central) converted
- Sources: SBND ntuple (bunch position / width per chunk, z = step / error), AD RWM traces, Booster ELOG, SBND ECL beam checklist, DAQ run log
"""

from datetime import datetime, timedelta

import registry as _reg

t04b_brf06_trip = datetime(2026, 1, 29, 0, 28, 26)            # [SPLIT] Booster ELOG BRF06 trip 01-28 18:28 CST; SBND offset +0.59 → +0.32 ns
t04b_booster_rf_reset = datetime(2026, 1, 29, 2, 43, 58)      # [SPLIT] permit reset after BRF05-15 sparks 20:38-20:43 CST; run 19868 first hit, offset back +0.64 ns
t04d_brf19_disconnected = datetime(2026, 2, 4, 17, 8, 59)     # [SPLIT] Booster ELOG "BRF19 was disconnected from the anode during the Linac downtime"; SBND offset −0.34 ns, width flat
t04e_feb07_offset_step = datetime(2026, 2, 7, 16, 11, 27)     # [SPLIT] run 19887 -- +0.30 ns offset step, width flat, no beam gap; cause unknown
t04c_beam_rate_normal = datetime(2026, 2, 9, 22, 26, 51)      # [SPLIT] run 19891 -- rate 1.8 → 17.9 hits/min after the "beam tests of the day" (DAQ 19887 / 19891)

t05_feb17_downtime_min = datetime(2026, 2, 17, 20, 30, 0)     # announced downtime start; ends the Jan-Feb high-σ period = wide / normal segment boundary
t05_feb17_beam_return = datetime(2026, 2, 17, 22, 34, 34)     # [SPLIT] run 19911 first hit after the 122.5-min beam gap of the announced downtime

t05a_run19913_start = datetime(2026, 2, 19, 6, 11, 2)         # [SPLIT] run 19913 first hit
t05a_run19913_narrow_end = datetime(2026, 2, 19, 8, 26, 0)    # [SPLIT] end of the narrow start of run 19913 (SBND σ 1.39 vs ~1.8 after, AD flat); not segmented
t05a_run19913_end = datetime(2026, 2, 19, 20, 2, 33)          # run 19913 last hit; the run stopped for the Booster study
t05b_narrow_bunch0_on = datetime(2026, 2, 19, 20, 19, 15)     # [SPLIT] run 19914 first hit -- Booster "Narrow Bunch Pre-Tune Study"; AD σ 1.34 → 0.87
t05b_narrow_bunch0_off = datetime(2026, 2, 19, 21, 22, 57)    # [SPLIT] run 19917 first hit = first hit after "Studies are complete"

t05b_beam_return_feb21 = datetime(2026, 2, 21, 17, 49, 20)    # [SPLIT] beam back after the 29 h downtime (run 19923 first hit)
t05b_beam_return_feb21_settled = datetime(2026, 2, 21, 19, 9, 20)  # [SPLIT] settling ends: offset +0.25 ns, width −0.21 ns

t06_narrow_bunch1_on = datetime(2026, 3, 30, 17, 37, 39)      # [SPLIT] narrow bunch 1 on: AD σ 1.21 → 0.70; SBND width −0.20 ns
t06_narrow_bunch1_gap_min = datetime(2026, 3, 30, 18, 19, 41) # back to normal inside run 20011: AD σ 0.71 → 1.22
t06_narrow_bunch1_gap_max = datetime(2026, 3, 30, 19, 55, 23) # narrow again: AD σ 1.38 → 0.74
t06_narrow_bunch1_beam_return = datetime(2026, 3, 30, 20, 44, 50)   # [SPLIT] run 20011 first hit after a 28.9-min beam gap inside narrow piece 2
t06_narrow_bunch1_off = datetime(2026, 3, 30, 21, 50, 41)     # [SPLIT] narrow bunch 1 off: AD σ 0.75 → 1.33; SBND width +0.26 ns
t06a_narrow_bunch2_on = datetime(2026, 3, 31, 14, 20, 28)     # [SPLIT] narrow bunch 2 on: AD σ 1.33 → 0.69; DAQ "~30 min before the end of 20014"
t06a_narrow_bunch2_sbnd_narrow = datetime(2026, 3, 31, 14, 30, 2)   # [SPLIT] first hit after 14:30 (run 20014): SBND narrow from here
t06a_narrow_bunch2_off = datetime(2026, 4, 2, 19, 11, 43)     # [SPLIT] run 20020 first hit -- narrow bunch 2 off: AD σ 0.83 → 1.24; SBND width +0.11 ns

t07b_bnf_network_start = datetime(2026, 4, 8, 22, 3, 21)      # DAQ 20055 "BNF went down due to a network switch issue" -- last SBND hit
t08a_power_outage_start = datetime(2026, 4, 27, 18, 8, 5)     # DAQ 20084 "Sitewide Power Outage" -- last SBND hit
t09b_downtime_may11_start = datetime(2026, 5, 11, 14, 58, 54) # beam downtime (run 20118 last hit); normal_a / normal_b segment boundary

t09c_beam_return = datetime(2026, 5, 20, 16, 51, 3)           # [SPLIT] beam back after a 5-min beam gap (run 20139), width high ~1 h after
t09c_beam_return_settled = datetime(2026, 5, 20, 18, 1, 3)    # [SPLIT] width settled (−0.28 ns)
t09c_laser_notcher_fixed = datetime(2026, 5, 21, 1, 39, 27)   # [SPLIT] Booster laser notcher off ~19:30 → back 20:29 CDT; SBND offset ≈ −0.2 ns before

t10_comed_outage_max = datetime(2026, 5, 30, 1, 16, 21)       # first SBND hit after the planned ComEd outages (run 20179)
t10a_cep400_beam_return = datetime(2026, 6, 3, 2, 40, 53)     # [SPLIT] beam back after the B:CEP400 fuse fix (Booster ELOG)
t10a_injection_tune = datetime(2026, 6, 3, 4, 10, 56)         # [SPLIT] Booster injection tune; width ≈ 2.2 → 2.4 ns
t10b_llrf_study_start = datetime(2026, 6, 3, 13, 51, 19)      # [SPLIT] Booster LLRF study, new frequency curve (ELOG); SBND offset 0.16 → 0.74 ns
t10b_new_freq_curve = datetime(2026, 6, 3, 17, 44, 48)        # [SPLIT] run 20200 first hit -- new LLRF frequency curve in operation, offset flat ≈ +1.0 ns
t10b_freq_curve_end = datetime(2026, 6, 4, 23, 37, 9)         # [SPLIT] run 20202 -- SBND offset 1.04 → 0.18 ns, width flat; cause unknown
t10e_frf_work_burst = datetime(2026, 6, 10, 16, 54, 3)        # [SPLIT] Booster FRF frequency-curve work: offset +1.8 → +2.6 ns for 13 min (too few entries → skipped)
t10e_frf_loaded = datetime(2026, 6, 10, 17, 34, 47)           # [SPLIT] first hit after a 28-min gap, new FRF loaded


# --- Hard chunk boundaries: (run, split time), chronological ---
# - Narrow-bunch edges at AD 10-min resolution, refined to the SBND hit with the largest width step inside the AD bracket
FINE_SPLITS = (
    (19867, t04b_brf06_trip),
    (19868, t04b_booster_rf_reset),
    (19874, t04d_brf19_disconnected),
    (19887, t04e_feb07_offset_step),
    (19891, t04c_beam_rate_normal),
    (19911, t05_feb17_beam_return),
    (19913, t05a_run19913_start),
    (19913, t05a_run19913_narrow_end),
    (19914, t05b_narrow_bunch0_on),
    (19917, t05b_narrow_bunch0_off),
    (19923, t05b_beam_return_feb21),
    (19923, t05b_beam_return_feb21_settled),
    (20011, t06_narrow_bunch1_on),
    (20011, t06_narrow_bunch1_beam_return),
    (20011, t06_narrow_bunch1_off),
    (20014, t06a_narrow_bunch2_on),
    (20014, t06a_narrow_bunch2_sbnd_narrow),
    (20020, t06a_narrow_bunch2_off),
    (20139, t09c_beam_return),
    (20139, t09c_beam_return_settled),
    (20139, t09c_laser_notcher_fixed),
    (20192, t10a_cep400_beam_return),
    (20192, t10a_injection_tune),
    (20192, t10b_llrf_study_start),
    (20200, t10b_new_freq_curve),
    (20202, t10b_freq_curve_end),
    (20218, t10e_frf_work_burst),
    (20218, t10e_frf_loaded),
)

# --- Bad runs (dropped before chunking) ---
BAD_RUNS = (
    19888,   # 02-09 "beam tests of the day", 791 entries
    19890,   # 14 entries (< 20)
    20180,   # 17 entries (< 20)
    20203,   # 8 entries (< 20)
    20204,   # 11 entries (< 20)
)

# --- Narrow-bunch (bunch rotation on) periods, SBND-level edges ---
NARROW_BUNCH = (
    (t05a_run19913_start, t05a_run19913_narrow_end),
    (t05b_narrow_bunch0_on, t05b_narrow_bunch0_off),
    (t06_narrow_bunch1_on, t06_narrow_bunch1_gap_min),
    (t06_narrow_bunch1_gap_max, t06_narrow_bunch1_off),
    (t06a_narrow_bunch2_on, t06a_narrow_bunch2_off),
)

# --- Stacking segments: wide / normal_a / normal_b / narrow ---
# - wide: chunk starts before the 02-17 downtime (Jan-Feb high σ)
# - normal: everything after, except narrow-bunch chunks; cut at the 05-11 downtime into `NORMAL_PIECES`
# - narrow: chunk overlaps a `NARROW_BUNCH` window
# - Own beam state with only a few chunks → no segment (`NO_SEGMENT_WINDOWS`)
NORMAL_PIECES = {
    'normal_a': (None, t09b_downtime_may11_start),
    'normal_b': (t09b_downtime_may11_start, None),
}
NO_SEGMENT_WINDOWS = {
    'run19913_narrow': (t05a_run19913_start, t05a_run19913_narrow_end),
}


def is_unsegmented(start, end):
    """True if a chunk [start, end] lies inside a `NO_SEGMENT_WINDOWS` entry."""
    return _reg.inside(start, end, NO_SEGMENT_WINDOWS.values())


def stack_segment(start, end):
    """'wide' / 'normal_a' / 'normal_b' / 'narrow' for a chunk [start, end] (narrow wins, then the 02-17 downtime by chunk start, then
    `NORMAL_PIECES` by chunk start); None inside `NO_SEGMENT_WINDOWS`."""
    if is_unsegmented(start, end):
        return None
    if any(start < b and end > a for a, b in NARROW_BUNCH):
        return 'narrow'
    if start < t05_feb17_downtime_min:
        return 'wide'
    return next(name for name, (lo, hi) in NORMAL_PIECES.items() if (lo is None or start >= lo) and (hi is None or start < hi))


# --- Production label columns (as Run 1) ---
# - segment: 1-4 = wide / normal_a / normal_b / narrow; None = no constants or not segmented
# - CHUNK_NOTES = (start, end, bunch_condition, logs): reviewed chunks + skipped chunks; a chunk inside a window takes its label
#   (Width Changing = bunch width, Step Changing = bunch position); otherwise fitted → 'Stable', skipped → '-'
SEGMENT_NUMBER = {'wide': 1, 'normal_a': 2, 'normal_b': 3, 'narrow': 4}
CHUNK_NOTES = (
    (datetime(2026, 2, 17, 10, 3, 18), datetime(2026, 2, 17, 20, 32, 6), 'Stable', 'Announced downtime from ~20:30 follows (checklist 19:25)'),
    (t05_feb17_beam_return, t05a_run19913_start, 'Stable', 'Beam back after the announced 02-17 downtime (checklist)'),
    (t05a_run19913_start, t05a_run19913_narrow_end, 'Stable', 'Narrow bunch study'),
    (t05a_run19913_narrow_end, t05a_run19913_end, 'Stable', 'Normal beam; run stopped for the Booster study at 20:03 (DAQ log, Booster ELOG)'),
    (t06_narrow_bunch1_on, t06_narrow_bunch1_off, 'Width Changing', 'Narrow bunch on / off / on inside run 20011 (ECL)'),
    (t06a_narrow_bunch2_on, t06a_narrow_bunch2_sbnd_narrow, 'Width Changing', 'Narrow bunch switch-on, "narrowing began ~30 min before the end" (DAQ 20014)'),
    (t06a_narrow_bunch2_sbnd_narrow, t06a_narrow_bunch2_off, 'Stable', '"Narrow beam bunch study" (DAQ 20015)'),
    (t06a_narrow_bunch2_off, datetime(2026, 4, 3, 13, 20, 26), 'Stable', 'Back to normal configuration (checklist 19:33)'),
    # skipped chunks: labelled where SBND was reviewed, else '-' + log
    (datetime(2026, 1, 27, 13, 42, 47), datetime(2026, 1, 27, 14, 19, 2), '-', 'Beam down > 50 min, MI8 horizontal dipole trim PS failure (checklist)'),
    (datetime(2026, 1, 27, 15, 20, 30), datetime(2026, 1, 27, 16, 42, 35), '-', 'Linac issues, beam outside nominal conditions (checklist)'),
    (t04b_brf06_trip, t04b_booster_rf_reset, 'Step Changing', 'BRF06 trip, offset down until the Booster RF reset (Booster ELOG)'),
    (datetime(2026, 2, 9, 21, 14, 39), t04c_beam_rate_normal, '-', 'Beam tests of the day, beam normal from ~22:30 (DAQ 19887 / 19891)'),
    (datetime(2026, 2, 11, 19, 17, 53), datetime(2026, 2, 11, 20, 58, 16), '-', 'Flat-injection ramp curve changed (~11:50 CST) before; Booster RF down for a short downtime at 15:00 CST (Booster ELOG)'),
    (t05b_narrow_bunch0_on, datetime(2026, 2, 19, 20, 58, 36), 'Stable', 'Narrow Bunch Pre-Tune Study (Booster ELOG)'),
    (t05b_narrow_bunch0_off, datetime(2026, 2, 19, 21, 25, 24), '-', 'First run after "Studies are complete" (Booster ELOG, DAQ 19917)'),
    (t05b_beam_return_feb21, t05b_beam_return_feb21_settled, 'Width + Step Changing', 'Beam back after the 29 h downtime, DAQ restarted by experts (checklist, DAQ 19923)'),
    (datetime(2026, 3, 23, 16, 54, 48), datetime(2026, 3, 23, 19, 23, 0), '-', 'Beam back after Linac work (~11:45 CDT); high Booster 026 losses intermittently disabling beam from ~14:15 CDT (checklist, Booster ELOG)'),
    (datetime(2026, 4, 3, 17, 10, 0), datetime(2026, 4, 3, 19, 31, 28), '-', 'Beam back ~noon, DAQ rack issue, Linac issue after (checklist)'),
    (datetime(2026, 4, 5, 8, 27, 11), datetime(2026, 4, 5, 10, 36, 0), '-', 'Linac issue, beam back but intermittent (checklist)'),
    (datetime(2026, 4, 8, 19, 56, 18), t07b_bnf_network_start, '-', 'BNF network switch issue, beam control lost (DAQ 20055)'),
    (datetime(2026, 4, 27, 17, 2, 37), t08a_power_outage_start, '-', 'Beam back after Booster access, sitewide power outage after (checklist, DAQ 20084)'),
    (datetime(2026, 5, 19, 20, 56, 6), datetime(2026, 5, 19, 21, 45, 45), '-', 'Driver issue on BNB line, beam in and out (checklist)'),
    (t09c_beam_return, t09c_beam_return_settled, 'Width Changing', 'BNB back after downtime (checklist)'),
    (datetime(2026, 5, 24, 23, 16, 59), datetime(2026, 5, 24, 23, 23, 21), '-', 'BNB back after > 5 h, stabilising (checklist)'),
    (t10_comed_outage_max, datetime(2026, 5, 30, 1, 18, 7), '-', 'Beam back after the ComEd outage, CRT ramp-up test run (DAQ 20179)'),
    (t10a_cep400_beam_return, t10a_injection_tune, 'Width Changing', 'Beam back after a Booster PS fault, until the injection tune (checklist, Booster ELOG)'),
    (t10b_llrf_study_start, datetime(2026, 6, 3, 16, 10, 7), 'Step Changing', 'Booster LLRF studies, new frequency curve (Booster ELOG)'),
    (t10e_frf_work_burst, t10e_frf_loaded, 'Step Changing', 'Booster FRF frequency-curve work (Booster ELOG)'),
    (t10e_frf_loaded, datetime(2026, 6, 11, 6, 26, 43), 'Stable', 'New FRF loaded, responses reverted over several supercycles (Booster ELOG)'),
)


def chunk_labels(start, end, skip_reason=''):
    """Production label columns of one chunk [start, end] (`registry.chunk_labels`)."""
    return _reg.chunk_labels(CHUNK_NOTES, SEGMENT_NUMBER, stack_segment, start, end, skip_reason)


# ===== Figures (overview / segmentation notebooks): transition chunks + beam-condition overlays =====
# - Not used by the calibration itself

t01_downtime_nov_start = datetime(2025, 11, 5, 16, 40, 0)      # 2-week downtime (AD grid gap); DAQ 19689 "power outage in Batavia site"
t01_downtime_nov_end = datetime(2025, 11, 19, 4, 40, 0)        # DAQ 19712 "The beam is back!!"
t01a_bad_quality_beam_min = datetime(2025, 11, 19, 6, 21, 14)  # DAQ 19714 "This run has bad quality beam" (whole run)
t02_linac_low_intensity = datetime(2025, 12, 19, 6, 0, 0)      # LINAC station trouble, MCR "lower intensity" (AD σ step −0.13 ns)
t04c_beam_tests_min = datetime(2026, 2, 9, 15, 56, 37)         # DAQ 19887 "stopped ... for the beam tests of the day" (runs 19888-19890)
t04c_beam_tests_max = t04c_beam_rate_normal                   # DAQ 19891 "Beam returned to normal", refined to the SBND rate change
t05c_downtime_feb23_start = datetime(2026, 2, 23, 14, 53, 35)  # beam downtime (14 days; ARARA cable work, APSAIA) -- last SBND hit, run 19923
t05c_downtime_feb23_end = datetime(2026, 3, 9, 18, 27, 12)     # run 19947 first hit, "BNB returning"
t05d_downtime_mar15_start = datetime(2026, 3, 15, 9, 0, 0)     # beam downtime (HVAC maintenance, DAQ work) -- last SBND hit, run 19952
t05d_downtime_mar15_end = datetime(2026, 3, 21, 0, 54, 47)     # run 19996 first hit, "BNB is back!"
t07_linac_outage_min = datetime(2026, 4, 5, 6, 30, 0)          # LINAC outage (AD grid gaps; σ −0.16 ns across)
t07_linac_outage_max = datetime(2026, 4, 6, 4, 10, 0)          # "beam back but intermittent", "may not be full intensity"
t07b_bnf_network_end = datetime(2026, 4, 11, 17, 26, 34)       # first hit, run 20067
t08a_power_outage_end = datetime(2026, 5, 1, 20, 27, 4)        # first hit, run 20102
t09_substation_trip = datetime(2026, 5, 2, 0, 30, 0)           # checklist: substation trip, beam down
t09b_downtime_may11_end = datetime(2026, 5, 13, 19, 4, 35)     # first hit, run 20128
t09d_linac_trip = datetime(2026, 5, 22, 2, 58, 14)             # DAQ 20142 "BNB downtime, linac trip and no ETA" -- last hit before the 3-h gap
t10_comed_outage_min = datetime(2026, 5, 26, 14, 59, 51)       # planned beam-off for ComEd outages (DAQ 20153) -- last SBND hit
t10e_intermittent_beam_min = datetime(2026, 6, 10, 13, 2, 45)  # DAQ 20218 "Intermittent beam (one hour down)"

SEGMENT_A_RUNS = (19683, 19854)   # start of Run 2, not processed yet (red band)

# --- Transition windows: beam changing between two stable states (chunks drawn grey) ---
TRANSITION_WINDOWS = (
    (t05b_beam_return_feb21, t05b_beam_return_feb21_settled),
    (t09c_beam_return, t09c_beam_return_settled),
    (t10a_cep400_beam_return, t10a_injection_tune),
    (t10b_llrf_study_start, t10b_new_freq_curve),
    (t10e_frf_work_burst, t10e_frf_loaded),
)
UNSTABLE_SEGMENTS = ()
STUDY_WINDOWS = (
    (t04c_beam_tests_min, t04c_beam_tests_max),
)
BEAM_FAULTS = (
    t01a_bad_quality_beam_min, t02_linac_low_intensity, t07_linac_outage_min, t07_linac_outage_max,
    t09_substation_trip, t09d_linac_trip, t10e_intermittent_beam_min,
)
DOWNTIME_WINDOWS = (
    (t01_downtime_nov_start, t01_downtime_nov_end),
    (datetime(2026, 2, 5, 12, 29, 59), datetime(2026, 2, 6, 17, 3, 7)),     # access + beam off (runs 19876 → 19884)
    (datetime(2026, 2, 20, 12, 57, 36), datetime(2026, 2, 21, 17, 49, 20)), # MCR > 2 h, extended (19921 → 19923)
    (t05c_downtime_feb23_start, t05c_downtime_feb23_end),
    (t05d_downtime_mar15_start, t05d_downtime_mar15_end),
    (t07b_bnf_network_start, t07b_bnf_network_end),
    (t08a_power_outage_start, t08a_power_outage_end),
    (t09b_downtime_may11_start, t09b_downtime_may11_end),
    (t10_comed_outage_min, t10_comed_outage_max),
)
DOWNTIME_MIN_DURATION = timedelta(hours=24)   # shaded only if longer


def is_transition(start, end):
    """True if a chunk [start, end] lies inside one `TRANSITION_WINDOWS` entry."""
    return _reg.inside(start, end, TRANSITION_WINDOWS)


def transition_mask(df):
    """Boolean array over df rows (needs 'start' / 'end')."""
    return _reg.mask(df, is_transition)


def plot_narrow_bunch(ax):
    """Purple band over each `NARROW_BUNCH` period."""
    for a, b in NARROW_BUNCH:
        ax.axvspan(a, b, color='tab:purple', alpha=0.3, lw=0)


def plot_studies(ax, label=False):
    """Orange band over each `STUDY_WINDOWS` entry; `label` adds one 'Accelerator Studies' legend entry."""
    for k, (a, b) in enumerate(STUDY_WINDOWS):
        ax.axvspan(a, b, color='tab:orange', alpha=0.3, lw=0, label='Accelerator Studies' if label and k == 0 else None)


def plot_beam_faults(ax):
    """Red dotted line at each `BEAM_FAULTS` time."""
    for t in BEAM_FAULTS:
        ax.axvline(t, linestyle=':', lw=1, color='tab:red')


def plot_extended_downtime(ax):
    """Grey band over each `DOWNTIME_WINDOWS` entry longer than `DOWNTIME_MIN_DURATION`."""
    for a, b in DOWNTIME_WINDOWS:
        if b - a > DOWNTIME_MIN_DURATION:
            ax.axvspan(a, b, color='grey', alpha=0.3, lw=0)
