"""Plotting helpers of the figure notebooks (`calibration/overview.ipynb`, `calibration/segmentation.ipynb`).

- Importing applies the common rcParams (`STYLE`); the notebooks then set the figure font sizes
- Per-chunk series: one flat line + error band per chunk; grey = transition chunk of the run registry (`ev=`)
- Segment colours (`SEG_COLORS`) and the legend handles of the beam-condition markers (`CONDITION_HANDLES*`)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import NullLocator
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from calib_fit_helpers import asymmetric_gaussian
import event_timestamps_run1 as _ev1
import event_timestamps_run2 as _ev2

STYLE = {
    'font.size': 13,
    'axes.titlesize': 15,
    'axes.labelsize': 13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 16,
}
plt.rcParams.update(STYLE)

ASYM_COLORS = {'L': 'blue', 'R': 'red'}   # asym fit halves: σL blue, σR red
MAX_GAP_HOURS = 24.0                      # centre line broken across longer gaps

UNSTABLE_COLOR = 'k'   # `ev.UNSTABLE_SEGMENTS`: black, as unsegmented chunks
SEG_COLORS = {'before_a': 'tab:blue', 'before_c': 'tab:olive', 'during': 'tab:orange', 'after': 'tab:purple',
              'wide': 'tab:cyan', 'normal_a': 'tab:green', 'normal_b': 'tab:brown', 'narrow': 'tab:pink',
              **{s: UNSTABLE_COLOR for s in _ev1.UNSTABLE_SEGMENTS + _ev2.UNSTABLE_SEGMENTS}}

# legend handles of the beam-condition markers (`ev.plot_events`, `ev.plot_extended_downtime`, ...)
CONDITION_HANDLES = [
    Patch(color='grey', alpha=0.3, label='Beam Downtime >24h'),
    Patch(color='tab:purple', alpha=0.3, label='Bunch Rotation Commissioning'),
    Line2D([], [], color='tab:purple', linestyle='-', lw=1, label='Bunch Rotation Effective'),
    Line2D([], [], color='tab:red', linestyle='--', lw=1, label='LLRF Works'),
    Line2D([], [], color='tab:red', linestyle=':', lw=1, label='Water-Pump Offset Steps'),
    Patch(color='tab:orange', alpha=0.3, label='Accelerator Studies'),
    Line2D([], [], color='tab:grey', lw=2, label='Transition Chunk (not segmented)'),
]
CONDITION_HANDLES_RUN2 = [
    Patch(color='grey', alpha=0.3, label='Beam Downtime >24h'),
    Patch(color='tab:purple', alpha=0.3, label='Bunch Rotation ON / Narrow Bunch'),
    Patch(color='tab:orange', alpha=0.3, label='Accelerator Studies'),
    Line2D([], [], color='tab:red', linestyle=':', lw=1, label='Beam Faults / Reduced Beam'),
]


def plot_center_line(ax, df_sub, y_col, color, start_col='start', end_col='end'):
    """Line through chunk centres, broken across gaps > `MAX_GAP_HOURS`."""
    d = df_sub.dropna(subset=[y_col]).sort_values(start_col).reset_index(drop=True)
    if len(d) < 2:
        return
    centers = d[start_col] + (d[end_col] - d[start_col]) / 2
    gap_hours = (d[start_col] - d[end_col].shift(1)) / pd.Timedelta(hours=1)
    group = (gap_hours > MAX_GAP_HOURS).cumsum()
    for _, g in pd.DataFrame({'center': centers, 'y': d[y_col], 'group': group}).groupby('group'):
        ax.plot(g['center'], g['y'], color=color, lw=1, alpha=0.8, zorder=4)


def style_time_axis(ax, xlim=None, xlabel='Date-Month-Year', fmt='%d-%b-%y', minor_grid_alpha=0.4):
    """Date formatter, minor ticks, major + dotted minor grid."""
    if xlim is not None:
        ax.set_xlim(*xlim)
    if xlabel:
        ax.set_xlabel(xlabel)
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    ax.minorticks_on()
    ax.grid(True, which='major', alpha=0.3)
    ax.grid(True, which='minor', alpha=minor_grid_alpha, linestyle=':', linewidth=0.7)


def time_axis(ax, xlim, fmt='%b\n%Y'):
    """Month ticks on a date axis (`style_time_axis` + one major tick per month, no minor ticks)."""
    style_time_axis(ax, xlim, xlabel=None, fmt=fmt)
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.tick_params(which='both', top=False, right=True)
    ax.tick_params(which='major', length=4, width=0.8)
    ax.tick_params(which='minor', length=2)


def plot_chunk_series(ax, df, col, err_col=None, color=None, center=True, center_color='tab:blue', alpha=0.3, ev=_ev1):
    """One flat line + error band (`err_col`, default `col` + '_err') per chunk.
    - `color=None`: blue stable / grey transition (`ev.is_transition`); else one fixed colour
    - `center`: centre line over the stable (non-transition) chunks
    - `ev`: run registry of the chunks (transition windows)
    """
    err_col = err_col or col + '_err'
    for _, row in df.iterrows():
        if pd.isna(row[col]):
            continue
        start, end = row['start'], row['end']
        x_date = pd.date_range(start=start, end=end, periods=5)
        c = color or ('tab:grey' if ev.is_transition(start, end) else 'tab:blue')
        ax.plot(x_date, [row[col]] * len(x_date), color=c, lw=1.5, zorder=3)
        ax.fill_between(x_date, row[col] - row[err_col], row[col] + row[err_col], color=c, alpha=alpha, linewidth=0, zorder=3)
    if center:
        stable = df if color else df[~ev.transition_mask(df)]
        plot_center_line(ax, stable, col, color or center_color)


def draw_asym_fit(ax, f, label):
    """Asym Gaussian of a `stacking.stacked_bunch_fit` result `f`: σL half (x < x0) blue, σR half red.
    - `label(side, sigma, err)` → legend text per half (side 'L' / 'R')
    """
    pa, ea, xx = f['pa'], f['ea'], f['xx']
    for side, i in (('L', 3), ('R', 4)):
        xs = xx[xx < pa[2]] if side == 'L' else xx[xx >= pa[2]]
        ax.plot(xs, f['k'] * asymmetric_gaussian(xs, *pa), color=ASYM_COLORS[side], lw=1.5, label=label(side, abs(pa[i]), ea[i]))
