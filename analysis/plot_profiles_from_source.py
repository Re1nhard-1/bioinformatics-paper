from pathlib import Path
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scientific_style import configure_style, figure_mm, COLORS

ROOT = Path(__file__).resolve().parents[1]

def save(fig, output, name):
    for suffix in ['pdf', 'png', 'svg']:
        fig.savefig(output / f'{name}.{suffix}', dpi=300)
    plt.close(fig)

def original(output):
    data = pd.read_csv(ROOT / 'results/music_gradient_figures/20260917T105111607303Z/fixed_budget_gradient.source_data.csv', dtype={'held_out': str})
    configure_style()
    fig = figure_mm(183, 165)
    axes = fig.subplots(2, 2)
    fig.subplots_adjust(left=.11, right=.97, bottom=.18, top=.83, hspace=.58, wspace=.30)
    fig.text(.11, .965, 'Reference allocation at a fixed budget of 60 cells per type', weight='bold')
    fig.text(.11, .933, 'Grey: individual held-out donors; coloured: equal-donor mean')
    donors = sorted(data.loc[data.summary.eq('donor'), 'held_out'].unique(), key=int)
    markers = ['o', 's', '^', 'v', 'D', 'P']
    fig.legend(handles=[Line2D([], [], color='#8a8a8a', marker=m, mfc='white', lw=.6, ms=3, label=d) for d, m in zip(donors, markers)], ncol=6, loc='upper left', bbox_to_anchor=(.11, .9), columnspacing=1.7)
    for i, ax in enumerate(axes.flat):
        rows = data[data.panel.eq('abcd'[i])]
        for donor, marker in zip(donors, markers):
            y = rows[rows.held_out.eq(donor)].sort_values('level').value_pp
            ax.plot(range(4), y, color='#969696', marker=marker, mfc='white', mec='#777777', lw=.55, ms=3, mew=.5)
        means = rows[rows.summary.eq('equal_donor_mean')].sort_values('level').value_pp
        color = COLORS['vermillion'] if i % 2 == 0 else COLORS['blue']
        ax.plot(range(4), means, color=color, marker='s' if i % 2 == 0 else 'o', lw=1.1, ms=4)
        bounds = data[data.quantity.eq('absolute_mae' if i < 2 else 'mae_minus_balanced')].value_pp
        if i < 2:
            limits = (np.floor(bounds.min()*2)/2-.2, np.ceil(bounds.max()*2)/2+.2)
        else:
            limits = (min(-.08, bounds.min()-.08), max(.1, bounds.max()+.08))
        ax.set(xticks=range(4), xticklabels=['20:20:20', '30:15:15', '40:10:10', '50:5:5'], xlim=(-.15, 3.15), ylim=limits)
        ax.set_title('Weighted MuSiC' if i % 2 == 0 else 'Ordinary NNLS', loc='left', pad=10)
        ax.text(-.18, 1.065, 'abcd'[i], transform=ax.transAxes, fontsize=8, weight='bold')
        ax.set_ylabel(('Absolute MAE' if i < 2 else 'MAE minus balanced') + '\n(percentage points)')
        ax.set_yticks([v for v in ax.get_yticks() if limits[0] <= v <= limits[1]])
        if i >= 2:
            ax.axhline(0, color=COLORS['neutral'], lw=.6, ls='--')
            ax.set_xlabel('Cells per donor within each type', labelpad=7)
    fig.text(.11, .081, 'All six donor profiles are shown; segments connect tested allocations only.')
    fig.text(.11, .048, 'Shared reference donors; synthetic mixtures; three draws; no population error bars.')
    save(fig, output, 'Figure_S1')

def external(output):
    data = pd.read_csv(ROOT / 'results/music_external_figures/20260917T133540882209Z/external_absolute_profiles.source_data.csv')
    configure_style(7)
    fig = figure_mm(183, 134)
    axes = fig.subplots(2, 3)
    fig.subplots_adjust(left=.085, right=.985, bottom=.13, top=.92, wspace=.24, hspace=.40)
    levels = ['balanced', 'ratio2', 'ratio4', 'ratio10']
    ymax = data.mae_pp.max()*1.08
    for i, ax in enumerate(axes.flat):
        rows = data[data.panel.eq('abcdef'[i])]
        for _, group in rows[rows.summary.eq('donor')].groupby('held_out', sort=False):
            ax.plot(range(4), group.set_index('level').loc[levels, 'mae_pp'], color='#B5B5B5', linewidth=.5, alpha=.85)
        mean = rows[rows.summary.eq('mean')].set_index('level').loc[levels, 'mae_pp']
        ax.plot(range(4), mean, color=COLORS['vermillion'] if i < 3 else COLORS['blue'], marker='s' if i < 3 else 'o', linewidth=1, markersize=3)
        method = 'Weighted MuSiC' if i < 3 else 'Ordinary NNLS'
        ax.set(xlim=(-.15, 3.15), ylim=(0, ymax), xticks=range(4), xticklabels=['1:1:1','2:1:1','4:1:1','10:1:1'], title=f'{method}\n{[60,120,300][i%3]} reference cells / type')
        ticks = MaxNLocator(nbins=4).tick_values(0, ymax)
        ax.set_yticks(ticks[(ticks >= 0) & (ticks <= ymax)])
        ax.text(-.17, 1.08, 'abcdef'[i], transform=ax.transAxes, fontsize=8, fontweight='bold', va='top')
        if i % 3 == 0: ax.set_ylabel('Mean absolute error (pp)')
        if i >= 3: ax.set_xlabel('Donor allocation ratio')
    fig.text(.085, .025, 'Grey lines: all 14 donors. Coloured line: equal-donor mean. Shared y scale; no population error bars.', fontsize=7)
    save(fig, output, 'Figure_S2')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT.parent/'figures')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    original(args.output)
    external(args.output)
