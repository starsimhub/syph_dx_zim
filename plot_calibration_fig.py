"""
Publication-quality calibration figure for supplementary materials.
Uses the calibration stats from the best-fit parameter sets.
"""
import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from utils import set_font

# Settings
RESULTS_DIR = 'results'
DATA_DIR = 'data'
FIGURES_DIR = 'figures'
START_YEAR = 1990
END_YEAR = 2025

# Colors
HIV_COLOR = '#2171b5'
SYPH_COLOR = '#cb181d'
DATA_COLOR = 'k'
BAND_ALPHA = 0.2


def load_data():
    cs = sc.loadobj(f'{RESULTS_DIR}/zimbabwe_calib_stats_all.df')
    data_hiv = pd.read_csv(f'{DATA_DIR}/zimbabwe_hiv_data.csv')
    data_syph = pd.read_csv(f'{DATA_DIR}/zimbabwe_syph_data.csv')
    return cs, data_hiv, data_syph


def get_stats(cs, col, lo='10%', hi='90%'):
    return cs[(col, '50%')], cs[(col, lo)], cs[(col, hi)]


def plot_panel(ax, cs, col, data_df=None, data_col=None, title='', color='C0',
               ylabel=None, si_ticks=False, ylim_bottom=True):
    """Plot a single calibration panel with model band and data points"""
    med, lo, hi = get_stats(cs, col)
    years = cs.index

    ax.fill_between(years, lo, hi, alpha=BAND_ALPHA, color=color, linewidth=0)
    ax.plot(years, med, color=color, linewidth=1.5, label='Model (median)')

    if data_df is not None and data_col is not None:
        d = data_df[['time', data_col]].dropna()
        if len(d):
            ax.scatter(d.time, d[data_col], color=DATA_COLOR, s=15, zorder=5,
                       label='Data', edgecolors='none')

    ax.set_title(title, fontsize=11, pad=4)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    ax.set_xlim(START_YEAR, END_YEAR)
    if ylim_bottom:
        ax.set_ylim(bottom=0)
    if si_ticks:
        sc.SIticks(ax, axis='y')
    ax.tick_params(labelsize=8)


def plot_calibration_figure():
    set_font(size=11)
    cs, data_hiv, data_syph = load_data()

    fig = plt.figure(figsize=(14, 14))
    gs = GridSpec(3, 3, left=0.08, right=0.97, bottom=0.05, top=0.94,
                  wspace=0.30, hspace=0.35)

    # --- Row 1: HIV ---
    ax = fig.add_subplot(gs[0, 0])
    plot_panel(ax, cs, 'hiv.prevalence_15_49', data_hiv, 'hiv.prevalence_15_49',
               'HIV prevalence (15-49)', HIV_COLOR, ylabel='Prevalence')
    ax.legend(fontsize=7, loc='upper right', frameon=False)

    ax = fig.add_subplot(gs[0, 1])
    plot_panel(ax, cs, 'hiv.new_infections', data_hiv, 'hiv.new_infections',
               'HIV new infections', HIV_COLOR, si_ticks=True)

    ax = fig.add_subplot(gs[0, 2])
    plot_panel(ax, cs, 'hiv.n_on_art', data_hiv, 'hiv.n_on_art',
               'People on ART', HIV_COLOR, si_ticks=True)

    # --- Row 2: HIV continued ---
    ax = fig.add_subplot(gs[1, 0])
    plot_panel(ax, cs, 'hiv.n_infected', data_hiv, 'hiv.n_infected',
               'People living with HIV', HIV_COLOR, si_ticks=True)

    ax = fig.add_subplot(gs[1, 1])
    plot_panel(ax, cs, 'hiv.new_deaths', data_hiv, 'hiv.new_deaths',
               'HIV-related deaths', HIV_COLOR, si_ticks=True)

    ax = fig.add_subplot(gs[1, 2])
    plot_panel(ax, cs, 'n_alive', data_hiv, 'n_alive',
               'Total population', '#555555', si_ticks=True)

    # --- Row 3: Syphilis ---
    ax = fig.add_subplot(gs[2, 0])
    plot_panel(ax, cs, 'syph.active_prevalence', data_syph, 'syph.active_prevalence',
               'Active syphilis prevalence', SYPH_COLOR, ylabel='Prevalence')

    ax = fig.add_subplot(gs[2, 1])
    plot_panel(ax, cs, 'syph.n_active', data_syph, 'syph.n_active',
               'Active syphilis cases', SYPH_COLOR, si_ticks=True)

    ax = fig.add_subplot(gs[2, 2])
    plot_panel(ax, cs, 'syph.new_infections', data_syph, 'syph.new_infections',
               'Syphilis new infections', SYPH_COLOR, si_ticks=True)

    # Panel labels
    labels = 'ABCDEFGHI'
    for i, ax in enumerate(fig.axes):
        ax.text(-0.12, 1.08, labels[i], transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top')

    fig.suptitle('Figure S2: Model calibration to Zimbabwe HIV and syphilis epidemiological data',
                 fontsize=12, y=0.98)

    plt.savefig(f'{FIGURES_DIR}/fig_s2_calibration.png', dpi=300, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig_s2_calibration.png')


if __name__ == '__main__':
    plot_calibration_figure()
