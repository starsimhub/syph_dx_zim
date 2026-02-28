"""
Additional syphilis epidemiology diagnostics (optional/exploratory, not in manuscript).

Panels:
    A: All syphilis prevalence by age and sex work status
    B: Ever exposed to syphilis (time series by sex/FSW)
    C: Syphilis prevalence by HIV status (bars)
"""

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font
from plot_fig1_epi import (
    plot_prev_by_sw, plot_syph_hiv_ratio_ts,
    F_COLOR, M_COLOR, F_COLOR_LIGHT, M_COLOR_LIGHT,
    HIV_POS_COLOR, HIV_NEG_COLOR,
    ZIMPHIA_SYPH, ZIMPHIA_SYPH_EVER_BY_AGE, WHO_FSW_SYPH_PREV,
)

BAND_ALPHA = 0.2


def plot_ever_exposed_ts(cs, ax, start_year=1995):
    """Ever exposed to syphilis over time, by sex and FSW status"""
    years = cs.index[cs.index >= start_year]

    for col, color, label in [
        ('epi_ts.syph_ever_exposed_fsw', F_COLOR, 'FSW'),
        ('epi_ts.syph_ever_exposed_f', F_COLOR_LIGHT, 'All girls/women'),
        ('epi_ts.syph_ever_exposed_m', M_COLOR_LIGHT, 'All boys/men'),
    ]:
        med = cs.loc[years, (col, '50%')] * 100
        lo = cs.loc[years, (col, '10%')] * 100
        hi = cs.loc[years, (col, '90%')] * 100
        ax.fill_between(years, lo, hi, alpha=0.15, color=color)
        ax.plot(years, med, color=color, label=label, linewidth=2)

    # ZIMPHIA 2016 ever-exposed data points
    kw = dict(marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5)
    ax.scatter([2016], [ZIMPHIA_SYPH['ever_exposed_f'] * 100], color=F_COLOR_LIGHT,
               label='ZIMPHIA 2016', **kw)
    ax.scatter([2016], [ZIMPHIA_SYPH['ever_exposed_m'] * 100], color=M_COLOR_LIGHT, **kw)

    # WHO GHO FSW syphilis seroprevalence (~30%, treponemal = ever-exposed)
    ax.scatter([2016], [WHO_FSW_SYPH_PREV * 100], color=F_COLOR,
               marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5, label='WHO GHO FSW')

    ax.legend(frameon=False, fontsize=11, loc='upper left')
    ax.set_ylabel('Prevalence (%)')
    ax.set_title('Ever exposed to syphilis')
    ax.set_ylim(bottom=0)
    ax.set_xlim(start_year, 2025)


def plot_syph_prev_by_hiv_bars(cs, ax, year=2016):
    """Overall syphilis prevalence by HIV status — bars with error bars + ZIMPHIA"""
    x = np.array([0, 1])
    labels = ['HIV+', 'HIV\u2212']
    colors = [HIV_POS_COLOR, HIV_NEG_COLOR]
    cols = ['coinfection_stats.syph_prev_has_hiv', 'coinfection_stats.syph_prev_no_hiv']

    meds = []
    los = []
    his = []
    for col in cols:
        med = cs.loc[year, (col, '50%')] * 100
        lo = cs.loc[year, (col, '10%')] * 100
        hi = cs.loc[year, (col, '90%')] * 100
        meds.append(med)
        los.append(med - lo)
        his.append(hi - med)

    ax.bar(x, meds, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5, width=0.6)
    ax.errorbar(x, meds, yerr=[los, his], fmt='none', ecolor='k', capsize=5, linewidth=1.5)

    # ZIMPHIA data as diamonds
    kw = dict(marker='D', s=80, zorder=5, edgecolors='k', linewidths=0.5)
    ax.scatter(0, ZIMPHIA_SYPH['active_hiv_pos'] * 100, color=HIV_POS_COLOR,
               label='ZIMPHIA 2016', **kw)
    ax.scatter(1, ZIMPHIA_SYPH['active_hiv_neg'] * 100, color=HIV_NEG_COLOR, **kw)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Syphilis\nprevalence (%)')
    ax.set_title('Syphilis prevalence\nby HIV status')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=11, loc='upper right')


if __name__ == '__main__':

    sw_prev_df = sc.loadobj('results/sw_prev_df.df')
    cs = sc.loadobj('results/zimbabwe_calib_stats_all.df')

    set_font(size=18)
    fig = pl.figure(figsize=(22, 7))
    gs = GridSpec(1, 3, left=0.06, right=0.98, bottom=0.12, top=0.88,
                  wspace=0.28)

    ax = fig.add_subplot(gs[0, 0])
    plot_prev_by_sw(sw_prev_df, ax=ax, prev_col='prevalence',
                    title='All syphilis prevalence\nby sex work status')

    ax = fig.add_subplot(gs[0, 1])
    plot_ever_exposed_ts(cs, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_syph_prev_by_hiv_bars(cs, ax)

    for i, ax in enumerate(fig.axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig('figures/fig_s1_epi.png', dpi=200, bbox_inches='tight')
    print('Saved figures/fig_s1_epi.png')
