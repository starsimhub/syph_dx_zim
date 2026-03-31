"""
Health impact analysis: GUD POC diagnostic vs standard of care.

2x4 panel figure comparing SOC and GUD scenarios across 8 key metrics,
shown as annualized time series with intervention year marked.
"""

import numpy as np
import sciris as sc
import pandas as pd
import pylab as pl
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'
INTV_YEAR = 2027
START_YEAR = 2015
END_YEAR = 2040

# Stock variables are averaged per year; flows are summed
STOCK_METRICS = {'syph_n_active'}


def load_scenario(scenario):
    """Load treatment outcomes for a scenario."""
    return sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df')


def annualize(df):
    """Convert monthly data to annual: sum flows, mean stocks."""
    df = df.copy()
    df['int_year'] = np.floor(df['year']).astype(int)
    stocks = df[df.metric.isin(STOCK_METRICS)]
    flows = df[~df.metric.isin(STOCK_METRICS)]
    ann_stocks = stocks.groupby(['scenario', 'par_idx', 'int_year', 'metric'])['value'].mean().reset_index()
    ann_flows = flows.groupby(['scenario', 'par_idx', 'int_year', 'metric'])['value'].sum().reset_index()
    out = pd.concat([ann_stocks, ann_flows])
    out.rename(columns={'int_year': 'year'}, inplace=True)
    return out


def get_ts(df, metric, start=START_YEAR, end=END_YEAR):
    """Get annual time series: median and 10th/90th percentiles across par_idx."""
    sub = df[(df.metric == metric) & (df.year >= start) & (df.year <= end)]
    grouped = sub.groupby('year')['value']
    return pd.DataFrame({
        'year': grouped.median().index,
        'median': grouped.median().values,
        'lo': grouped.quantile(0.1).values,
        'hi': grouped.quantile(0.9).values,
    })


def get_cumulative(df, metric, start=INTV_YEAR, end=END_YEAR):
    """Get cumulative value post-intervention, per par_idx."""
    sub = df[(df.metric == metric) & (df.year >= start) & (df.year <= end)]
    return sub.groupby('par_idx')['value'].sum()


def plot_health_impact():
    """2x4 figure: annualized time series for SOC vs GUD across 8 health metrics."""
    soc = annualize(load_scenario('soc'))
    gud = annualize(load_scenario('gud'))

    set_font(size=18)

    metrics = {
        'syph_n_active':                  'Active syphilis cases',
        'syph_new_infections':            'Incident syphilis cases',
        'syph_new_congenital':            'Congenital syphilis cases',
        'syph_idalys_dalys':              'Syphilis DALYs',
        'syph_new_treated':               'Treated',
        'syph_new_treated_success':       'Correctly treated',
        'syph_new_treated_unnecessary':   'Overtreated',
        'syph_new_false_neg':             'Missed diagnoses',
    }

    fig, axes = pl.subplots(2, 4, figsize=(24, 10))
    axes = axes.ravel()

    colors = {'soc': '#555555', 'gud': '#e41a1c'}
    labels = {'soc': 'SOC', 'gud': 'GUD POC'}

    for i, (metric, title) in enumerate(metrics.items()):
        ax = axes[i]
        for scen, df in [('soc', soc), ('gud', gud)]:
            ts = get_ts(df, metric)
            if len(ts) == 0:
                continue
            ax.plot(ts.year, ts['median'], color=colors[scen], label=labels[scen], linewidth=2)
            ax.fill_between(ts.year, ts.lo, ts.hi, color=colors[scen], alpha=0.2)
        ax.axvline(x=INTV_YEAR, color='k', ls='--', alpha=0.5)
        ax.set_title(title)
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)
        if i == 0:
            ax.legend(frameon=False)

    pl.tight_layout()
    fname = f'{FIGURES_DIR}/gud_health_impact.png'
    fig.savefig(fname, dpi=200, bbox_inches='tight')
    print(f'Saved {fname}')

    # Print summary table
    print(f'\nCumulative outcomes {INTV_YEAR}\u2013{END_YEAR} (median):')
    print(f'{"Metric":<30s} {"SOC":>10s} {"GUD":>10s} {"Change":>10s}')
    print('-' * 62)
    for metric, title in metrics.items():
        soc_val = get_cumulative(soc, metric).median()
        gud_val = get_cumulative(gud, metric).median()
        pct = (gud_val - soc_val) / max(soc_val, 1) * 100
        print(f'{title:<30s} {soc_val:>10,.0f} {gud_val:>10,.0f} {pct:>+9.1f}%')

    return fig


if __name__ == '__main__':
    fig = plot_health_impact()
    print('Done!')
