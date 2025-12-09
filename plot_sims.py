"""
Plotting functions for HIV-syphilis coinfection model calibration and results.

This module provides functions to visualize model outputs including:
- Individual simulation runs
- Calibration results with uncertainty bounds
- HIV and syphilis epidemiological trends
- Coinfection dynamics
"""

import numpy as np
import pandas as pd
import pylab as pl
import sciris as sc
import seaborn as sns
from utils import set_font, get_y, plot_single, percentile_pairs


# Constants
LOCATION = 'zimbabwe'
DATA_DIR = 'data'
RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'


def plot_sims(df, dislist='all', which='single', **kwargs):
    """
    Main plotting dispatcher function.

    Args:
        df (pd.DataFrame): Simulation results dataframe
        dislist (str): Disease type ('hiv', 'all')
        which (str): Plot type ('single' for individual runs, 'multi' for quantiles)
        **kwargs: Additional arguments passed to plotting functions

    Returns:
        matplotlib.figure.Figure: Generated figure
    """
    if dislist == 'hiv':
        fig = plot_hiv_sims(df, which=which, **kwargs)
    elif dislist == 'all' and which == 'single':
        fig = plot_coinfection(df, **kwargs)
    elif dislist == 'all' and which == 'multi':
        fig = plot_coinfection_quantiles(df, **kwargs)
    else:
        raise ValueError(f"Invalid combination: dislist='{dislist}', which='{which}'")

    return fig


def plot_calibrations(dislist='hiv', **kwargs):
    """
    Plot calibration results and print posterior parameter summaries.

    Args:
        dislist (str): Disease type to plot ('hiv' or 'all')
    """
    # Load calibration results
    df_filename = f'{RESULTS_DIR}/{LOCATION}_calib_stats_{dislist}.df'
    par_filename = f'{RESULTS_DIR}/{LOCATION}_par_stats_{dislist}.df'
    df_stats = sc.loadobj(df_filename)
    par_stats = sc.loadobj(par_filename)

    # Plot settings
    plot_kwargs = dict(
        start_year=1985,
        end_year=2025,
        which='multi',
        percentile_pairs=percentile_pairs,
        title=f'{dislist}_calib',
    )
    plot_kwargs.update(kwargs)

    # Generate plots
    if dislist == 'hiv':
        plot_hiv_sims(df_stats, **plot_kwargs)
    elif dislist == 'all':
        plot_coinfection_quantiles(df_stats, **plot_kwargs)

    # Print posterior parameter summaries
    print(f'\n{dislist.upper()} Calibration - Posterior Parameter Estimates:')
    print('=' * 60)
    pars = [p for p in par_stats.columns if p not in ['index', 'mismatch']]
    for p in pars:
        mean = par_stats[p]['mean']
        ci_low = par_stats[p]['5%']
        ci_high = par_stats[p]['95%']
        print(f'{p:30s}: {mean:6.3f} (95% CI: {ci_low:6.3f}–{ci_high:6.3f})')
    print('=' * 60)

    return


def plot_coinfection(df, location=LOCATION, start_year=2000, end_year=2040,
                     which='single', percentile_pairs=[[.1, .99]],
                     title='syph_coinf_plots', alpha=0.7):
    """
    Plot HIV-syphilis coinfection dynamics.

    Args:
        df (pd.DataFrame): Simulation results (single run or multi-run with quantiles)
        location (str): Location name for data files
        start_year (int): First year to plot
        end_year (int): Last year to plot
        which (str): 'single' for individual runs, 'multi' for quantiles
        percentile_pairs (list): List of [low, high] percentile pairs (for 'multi')
        title (str): Plot title/filename
        alpha (float): Line transparency (for 'single')

    Returns:
        matplotlib.figure.Figure: Generated figure
    """
    set_font(size=20)
    fig, axes = pl.subplots(3, 4, figsize=(18, 15))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    # Load data
    syph_data = pd.read_csv(f'{DATA_DIR}/{location}_syph_data.csv')
    hiv_data = pd.read_csv(f'{DATA_DIR}/{location}_hiv_data.csv')
    projections_data = pd.read_csv(f'{DATA_DIR}/{location}_projections.csv')
    gbd_estimates_new = pd.read_csv(f'{DATA_DIR}/{location}_gbd_estimates_new.csv')

    # Subset data to plotting years
    syph_data = syph_data.loc[(syph_data.year >= start_year) & (syph_data.year <= end_year)]
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]
    projections_data = projections_data.loc[(projections_data.year >= start_year) & (projections_data.year <= end_year)]
    gbd_estimates_new = gbd_estimates_new.loc[(gbd_estimates_new.year >= start_year) & (gbd_estimates_new.year <= end_year)]

    # Subset model results
    if which == 'single':
        dfplot = df.loc[df.timevec >= start_year]
        x = dfplot['timevec']
    else:  # multi
        dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]
        x = np.unique(dfplot.index)

    pn = 0

    # Panel 1: Population size
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['n_alive'], color='k', label='UNAIDS')
    if which == 'multi':
        ax.scatter(gbd_estimates_new.year, gbd_estimates_new['n_alive'],
                   color='darkviolet', alpha=0.6, label='GBD')
    resname = 'n_alive'
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x if which == 'single' else x[:-1],
                    y if which == 'single' else y[:-1],
                    label='Model')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                           facecolor=line.get_color())
    ax.set_title('Population size')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=10)
    sc.SIticks(ax)
    pn += 1

    # Panel 2: Syphilis prevalence
    ax = axes[pn]
    # ANC surveillance data
    syph_prev_data_early = np.array([0.045, 0.028, 0.035])
    syph_prev_time_early = np.array([2000, 2006, 2008])
    syph_prev_data_recent = np.array([0.0431, 0.0167, 0.0185, 0.0225, 0.0214,
                                      0.0190, 0.0237, 0.0193, 0.0251, 0.0229,
                                      0.0200, 0.0202, 0.0181])
    syph_prev_time_recent = np.array([2010, 2011, 2012, 2013, 2014, 2015,
                                      2016, 2017, 2018, 2019, 2020, 2021, 2022])
    ax.scatter(syph_prev_time_early, syph_prev_data_early * 100,
               color='k', marker='*', label='Early data')
    ax.scatter(syph_prev_time_recent, syph_prev_data_recent * 100,
               color='k', marker='d', label='Recent data')

    # Model results
    if which == 'single':
        ax.plot(x, dfplot['syph.active_prevalence'] * 100,
                label='Active', alpha=alpha)
        ax.plot(x, dfplot['syph.detected_pregnant_prevalence'] * 100,
                label='ANC', alpha=alpha)
    else:  # multi
        resnames = {'Active': 'syphilis_active_prevalence',
                    'ANC': 'syphilis_detected_pregnant_prevalence'}
        for rlabel, rname in resnames.items():
            scale = 100 * 0.7  # Scaling factor
            y = dfplot[(rname, '50%')]
            line, = ax.plot(x, y * scale, label=rlabel)
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl * scale, yu * scale, alpha=alphas[idx],
                               facecolor=line.get_color())
    ax.legend(frameon=False, loc='upper right', fontsize=10)
    ax.set_title('Syphilis prevalence (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Panel 3: Syphilis prevalence by HIV status
    ax = axes[pn]
    if which == 'single':
        ax.plot(x, dfplot['coinfection_stats.syph_prev_no_hiv'] * 100,
                label='HIV−', alpha=alpha)
        ax.plot(x, dfplot['coinfection_stats.syph_prev_has_hiv'] * 100,
                label='HIV+', alpha=alpha)
    else:  # multi
        resnames = {'HIV−': 'coinfection_stats_syph_prev_no_hiv',
                    'HIV+': 'coinfection_stats_syph_prev_has_hiv'}
        for rlabel, rname in resnames.items():
            y = dfplot[(rname, '50%')]
            line, = ax.plot(x, y * 100, label=rlabel)
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx],
                               facecolor=line.get_color())
    ax.legend(frameon=False)
    ax.set_title('Syphilis prevalence\nby HIV status (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Panel 4: Syphilis infections
    ax = axes[pn]
    if which == 'single':
        ax = plot_single(ax, syph_data, dfplot, 'syph.new_infections',
                        'syph.new_infections', annualize=False)
    else:  # multi
        resname = 'syphilis_new_infections'
        ax.scatter(syph_data.year, syph_data['syphilis.new_infections'],
                   label='Data', color='k')
        ax.scatter(gbd_estimates_new.year, gbd_estimates_new['syphilis.new_infections'],
                   color='darkviolet', alpha=0.6, label='GBD')
        y = dfplot[(resname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label='Model')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                           facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=10)
    ax.set_title('Syphilis infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 5: Cumulative congenital syphilis cases
    ax = axes[pn]
    if which == 'single':
        ax = plot_single(ax, syph_data, dfplot, 'syph.new_congenital',
                        'syph.new_congenital', annualize=False, smooth=True)
        ax.set_title('Congenital syphilis cases')
    else:  # multi
        resname = 'syphilis_cum_congenital'
        ydata = syph_data['syphilis.cum_congenital'] - syph_data['syphilis.cum_congenital'].iloc[0]
        ax.scatter(syph_data.year, ydata, label='Data', color='k')
        y = dfplot[(resname, '50%')].values
        y = y - y[0]
        line, = ax.plot(x, y, label='Model')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")].values - y[0]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")].values - y[0]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=10)
        ax.set_title(f'Cumulative CS cases, {start_year}–')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 6: Cumulative congenital syphilis deaths
    ax = axes[pn]
    if which == 'single':
        ax = plot_single(ax, syph_data, dfplot, 'syph.new_congenital_deaths',
                        'syph.new_congenital_deaths', annualize=False, smooth=True)
        ax.set_title('Congenital syphilis deaths')
    else:  # multi
        resname = 'syphilis_cum_congenital_deaths'
        ydata = syph_data['syphilis.cum_congenital_deaths'] - syph_data['syphilis.cum_congenital_deaths'].iloc[0]
        ax.scatter(syph_data.year, ydata, label='Data', color='k')
        y = dfplot[(resname, '50%')].values
        y = y - y[0]
        line, = ax.plot(x, y, label='Model')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")].values - y[0]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")].values - y[0]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=10)
        ax.set_title(f'Cumulative CS deaths, {start_year}–')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 7: Syphilis treatments
    ax = axes[pn]
    if which == 'single':
        # For single sim, would need appropriate data column
        ax.text(0.5, 0.5, 'Treatments\n(multi-sim only)',
                ha='center', va='center', transform=ax.transAxes)
    else:  # multi
        resname = 'syphilis_new_treated'
        y = dfplot[(resname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label='Treatments')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                           facecolor=line.get_color())
    ax.set_title('Syphilis treatments')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 8: Unnecessary syphilis treatments
    ax = axes[pn]
    if which == 'single':
        # For single sim, would need appropriate data column
        ax.text(0.5, 0.5, 'Overtreatment\n(multi-sim only)',
                ha='center', va='center', transform=ax.transAxes)
    else:  # multi
        resname = 'syphilis_new_treated_unnecessary'
        y = dfplot[(resname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label='Overtreatment')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                           facecolor=line.get_color())
    ax.set_title('Unnecessary syphilis treatments')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 9: HIV infections
    ax = axes[pn]
    if which == 'single':
        ax = plot_single(ax, hiv_data, dfplot, 'hiv.new_infections',
                        'hiv.new_infections', annualize=False)
    else:  # multi
        resname = 'hiv_new_infections'
        ax.scatter(hiv_data.year, hiv_data['hiv.new_infections'],
                   label='UNAIDS', color='k')
        ax.scatter(projections_data.year, projections_data['hiv.new_infections'],
                   label='GBD (old)', alpha=0.3, color='k')
        ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.new_infections'],
                   color='darkviolet', alpha=0.6, label='GBD (new)')
        y = dfplot[(resname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label='Model')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                           facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=8)
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 10: HIV deaths
    ax = axes[pn]
    if which == 'single':
        ax = plot_single(ax, hiv_data, dfplot, 'hiv.new_deaths',
                        'hiv.new_deaths', annualize=False)
    else:  # multi
        resname = 'hiv_new_deaths'
        ax.scatter(hiv_data.year, hiv_data['hiv.new_deaths'],
                   label='UNAIDS', color='k')
        ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.new_deaths'],
                   color='darkviolet', alpha=0.6, label='GBD')
        y = dfplot[(resname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label='Model')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                           facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=10)
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 11: People living with HIV (total, diagnosed, treated)
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv.n_infected'], color='k')
    if which == 'multi':
        ax.scatter(projections_data.year, projections_data['hiv.n_infected'],
                   alpha=0.3, color='k')
        ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.n_infected'],
                   color='darkviolet', alpha=0.6)

    if which == 'single':
        ax.plot(x, dfplot['hiv.n_infected'], label='PLHIV', alpha=alpha)
        ax.plot(x, dfplot['hiv.n_diagnosed'], label='Diagnosed', alpha=alpha)
        ax.plot(x, dfplot['hiv.n_on_art'], label='On ART', alpha=alpha)
    else:  # multi
        resnames = {'Total': 'hiv_n_infected', 'Diagnosed': 'hiv_n_diagnosed',
                    'On ART': 'hiv_n_on_art'}
        for rlabel, rname in resnames.items():
            y = dfplot[(rname, '50%')]
            line, = ax.plot(x[:-1], y[:-1], label=rlabel)
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx],
                               facecolor=line.get_color())
    ax.set_title('PLHIV: total, diagnosed, on ART')
    ax.legend(frameon=False, fontsize=10)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 12: HIV prevalence
    ax = axes[pn]
    if which == 'single':
        ax.scatter(hiv_data.time, hiv_data['hiv.prevalence_15_49'] * 100,
                   color='k', label='UNAIDS')
        ax.plot(x, dfplot['hiv.prevalence_15_49'] * 100,
                label='Model', alpha=alpha)
    else:  # multi
        resname = 'hiv_prevalence'
        ax.scatter(hiv_data.year, hiv_data['hiv.prevalence'] * 100,
                   label='UNAIDS', color='k')
        ax.scatter(projections_data.year, projections_data['hiv.prevalence'] * 100,
                   label='GBD (old)', color='k', alpha=0.3)
        ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.prevalence'] * 100,
                   label='GBD (new)', color='darkviolet', alpha=0.6)
        y = dfplot[(resname, '50%')]
        line, = ax.plot(x, y * 100, label='Model')
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx],
                           facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=8)
    ax.set_title('HIV prevalence (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Finalize
    sc.figlayout()
    sc.savefig(f'{FIGURES_DIR}/{title}_{which}.png', dpi=100)

    return fig


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='hiv_plots'):
    """ Create quantile or individual plots of HIV epi dynamics """
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(18, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    hiv_data = pd.read_csv(f'{DATA_DIR}/{LOCATION}_hiv_data.csv')
    hiv_data = hiv_data.loc[(hiv_data.time >= start_year) & (hiv_data.time <= end_year)]
    dfplot = df.loc[(df.index >= start_year) & (df.index <= end_year)]

    pn = 0
    x = dfplot.index

    # Population size
    ax = axes[pn]
    resname = 'n_alive'
    ax.scatter(hiv_data.time, hiv_data[resname], color='k', label='Data')
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='Modeled')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Population size')
    ax.legend(frameon=False)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    pn += 1

    # PLHIV
    ax = axes[pn]
    resname = 'hiv.n_infected'
    ax.scatter(hiv_data.time, hiv_data[resname], label='Data', color='k')
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='PLHIV')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    resname = 'hiv.prevalence_15_49'
    ax.scatter(hiv_data.time, hiv_data[resname] * 100, label='Data', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y*100, label='Prevalence')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Infections
    ax = axes[pn]
    resname = 'hiv.new_infections'
    ax.scatter(hiv_data.time, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='HIV infections')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    resname = 'hiv.new_deaths'
    ax.scatter(hiv_data.time, hiv_data[resname], label='UNAIDS', color='k')
    x = dfplot.index
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='HIV deaths')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # 90-90-90
    ax = axes[pn]
    ax.scatter(hiv_data.time, hiv_data['hiv.n_infected'], color='k')  # label='UNAIDS',
    resnames = {'PLHIV': 'hiv.n_infected', 'Dx': 'hiv.n_diagnosed', 'Treated': 'hiv.n_on_art'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y, label=rlabel)
        # if which == 'multi':
        #     for idx, percentile_pair in enumerate(percentile_pairs):
        #         yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
        #         yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
        #         ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Diagnosed and treated')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig(f"{FIGURES_DIR}/" + title + str(start_year) + "_" + which + ".png", dpi=100)

    return fig

