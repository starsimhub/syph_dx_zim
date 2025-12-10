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
        fig = plot_coinfection(df, which=which, **kwargs)
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
        start_year=1990,
        end_year=2025,
        which='multi',
        percentile_pairs=percentile_pairs,
        title=f'{dislist}_calib',
    )
    plot_kwargs.update(kwargs)

    # Generate plots
    if dislist == 'hiv':
        plot_hiv_sims(df_stats, **plot_kwargs)
    elif dislist == 'syph':
        plot_coinfection(df_stats, **plot_kwargs)

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

    # Subset data to plotting years
    syph_data = syph_data.loc[(syph_data.time >= start_year) & (syph_data.time <= end_year)]
    hiv_data = hiv_data.loc[(hiv_data.time >= start_year) & (hiv_data.time <= end_year)]

    # Subset model results
    dfplot = df.loc[(df.timevec >= start_year) & (df.timevec <= end_year)]
    x = dfplot['timevec']
    # if which == 'single':
    #     dfplot = df.loc[(df.timevec >= start_year) & (df.timevec <= end_year)]
    #     x = dfplot['timevec']
    # else:  # multi
    #     dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]
    #     x = np.unique(dfplot.index)

    pn = 0

    # Panel 1: Population size
    ax = axes[pn]
    ax.scatter(hiv_data.time, hiv_data['n_alive'], color='k', label='UNAIDS')
    resname = 'n_alive'
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y, label='Model')
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
    resnames = {'Active': 'syph.active_prevalence',
                'ANC': 'syph.detected_pregnant_prevalence'}
    for rlabel, rname in resnames.items():
        scale = 100  # to convert to percentage
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y * scale, label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl * scale, yu * scale, alpha=alphas[idx], facecolor=line.get_color())
    ax.legend(frameon=False, loc='upper right', fontsize=10)
    ax.set_title('Syphilis prevalence (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Panel 3: Syphilis prevalence by HIV status
    ax = axes[pn]
    resnames = {'HIV−': 'coinfection_stats.syph_prev_no_hiv',
                'HIV+': 'coinfection_stats.syph_prev_has_hiv'}
    for rlabel, rname in resnames.items():
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x, y * 100, label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.legend(frameon=False)
    ax.set_title('Syphilis prevalence\nby HIV status (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Panel 4-5: Syphilis burden and new infections
    resnames = ['syph.n_infected', 'syph.new_infections']
    for resname in resnames:
        ax = axes[pn]
        ax.scatter(syph_data.time, syph_data[resname], label='GBD', color='k')
        y = get_y(dfplot, which, resname)
        line, = ax.plot(x[:-1], y[:-1], label='Model')
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=10)
        subtitle = 'Syphilis burden' if resname == 'syph.n_infected' else 'New syphilis infections'
        ax.set_title(subtitle)
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)
        pn += 1

    # Panel 6-7: Cumulative congenital syphilis cases
    for resname in ['syph.new_congenital', 'syph.new_congenital_deaths']:
        ax = axes[pn]
        ydata = syph_data[resname]
        ax.scatter(syph_data.time, ydata, label='Data', color='k')
        y = get_y(dfplot, which, resname).values
        y = y - y[0]
        line, = ax.plot(x, y, label='Model')
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")].values
                yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")].values
                ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=10)
        subtitle = 'CS cases' if resname == 'syph.new_congenital' else 'CS deaths'
        ax.set_title(f'{subtitle}, {start_year}–')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)
        pn += 1

    # Panel 8: Syphilis treatments
    for resname in ['syph.new_treated']:  #, 'syph.new_treated_unnecessary']:
        ax = axes[pn]
        y = get_y(dfplot, which, resname)
        line, = ax.plot(x[:-1], y[:-1], label='Treatments')
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
        ax.set_title('Syphilis treatments')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)
        pn += 1

    # Panel 9-10: HIV infections and deaths
    for resname in ['hiv.new_infections', 'hiv.new_deaths']:
        ax = axes[pn]
        ax.scatter(hiv_data.time, hiv_data[resname], label='UNAIDS', color='k')
        y = get_y(dfplot, which, resname)
        line, = ax.plot(x[:-1], y[:-1], label='Model')
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
        ax.legend(frameon=False, fontsize=8)
        ax.set_title('HIV infections')
        ax.set_ylim(bottom=0)
        sc.SIticks(ax)
        pn += 1

    # Panel 11: People living with HIV (total, diagnosed, treated)
    ax = axes[pn]
    ax.scatter(hiv_data.time, hiv_data['hiv.n_infected'], color='k')

    resnames = {'Total': 'hiv.n_infected', 'Diagnosed': 'hiv.n_diagnosed', 'On ART': 'hiv.n_on_art'}
    for rlabel, rname in resnames.items():
        y = get_y(dfplot, which, rname)
        line, = ax.plot(x[:-1], y[:-1], label=rlabel)
        if which == 'multi':
            for idx, percentile_pair in enumerate(percentile_pairs):
                yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
                yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
                ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV: total, diagnosed, on ART')
    ax.legend(frameon=False, fontsize=10)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax)
    pn += 1

    # Panel 12: HIV prevalence
    ax = axes[pn]
    resname = 'hiv.prevalence_15_49'
    ax.scatter(hiv_data.time, hiv_data[resname] * 100, color='k', label='UNAIDS')
    y = get_y(dfplot, which, resname)
    line, = ax.plot(x, y * 100, label='Model')
    if which == 'multi':
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
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

