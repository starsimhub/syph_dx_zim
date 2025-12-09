# %% Imports and settings
import sciris as sc
import pylab as pl
import numpy as np
import pandas as pd
from utils import set_font, get_y, plot_single

location = 'zimbabwe'


def plot_sims(df, dislist='all', which='single', **kwargs):
    if dislist == 'hiv':
        fig = plot_hiv_sims(df, which=which, **kwargs)
    elif dislist == 'all' and which=='single':
        fig = plot_coinfection(df, **kwargs)
    elif dislist == 'all' and which=='multi':
        fig = plot_coinfection_quantiles(df, **kwargs)
    else:
        raise ValueError(f"Unknown disease type: {dislist}")

    return fig


def plot_coinfection(sim_output, location='zimbabwe', start_year=1970, title='syph_coinf_plots', alpha=0.7):
    """ Create plots for syphilis and HIV"""
    set_font(size=14)
    fig, axes = pl.subplots(2, 5, figsize=(18, 7))
    axes = axes.ravel()

    syph_data = pd.read_csv(f'data/{location}_syph_data.csv')
    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')

    sim_output = sim_output.loc[sim_output.timevec >= start_year]
    x = sim_output['timevec']
    pn = 0

    # Population size
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['n_alive'], color='k', label='UNAIDS')
    y0 = sim_output['n_alive']
    ax.plot(x, y0, label='Modeled', alpha=alpha)
    ax.set_title('Population size')
    ax.legend(frameon=False)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    pn += 1

    # Syphilis prevalence at ANCs and with symptoms
    ax = axes[pn]
    syph_prev_data1 = np.array([0.045, 0.028, 0.035])
    syph_prev_time1 = np.array([2000, 2006, 2008])
    syph_prev_data2 = np.array([0.0431, 0.0167, 0.0185, 0.0225, 0.0214, 0.0190, 0.0237, 0.0193, 0.0251, 0.0229, 0.0200, 0.0202, 0.0181])
    syph_prev_time2 = np.array([2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022])
    ax.scatter(syph_prev_time1, syph_prev_data1 * 100, color='k', marker='d')  #, label='Korenromp 2017')
    ax.scatter(syph_prev_time2, syph_prev_data2 * 100, color='k', marker='d')  #, label='JR')
    y0 = sim_output['syph.active_prevalence']
    y1 = sim_output['syph.detected_pregnant_prevalence']
    y2 = sim_output['syph.delivery_prevalence']
    ax.plot(x, y0 * 100, label='Active', alpha=alpha)
    ax.plot(x, y1 * 100, label='Pregnant', alpha=alpha)
    ax.plot(x, y2 * 100, label='Delivery', alpha=alpha)
    ax.set_title('Syphilis prevalence (%)')
    ax.legend(frameon=False)
    pn += 1

    # Adult syphilis prevalence by HIV status
    ax = axes[pn]
    y1 = sim_output['coinfection_stats.syph_prev_no_hiv']
    y2 = sim_output['coinfection_stats.syph_prev_has_hiv']
    ax.plot(x, y1 * 100, label='HIV-', alpha=alpha)
    ax.plot(x, y2 * 100, label='HIV+', alpha=alpha)
    ax.set_title('Syphilis by HIV')
    ax.legend(frameon=False)
    pn += 1

    # Syphilis infections
    ax = axes[pn]
    ax = plot_single(ax, syph_data, sim_output, 'syph.new_infections', 'syph.new_infections', annualize=False)
    ax.set_title('New syphilis infections')
    pn += 1

    # Congenital outcomes
    ax = axes[pn]
    ax = plot_single(ax, syph_data, sim_output, 'syph.new_congenital', 'syph.new_congenital', annualize=False, smooth=True)
    ax.set_title('Congenital cases')
    pn += 1

    # Congenital outcomes
    ax = axes[pn]
    ax = plot_single(ax, syph_data, sim_output, 'syph.new_congenital_deaths', 'syph_new.congenital_deaths', annualize=False, smooth=True)
    ax.set_title('Congenital deaths')
    pn += 1

    # PLHIV
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv.n_infected'], label='UNAIDS estimates', color='k')
    y0 = sim_output['hiv.n_infected']
    y1 = sim_output['hiv.n_diagnosed']
    y2 = sim_output['hiv.n_on_art']

    ax.plot(x, y0, label='Modeled PLHIV', alpha=alpha)
    ax.plot(x, y1, label='Diagnosed PLHIV', alpha=alpha)
    ax.plot(x, y2, label='Treated PLHIV', alpha=alpha)
    ax.set_title('PLHIV: total, diagnosed, and treated')
    ax.legend(frameon=False)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv.prevalence'] * 100, color='k', label='UNAIDS')
    y0 = sim_output['hiv.prevalence']
    ax.plot(x, y0 * 100, label='Modeled', alpha=alpha)
    ax.set_title('HIV prevalence 15-49 (%)')
    ax.legend(frameon=False)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    pn += 1

    # HIV infections
    ax = axes[pn]
    ax = plot_single(ax, hiv_data, sim_output, 'hiv.new_infections', 'hiv.new_infections', annualize=False)
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    ax = plot_single(ax, hiv_data, sim_output, 'hiv.new_deaths', 'hiv.new_deaths', annualize=False)
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + str(start_year) + ".png", dpi=100)


    return fig


def plot_coinfection_quantiles(df, location='zimbabwe', start_year=2000, end_year=2040, percentile_pairs=[[.1, .99]], title='syph_coinf_plots'):
    """ Create quantile plots for syphilis and HIV"""
    set_font(size=20)
    fig, axes = pl.subplots(3, 4, figsize=(18, 15))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    syph_data = pd.read_csv(f'data/{location}_syph_data.csv')
    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    projections_data = pd.read_csv(f'data/{location}_projections.csv')
    gbd_estimates_new = pd.read_csv(f'data/{location}_gbd_estimates_new.csv')

    # Subset data and model results for particular years
    syph_data = syph_data.loc[(syph_data.year >= start_year) & (syph_data.year <= end_year)]
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]
    projections_data = projections_data.loc[(projections_data.year >= start_year) & (projections_data.year <= end_year)]
    gbd_estimates_new = gbd_estimates_new.loc[(gbd_estimates_new.year >= start_year) & (gbd_estimates_new.year <= end_year)]
    # df['year'] = np.floor(np.round(df.index, 1)).astype(int)
    dfplot = df.iloc[(df.index >= start_year) & (df.index <= end_year)]

    # Population size
    pn = 0
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['n_alive'], color='k')
    ax.scatter(gbd_estimates_new.year, gbd_estimates_new['n_alive'], color='darkviolet', alpha=0.6)
    resname = 'n_alive'
    x = np.unique(dfplot.index)
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1])
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color(), label=f"{percentile_pair[0]:.0%}" + ' - ' + f"{percentile_pair[1]:.0%}")
    ax.set_title('Population size')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False)
    sc.SIticks(ax=ax)
    pn += 1

    # Adult syphilis prevalence
    ax = axes[pn]
    resnames = {'Active': 'syphilis_active_prevalence', 'ANC': 'syphilis_detected_pregnant_prevalence'}
    syph_prev_data1 = np.array([0.045, 0.028, 0.035])
    syph_prev_time1 = np.array([2000, 2006, 2008])
    syph_prev_data2 = np.array([0.0431, 0.0167, 0.0185, 0.0225, 0.0214, 0.0190, 0.0237, 0.0193, 0.0251, 0.0229, 0.0200, 0.0202, 0.0181])
    syph_prev_time2 = np.array([2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022])
    ax.scatter(syph_prev_time1, syph_prev_data1 * 100, color='k', marker='*')  #, label='Korenromp 2017')
    ax.scatter(syph_prev_time2, syph_prev_data2 * 100, color='k', marker='d')  #, label='JR')
    for rlabel, rname in resnames.items():
        scale = 100*.7
        y = dfplot[(rname, '50%')]
        line, = ax.plot(x, y * scale, label=rlabel)
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * scale, yu * scale, alpha=alphas[idx], facecolor=line.get_color())
    ax.legend(frameon=False, loc='upper right')
    ax.set_title('Syphilis prevalence (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Adult syphilis prevalence by HIV status
    ax = axes[pn]
    resnames = {'HIV-': 'coinfection_stats_syph_prev_no_hiv', 'HIV+': 'coinfection_stats_syph_prev_has_hiv'}
    for rlabel, rname in resnames.items():
        y = dfplot[(rname, '50%')]
        line, = ax.plot(x, y * 100, label=rlabel)
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.legend(frameon=False)
    ax.set_title('Syphilis prevalence\n by HIV status (%)')
    ax.set_ylim(bottom=0)
    pn += 1

    # Syphilis infections
    ax = axes[pn]
    resname = 'syphilis_new_infections'
    ax.scatter(syph_data.year, syph_data['syphilis.new_infections'], label='Data', color='k')
    ax.scatter(gbd_estimates_new.year, gbd_estimates_new['syphilis.new_infections'], color='darkviolet', alpha=0.6)
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Syphilis infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Congenital outcomes
    ax = axes[pn]
    resname = 'syphilis_cum_congenital'
    ydata = syph_data['syphilis.cum_congenital'] - syph_data['syphilis.cum_congenital'].iloc[0]
    ax.scatter(syph_data.year, ydata, label='Data', color='k')
    x = dfplot.index
    y = dfplot[(resname, '50%')].values
    y = y - y[0]
    line, = ax.plot(x, y, label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")].values
        yl = yl - yl[0]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")].values
        yu = yu - yu[0]
        ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title(f'Total CS cases, {start_year}–')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Congenital deaths
    ax = axes[pn]
    resname = 'syphilis_cum_congenital_deaths'
    ydata = syph_data['syphilis.cum_congenital_deaths'] - syph_data['syphilis.cum_congenital_deaths'].iloc[0]
    ax.scatter(syph_data.year, ydata, label='Data', color='k')
    x = dfplot.index
    y = dfplot[(resname, '50%')].values
    y = y - y[0]
    line, = ax.plot(x, y, label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")].values
        yl = yl - yl[0]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")].values
        yu = yu - yu[0]
        ax.fill_between(x, yl, yu, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title(f'Total CS deaths, {start_year}–')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Treatment
    ax = axes[pn]
    resname = 'syphilis_new_treated'
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Syphilis treatments')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # Treatment
    ax = axes[pn]
    resname = 'syphilis_new_treated_unnecessary'
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('Syphilis false negatives')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV infections
    ax = axes[pn]
    resname = 'hiv_new_infections'
    ax.scatter(hiv_data.year, hiv_data['hiv.new_infections'], label='UNAIDS', color='k')
    ax.scatter(projections_data.year, projections_data['hiv.new_infections'], label='GBD', alpha=0.3, color='k')
    ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.new_infections'], color='darkviolet', alpha=0.6)
    # ax.scatter(hiv_optima_data.year, hiv_optima_data[resname], label='Optima Data', color='tab:red')
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV infections')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV deaths
    ax = axes[pn]
    resname = 'hiv_new_deaths'
    ax.scatter(hiv_data.year, hiv_data['hiv.new_deaths'], label='UNAIDS', color='k')
    ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.new_deaths'], color='darkviolet', alpha=0.6)
    # ax.scatter(hiv_optima_data.year, hiv_optima_data[resname], label='Optima Data', color='tab:red')
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x[:-1], y[:-1], label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV-related deaths')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # PLHIV
    ax = axes[pn]
    ax.scatter(hiv_data.year, hiv_data['hiv.n_infected'], color='k')  # label='UNAIDS',
    ax.scatter(projections_data.year, projections_data['hiv.n_infected'], alpha=0.3, color='k')  # label='GBD',
    ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.n_infected'], color='darkviolet', alpha=0.6)
    # ax.scatter(hiv_optima_data.year, hiv_optima_data['hiv.n_infected'], label='Optima Data', color='tab:red')
    resnames = {'Total': 'hiv_n_infected', 'Dx': 'hiv_n_diagnosed', 'Treated': 'hiv_n_on_art'}
    for rlabel, rname in resnames.items():
        x = dfplot.index
        y = dfplot[(rname, '50%')]
        line, = ax.plot(x[:-1], y[:-1], label=rlabel)
        for idx, percentile_pair in enumerate(percentile_pairs):
            yl = dfplot[(rname, f"{percentile_pair[0]:.0%}")]
            yu = dfplot[(rname, f"{percentile_pair[1]:.0%}")]
            ax.fill_between(x[:-1], yl[:-1], yu[:-1], alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('PLHIV: total, diagnosed, treated')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax=ax)
    pn += 1

    # HIV prevalence
    ax = axes[pn]
    resname = 'hiv_prevalence'
    ax.scatter(hiv_data.year, hiv_data['hiv.prevalence'] * 100, label='Data', color='k')
    ax.scatter(projections_data.year, projections_data['hiv.prevalence'] * 100, label='GBD', color='k', alpha=0.3)
    ax.scatter(gbd_estimates_new.year, gbd_estimates_new['hiv.prevalence'] * 100, label='GBD (new)', color='darkviolet', alpha=0.6)
    # ax.scatter(hiv_optima_data.year, hiv_optima_data[resname] * 100, label='Optima Data', color='tab:red')
    x = dfplot.index
    y = dfplot[(resname, '50%')]
    line, = ax.plot(x, y * 100, label='Total')
    for idx, percentile_pair in enumerate(percentile_pairs):
        yl = dfplot[(resname, f"{percentile_pair[0]:.0%}")]
        yu = dfplot[(resname, f"{percentile_pair[1]:.0%}")]
        ax.fill_between(x, yl * 100, yu * 100, alpha=alphas[idx], facecolor=line.get_color())
    ax.set_title('HIV prevalence (%)')
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    pn += 1

    sc.figlayout()
    sc.savefig("figures/" + title + ".png", dpi=100)

    return fig


def make_calibration_plot(modeldf=None, to_plot=None, location='zimbabwe'):
    """ Make plot"""
    set_font(size=14)
    fig, axes = pl.subplots(2, 4, figsize=(15, 7))
    axes = axes.ravel()
    data = pd.read_csv(f'data/{location}_data.csv')

    for pn, resname in enumerate(to_plot):
        modeldf_stats = modeldf.groupby(['year']).describe(percentiles=[.1, .5, .9])
        x = modeldf_stats.index
        medians = modeldf_stats[(resname, '50%')]
        medians.name = resname
        quartiles1 = modeldf_stats[(resname, '10%')]
        quartiles3 = modeldf_stats[(resname, '90%')]

        ax = axes[pn]
        sns.lineplot(data=medians, ax=ax)
        ax.fill_between(x, quartiles1, quartiles3, alpha=0.3)
        ax.scatter(data.year, data[resname], color='k')
        ax.set_title(resname)
        ax.set_ylim(bottom=0)
        if 'prevalence' not in resname: sc.SIticks(ax)

    sc.figlayout()
    sc.savefig("figures/syph_calib_plots.png", dpi=100)


def plot_hiv_sims(df, start_year=2000, end_year=2025, which='single', percentile_pairs=[[.1, .99]], title='hiv_plots'):
    """ Create quantile or individual plots of HIV epi dynamics """
    set_font(size=20)
    fig, axes = pl.subplots(2, 3, figsize=(18, 7))
    axes = axes.ravel()
    alphas = np.linspace(0.2, 0.5, len(percentile_pairs))

    hiv_data = pd.read_csv(f'data/{location}_hiv_data.csv')
    hiv_data = hiv_data.loc[(hiv_data.year >= start_year) & (hiv_data.year <= end_year)]
    dfplot = df.loc[(df.index >= start_year) & (df.index <= end_year)]

    pn = 0
    x = dfplot.index

    # Population size
    ax = axes[pn]
    resname = 'n_alive'
    ax.scatter(hiv_data.year, hiv_data[resname], color='k', label='Data')
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
    ax.scatter(hiv_data.year, hiv_data[resname], label='Data', color='k')
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
    ax.scatter(hiv_data.year, hiv_data[resname] * 100, label='Data', color='k')
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
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
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
    ax.scatter(hiv_data.year, hiv_data[resname], label='UNAIDS', color='k')
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
    ax.scatter(hiv_data.year, hiv_data['hiv.n_infected'], color='k')  # label='UNAIDS',
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
    sc.savefig("figures/" + title + str(start_year) + "_" + which + ".png", dpi=100)

    return fig

