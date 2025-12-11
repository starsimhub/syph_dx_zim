"""
Plot syphilis and HIV epidemiology in Zimbabwe: Fig 1 in manuscript.

If updates are needed:
    1. Make required updates to the model (model.py) or data folder
    2. Run run_calibration to generate the files 'results/zimbabwe_all_calib_stats.df'
    3. Run run_plot_data.py to generate the epi result files:
            epi_df = 'results/epi_df.df'
            sw_df = 'results/sw_df.df'
            coinf_df = 'results/coinf_df.df'
    4. Run this script to generate the figure
"""

# Import packages
import sciris as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font


# %% Plotting functions
def plot_syph_prev_by_hiv(coinf_df, ax=None, start_year=2000):
    """Plot time series of active syphilis prevalence by HIV status"""
    set_font(size=20)
    colors = ['#ee7989', '#4682b4']

    # Subset to start year
    bi = sc.findfirst(coinf_df.index, start_year)
    x = coinf_df.index[bi:]

    # Plot syphilis prevalence stratified by HIV status
    y_hiv_neg = coinf_df['syph_prev_no_hiv'][bi:] * 100
    y_hiv_pos = coinf_df['syph_prev_has_hiv'][bi:] * 100

    ax.plot(x, y_hiv_neg, color=colors[0], label='HIV−', linewidth=2)
    ax.plot(x, y_hiv_pos, color=colors[1], label='HIV+', linewidth=2)

    ax.legend(frameon=False)
    ax.set_ylabel('Prevalence (%)')
    ax.set_xlabel('Year')
    ax.set_title('Active syphilis prevalence by HIV status')
    ax.set_ylim(bottom=0)

    return ax


def plot_infections_by_sex(epi_df, ax=None, start_year=2000):
    """Plot time series of annual new infections by sex for syphilis and HIV"""
    set_font(size=20)
    colors_syph = ['#d46e9c', '#8b4789']  # Pink/purple for syphilis
    colors_hiv = ['#2f734a', '#1a4d2e']   # Green for HIV

    # Subset to start year
    bi = sc.findfirst(coinf_df.index, start_year)
    x = coinf_df.index[bi:]

    # Plot syphilis infections
    # syph_df = epi_df.loc[(epi_df.disease == 'syph')].copy()
    #     # sns.barplot(data=thisdf, x="age", y="new_infections", hue="sex", ax=ax, palette=scolors)
    #     thisdf['prevalence'] *= 100
    #     sns.barplot(data=thisdf, x="age", y="prevalence", hue="sex", ax=ax, palette=scolors)
    #
    # ax.plot(x, coinf_df['new_infections_f_syphilis'][bi:],
    #         color=colors_syph[0], label='Syphilis F', linewidth=2)
    # ax.plot(x, coinf_df['new_infections_m_syphilis'][bi:],
    #         color=colors_syph[1], label='Syphilis M', linewidth=2, linestyle='--')
    #
    # # Plot HIV infections
    # ax.plot(x, coinf_df['new_infections_f_hiv'][bi:],
    #         color=colors_hiv[0], label='HIV F', linewidth=2)
    # ax.plot(x, coinf_df['new_infections_m_hiv'][bi:],
    #         color=colors_hiv[1], label='HIV M', linewidth=2, linestyle='--')

    ax.legend(frameon=False, ncol=2, fontsize=14)
    ax.set_ylabel('Annual infections')
    ax.set_xlabel('Year')
    ax.set_title('New syphilis and HIV infections by sex')
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')

    return ax


def plot_prevalence_by_age(epi_df, disease='syphilis', ax=None):
    """Plot prevalence by age and sex"""
    set_font(size=20)
    scolors = ['#ee7989', '#4682b4']  # Pink for F, blue for M

    # Filter data for specified disease and exclude extreme age groups
    thisdf = epi_df.loc[(epi_df.disease == disease) &
                        (epi_df.age != '0-15') &
                        (epi_df.age != '65+')].copy()

    # Convert prevalence to percentage
    thisdf['prevalence'] *= 100

    # Create bar plot
    sns.barplot(data=thisdf, x="age", y="prevalence", hue="sex",
                ax=ax, palette=scolors)

    ax.set_title(f'{disease.capitalize()} prevalence by age')
    ax.set_ylabel('Prevalence (%)')
    ax.set_xlabel('Age group')
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, title='Sex')

    return ax


def plot_infections_by_sw(sw_df, disease=None, ax=None, start_year=2000, end_year=2019):
    """Plot infections acquired and transmitted by sex work status"""
    set_font(size=20)
    width = 0.6

    groups = {'fsw': 'FSW', 'client': 'Client', 'non_fsw': 'F', 'non_client': 'M'}
    colors = ['#d46e9c', '#2f734a', '#d46e9c', '#2f734a']
    alphas = [0.9, 0.9, 0.3, 0.3]

    si = sc.findfirst(sw_df.index, start_year)
    ei = sc.findfirst(sw_df.index, end_year)

    # New infections acquired by sex and sex work status
    bottom = np.zeros(2)
    x = np.array([0.5, 1.5])
    g = 0
    for group, glabel in groups.items():
        vals = np.array([
            sw_df[f'new_infections_{group}_{disease}'][si:ei].mean(),
            sw_df[f'new_transmissions_{group}_{disease}'][si:ei].mean(),
        ])
        p = ax.barh(x, vals, width, label=glabel, left=bottom,
                    color=colors[g], alpha=alphas[g])
        ax.bar_label(p, labels=[glabel, glabel], label_type='center')
        bottom += vals
        g += 1

    ax.set_title(f'{disease.upper()} infections\n{start_year}–{end_year} average')
    sc.SIticks(ax, axis='x')
    ax.set_xlim(left=0)
    ax.set_yticks(x, ['Acquired', 'Transmitted'])
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Calculate and print sex work share of transmission
    total_trans = sum(sw_df[f'new_transmissions_{group}_{disease}'][si:ei].mean()
                      for group in groups.keys())
    sw_trans = (sw_df[f'new_transmissions_fsw_{disease}'][si:ei].mean() +
                sw_df[f'new_transmissions_client_{disease}'][si:ei].mean())
    print(f'{disease.upper()} sex work share of transmission: {sw_trans/total_trans:.1%}')

    return ax


# %% Run as script
if __name__ == '__main__':

    # Load data files
    epi_df = sc.loadobj(f'results/epi_df.df')
    sw_df = sc.loadobj(f'results/sw_df.df')
    coinf_df = sc.loadobj(f'results/coinf_df.df')

    # Initialize plot - 2x3 grid
    set_font(size=20)
    fig = pl.figure(figsize=(20, 12))
    gs = GridSpec(2, 3, left=0.08, right=0.98, bottom=0.08, top=0.95,
                  wspace=0.25, hspace=0.3)

    # Panel A: Syphilis prevalence by HIV status (time series)
    ax = fig.add_subplot(gs[0, 0])
    plot_syph_prev_by_hiv(coinf_df, ax=ax, start_year=2000)

    # Panel B: New infections by sex (time series)
    ax = fig.add_subplot(gs[0, 1])
    plot_infections_by_sex(coinf_df, ax=ax, start_year=2000)

    # Panel C: Syphilis prevalence by age
    ax = fig.add_subplot(gs[0, 2])
    plot_prevalence_by_age(epi_df, disease='syphilis', ax=ax)

    # Panel D: Syphilis infections by sex work status
    ax = fig.add_subplot(gs[1, 0])
    plot_infections_by_sw(sw_df, disease='syphilis', ax=ax)

    # Panel E: HIV infections by sex work status
    ax = fig.add_subplot(gs[1, 1])
    plot_infections_by_sw(sw_df, disease='hiv', ax=ax)

    # Panel F: HIV prevalence by age
    ax = fig.add_subplot(gs[1, 2])
    plot_prevalence_by_age(epi_df, disease='hiv', ax=ax)

    # Add panel labels
    pl.figtext(0.02, 0.95, 'A', fontsize=40, ha='center', va='center', weight='bold')
    pl.figtext(0.35, 0.95, 'B', fontsize=40, ha='center', va='center', weight='bold')
    pl.figtext(0.68, 0.95, 'C', fontsize=40, ha='center', va='center', weight='bold')
    pl.figtext(0.02, 0.48, 'D', fontsize=40, ha='center', va='center', weight='bold')
    pl.figtext(0.35, 0.48, 'E', fontsize=40, ha='center', va='center', weight='bold')
    pl.figtext(0.68, 0.48, 'F', fontsize=40, ha='center', va='center', weight='bold')

    # Save figure
    figname = 'fig1_syph_hiv_epi'
    pl.savefig(f"figures/{figname}.png", dpi=100)

    print('Done.')
