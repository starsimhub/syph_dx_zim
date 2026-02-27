"""
Plot Figure 2: Treatment outcomes under standard of care.

Panels:
    A: Stacked treatment outcomes over time (all pathways)
    B: Stacked bars — correct + overtreated + missed by pathway
    C: Overtreatment rate (%) by pathway
    D: Detection rate over time (treated / active prevalence)
    E: Stage at detection by pathway
    F: Missed infections over time by pathway
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

PATHWAY_LABELS = {
    'gud_syndromic': 'GUD\nsyndromic',
    'anc_screen': 'ANC\nscreening',
    'secondary_rash': 'Secondary\nrash',
    'kp_screen': 'KP dual\nRDT',
    'plhiv_screen': 'PLHIV\nscreening',
    'newborn': 'Newborn',
}
PATHWAY_COLORS = {
    'gud_syndromic': '#e41a1c',
    'anc_screen': '#377eb8',
    'secondary_rash': '#ff7f00',
    'kp_screen': '#984ea3',
    'plhiv_screen': '#ff69b4',
    'newborn': '#4daf4a',
}

OC_COLORS = {
    'success': '#4daf4a',
    'unnecessary': '#e41a1c',
    'missed': '#984ea3',
}

PATHWAYS = ['gud_syndromic', 'anc_screen', 'kp_screen', 'plhiv_screen', 'newborn']


def load_data(scenario='soc'):
    return sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df')


def pivot_annual(df, start_year=2000, end_year=2035):
    df = df[(df.year >= start_year) & (df.year <= end_year)].copy()
    df['year_int'] = df.year.astype(int)
    return df.groupby(['scenario', 'year_int', 'metric']).value.agg(['mean', 'std', 'count']).reset_index()


def get_metric(df, metric, start_year=2000, end_year=2035):
    sub = df[(df.metric == metric) & (df.year_int >= start_year) & (df.year_int <= end_year)]
    return sub.set_index('year_int')['mean']


def plot_stacked_outcomes_ts(df, ax, start_year=2000, end_year=2035):
    """Panel A: Stacked area — correctly treated, overtreated, missed across ALL pathways over time"""
    # Sum across all pathways for each outcome
    success_total = None
    unnecessary_total = None
    missed_total = None
    for pw in PATHWAYS:
        s = get_metric(df, f'{pw}_success', start_year, end_year)
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        m = get_metric(df, f'{pw}_missed', start_year, end_year)
        if success_total is None:
            success_total, unnecessary_total, missed_total = s, u, m
        else:
            success_total = success_total.add(s, fill_value=0)
            unnecessary_total = unnecessary_total.add(u, fill_value=0)
            missed_total = missed_total.add(m, fill_value=0)

    years = success_total.index
    ax.stackplot(years, success_total, unnecessary_total, missed_total,
                 labels=['Correctly treated', 'Overtreated', 'Missed'],
                 colors=[OC_COLORS['success'], OC_COLORS['unnecessary'], OC_COLORS['missed']],
                 alpha=0.85)

    ax.set_ylabel('People per year')
    ax.set_title('Treatment outcomes over time')
    ax.legend(frameon=False, fontsize=13, loc='upper left')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_stacked_outcomes(df, ax, start_year=2000, end_year=2035):
    """Panel B: Stacked bars — correct (green) + overtreated (red) + missed (purple) by pathway"""
    x = np.arange(len(PATHWAYS))
    width = 0.6

    success_vals, unnecessary_vals, missed_vals = [], [], []
    for pw in PATHWAYS:
        success_vals.append(get_metric(df, f'{pw}_success', start_year, end_year).mean())
        unnecessary_vals.append(get_metric(df, f'{pw}_unnecessary', start_year, end_year).mean())
        missed_vals.append(get_metric(df, f'{pw}_missed', start_year, end_year).mean())

    s = np.array(success_vals)
    u = np.array(unnecessary_vals)
    m = np.array(missed_vals)

    ax.bar(x, s, width, label='Correctly treated', color=OC_COLORS['success'], alpha=0.85)
    ax.bar(x, u, width, bottom=s, label='Overtreated', color=OC_COLORS['unnecessary'], alpha=0.85)
    ax.bar(x, m, width, bottom=s + u, label='Missed', color=OC_COLORS['missed'], alpha=0.85)

    # Total labels
    for i in range(len(PATHWAYS)):
        total = s[i] + u[i] + m[i]
        if total > 0:
            ax.text(x[i], total + total * 0.02, f'{total:,.0f}',
                    ha='center', va='bottom', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in PATHWAYS])
    ax.set_ylabel('Annual average')
    ax.set_title('Treatment outcomes by pathway')
    ax.legend(frameon=False, fontsize=13)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_rate_bars(df, ax, start_year=2000, end_year=2035):
    """Panel C: Overtreatment rate (%) by pathway"""
    x = np.arange(len(PATHWAYS))

    rates = []
    for pw in PATHWAYS:
        s = get_metric(df, f'{pw}_success', start_year, end_year).mean()
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year).mean()
        total = s + u
        rates.append(u / total * 100 if total > 0 else 0)

    bars = ax.bar(x, rates, color=[PATHWAY_COLORS[pw] for pw in PATHWAYS], alpha=0.85)
    for bar, rate in zip(bars, rates):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{rate:.0f}%', ha='center', va='bottom', fontsize=16, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in PATHWAYS])
    ax.set_ylabel('Overtreatment rate (%)')
    ax.set_title('Overtreatment rate by pathway')
    ax.set_ylim(0, 110)


def plot_detection_rate(df, ax, start_year=2000, end_year=2035):
    """Panel D: Detection rate — fraction of cases treated each year (overall and active only)"""
    # Total treated (success + unnecessary) and correctly treated across all pathways
    success_total = None
    treated_total = None
    for pw in PATHWAYS:
        s = get_metric(df, f'{pw}_success', start_year, end_year)
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        if success_total is None:
            success_total = s
            treated_total = s.add(u, fill_value=0)
        else:
            success_total = success_total.add(s, fill_value=0)
            treated_total = treated_total.add(s, fill_value=0).add(u, fill_value=0)

    n_active = get_metric(df, 'n_active', start_year, end_year)
    n_infected = get_metric(df, 'n_infected', start_year, end_year)

    active_rate = (success_total / n_active * 100).clip(upper=100)
    overall_rate = (treated_total / n_infected * 100).clip(upper=100)

    ax.plot(overall_rate.index, overall_rate, color='#377eb8', linewidth=2.5, label='All infected')
    ax.plot(active_rate.index, active_rate, color='#e41a1c', linewidth=2.5, label='Active only')
    ax.set_ylabel('Cases treated (%)')
    ax.set_title('Detection rate')
    ax.legend(frameon=False, fontsize=13)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(0, 100)


STAGES = ['primary', 'secondary', 'early', 'late', 'tertiary']
STAGE_LABELS = ['Primary', 'Secondary', 'Early\nlatent', 'Late\nlatent', 'Tertiary']
STAGE_COLORS = ['#e41a1c', '#ff7f00', '#984ea3', '#377eb8', '#4daf4a']


def plot_stage_at_detection(df, ax, start_year=2000, end_year=2035):
    """Panel E: Stage at detection — what stage were correctly treated people in, by pathway"""
    det_pathways = [pw for pw in PATHWAYS if pw != 'newborn']
    x = np.arange(len(det_pathways))
    width = 0.6

    # Collect stage counts per pathway
    stage_vals = {stage: [] for stage in STAGES}
    for pw in det_pathways:
        for stage in STAGES:
            stage_vals[stage].append(get_metric(df, f'{pw}_stage_{stage}', start_year, end_year).mean())

    # Plot as stacked bars (proportions)
    totals = np.zeros(len(det_pathways))
    for stage in STAGES:
        totals += np.array(stage_vals[stage])

    bottom = np.zeros(len(det_pathways))
    for stage, label, color in zip(STAGES, STAGE_LABELS, STAGE_COLORS):
        vals = np.array(stage_vals[stage])
        pcts = np.where(totals > 0, vals / totals * 100, 0)
        ax.bar(x, pcts, width, bottom=bottom, label=label.replace('\n', ' '),
               color=color, alpha=0.85)
        bottom += pcts

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in det_pathways])
    ax.set_ylabel('Stage at detection (%)')
    ax.set_title('Stage at detection')
    ax.legend(frameon=False, fontsize=12, loc='upper right')
    ax.set_ylim(0, 105)


def plot_missed_ts(df, ax, start_year=2000, end_year=2035):
    """Panel F: Missed active infections over time, stacked by pathway"""
    missed_by_pw = {}
    for pw in PATHWAYS:
        missed_by_pw[pw] = get_metric(df, f'{pw}_missed', start_year, end_year)

    # Align all series to the same index
    years = missed_by_pw[PATHWAYS[0]].index
    stacks = [missed_by_pw[pw].reindex(years, fill_value=0) for pw in PATHWAYS]

    ax.stackplot(years, *stacks,
                 labels=[PATHWAY_LABELS[pw].replace('\n', ' ') for pw in PATHWAYS],
                 colors=[PATHWAY_COLORS[pw] for pw in PATHWAYS],
                 alpha=0.85)

    ax.set_ylabel('Missed cases per year')
    ax.set_title('Missed infections by pathway')
    ax.legend(frameon=False, fontsize=12, loc='upper left')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


if __name__ == '__main__':

    scenario = 'soc'
    df = load_data(scenario)
    annual = pivot_annual(df, start_year=2000, end_year=2035)

    set_font(size=18)
    fig = pl.figure(figsize=(22, 14))
    gs = GridSpec(2, 3, left=0.06, right=0.98, bottom=0.06, top=0.93,
                  wspace=0.28, hspace=0.38)

    ax = fig.add_subplot(gs[0, 0])
    plot_stacked_outcomes_ts(annual, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_stacked_outcomes(annual, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_overtreatment_rate_bars(annual, ax)

    ax = fig.add_subplot(gs[1, 0])
    plot_detection_rate(annual, ax)

    ax = fig.add_subplot(gs[1, 1])
    plot_stage_at_detection(annual, ax)

    ax = fig.add_subplot(gs[1, 2])
    plot_missed_ts(annual, ax)

    for i, ax in enumerate(fig.axes):
        ax.text(-0.10, 1.06, chr(65 + i), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top')

    pl.savefig(f'{FIGURES_DIR}/fig2_treatment_outcomes_soc.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig2_treatment_outcomes_soc.png')
