"""
Plot Figure 2: Treatment outcomes under standard of care.

Panels:
    A: Stacked treatment outcomes over time (all pathways)
    B: Stacked bars — correct + overtreated + missed by pathway
    C: Overtreatment rate (%) by pathway
    D: Care-seeking gap — symptomatic vs presenting vs treated
    E: Stage at detection by pathway
    F: Congenital syphilis outcomes over time
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.gridspec import GridSpec
from utils import set_font, get_metric

RESULTS_DIR = 'results'
FIGURES_DIR = 'figures'

END_YEAR = 2025
FIG3_START_YEAR = 2020  # Fig 3 averages over this period (close to Fig 4 baseline year)

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
    'plhiv_screen': '#ff7f00',
    'newborn': '#4daf4a',
}

OC_COLORS = {
    'success': '#4daf4a',
    'unnecessary': '#e41a1c',
    'missed': '#984ea3',
}

PATHWAYS = ['gud_syndromic', 'anc_screen', 'kp_screen', 'plhiv_screen']  # adult pathways only


def load_data(scenario='soc'):
    return sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df').copy()



def plot_stacked_outcomes_ts(df, ax, start_year=2000, end_year=END_YEAR):
    """Panel A: Stacked area — correctly treated, overtreated, missed across ALL pathways over time"""
    success_total = None
    unnecessary_total = None
    for pw in PATHWAYS:
        s = get_metric(df, f'{pw}_success', start_year, end_year)
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year)
        if success_total is None:
            success_total, unnecessary_total = s, u
        else:
            success_total = success_total.add(s, fill_value=0)
            unnecessary_total = unnecessary_total.add(u, fill_value=0)

    years = success_total.index
    ax.stackplot(years, success_total, unnecessary_total,
                 labels=['Correctly treated', 'Overtreated'],
                 colors=[OC_COLORS['success'], OC_COLORS['unnecessary']],
                 alpha=0.85)

    ax.set_ylabel('People per year')
    ax.set_title('Treatment outcomes over time')
    ax.legend(frameon=False, fontsize=13, loc='upper left')
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_stacked_outcomes(df, ax, start_year=FIG3_START_YEAR, end_year=END_YEAR):
    """Panel A: Stacked bars — correctly treated (solid) + overtreated (faded+hatched) by use case.
    Colors match Fig 4A: pathway color for both components."""
    from matplotlib.patches import Patch
    x = np.arange(len(PATHWAYS))
    width = 0.6

    for i, pw in enumerate(PATHWAYS):
        s = get_metric(df, f'{pw}_success',     start_year, end_year).mean()
        u = get_metric(df, f'{pw}_unnecessary', start_year, end_year).mean()
        color = PATHWAY_COLORS[pw]
        ax.bar(x[i], s, width, color=color, alpha=0.85)
        ax.bar(x[i], u, width, bottom=s, color=color, alpha=0.30,
               hatch='////', edgecolor=color, linewidth=0.5)
        total = s + u
        if total > 0:
            ax.text(x[i], total + total * 0.02, f'{total:,.0f}',
                    ha='center', va='bottom', fontsize=12)

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in PATHWAYS])
    ax.set_ylabel(f'Annual average ({start_year}–{end_year})')
    ax.set_title('(A) Treatment outcomes\nby use case')
    ax.legend(handles=[
        Patch(facecolor='#888888', alpha=0.85, label='Correctly treated'),
        Patch(facecolor='#888888', alpha=0.30, hatch='////', edgecolor='#888888',
              label='Overtreated'),
    ], frameon=False, fontsize=14)
    ax.set_ylim(bottom=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    sc.SIticks(ax, axis='y')


def plot_overtreatment_rate_bars(df, ax, start_year=FIG3_START_YEAR, end_year=END_YEAR):
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
    ax.set_title(f'(B) Overtreatment rate\nby pathway ({start_year}–{end_year} avg)')
    ax.set_ylim(0, 110)


def plot_care_seeking_gap(df, ax, start_year=2000, end_year=END_YEAR):
    """Panel D: Care-seeking gap — active cases vs treated (correct + overtreat) vs correctly treated"""
    n_active = get_metric(df, 'n_active', start_year, end_year)

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

    ax.fill_between(n_active.index, 0, n_active, alpha=0.15, color='#555555')
    ax.plot(n_active.index, n_active, color='#555555', linewidth=2.5, label='Active infections')
    ax.plot(treated_total.index, treated_total, color='#e41a1c', linewidth=2.5, label='Total treated')
    ax.plot(success_total.index, success_total, color='#4daf4a', linewidth=2.5, label='Correctly treated')

    ax.set_ylabel('People per year')
    ax.set_title('Care-seeking gap')
    ax.legend(frameon=False, fontsize=12)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


STAGES = ['primary', 'secondary', 'early', 'late', 'tertiary']
STAGE_LABELS = ['Primary', 'Secondary', 'Early latent', 'Late latent', 'Tertiary']
import matplotlib.cm as _cm
STAGE_COLORS = [_cm.magma(v) for v in [0.15, 0.35, 0.55, 0.72, 0.88]]


def plot_stage_at_detection(df, ax, start_year=FIG3_START_YEAR, end_year=END_YEAR):
    """Panel E: Stage at detection — what stage were correctly treated people in, by pathway"""
    det_pathways = [pw for pw in PATHWAYS if pw != 'newborn']
    x = np.arange(len(det_pathways))
    width = 0.6

    stage_vals = {stage: [] for stage in STAGES}
    for pw in det_pathways:
        for stage in STAGES:
            stage_vals[stage].append(get_metric(df, f'{pw}_stage_{stage}', start_year, end_year).mean())

    totals = np.zeros(len(det_pathways))
    for stage in STAGES:
        totals += np.array(stage_vals[stage])

    bottom = np.zeros(len(det_pathways))
    for stage, label, color in zip(STAGES, STAGE_LABELS, STAGE_COLORS):
        vals = np.array(stage_vals[stage])
        pcts = np.where(totals > 0, vals / totals * 100, 0)
        ax.bar(x, pcts, width, bottom=bottom, label=label,
               color=color, alpha=0.85)
        bottom += pcts

    ax.set_xticks(x)
    ax.set_xticklabels([PATHWAY_LABELS[pw] for pw in det_pathways])
    ax.set_ylabel('Stage at detection (%)')
    ax.set_title(f'(C) Stage at detection\nby pathway ({start_year}–{end_year} avg)')
    ax.legend(frameon=False, fontsize=14, loc='upper right', ncol=2)
    ax.set_ylim(0, 120)


def plot_congenital_outcomes_ts(df, ax, start_year=2000, end_year=END_YEAR):
    """Panel F: Congenital syphilis outcomes over time"""
    new_congenital = get_metric(df, 'new_congenital', start_year, end_year)

    # Newborn treatment outcomes
    nb_success = get_metric(df, 'newborn_success', start_year, end_year)
    nb_unnecessary = get_metric(df, 'newborn_unnecessary', start_year, end_year)

    if len(new_congenital) > 0:
        ax.fill_between(new_congenital.index, 0, new_congenital, alpha=0.2, color='#cb181d')
        ax.plot(new_congenital.index, new_congenital, color='#cb181d', linewidth=2.5, label='Congenital cases')

    if len(nb_success) > 0:
        nb_treated = nb_success.add(nb_unnecessary, fill_value=0)
        ax.plot(nb_treated.index, nb_treated, color='#377eb8', linewidth=2, linestyle='--', label='Newborns treated')
        ax.plot(nb_success.index, nb_success, color='#4daf4a', linewidth=2, linestyle='--', label='Correctly treated')

    ax.set_ylabel('Cases per year')
    ax.set_title('Congenital syphilis')
    ax.legend(frameon=False, fontsize=12)
    ax.set_xlim(start_year, end_year)
    ax.set_ylim(bottom=0)
    sc.SIticks(ax, axis='y')


def plot_gud_cascade(df, cs, ax, start_year=2015, end_year=END_YEAR):
    """GUD syndromic care-seeking cascade: per 100 primary syphilis infections"""

    # p_visible: weighted by model incidence sex ratio (correct denominator for "per 100 new primary cases")
    # p_symp_primary: 30% female, 80% male (updated in recalibration 2.2)
    new_inf_f = cs.loc[cs.index >= start_year, ('syph.new_infections_f', '50%')].mean()
    new_inf_m = cs.loc[cs.index >= start_year, ('syph.new_infections_m', '50%')].mean()
    new_inf = new_inf_f + new_inf_m
    p_visible = (new_inf_f * 0.3 + new_inf_m * 0.8) / new_inf if new_inf > 0 else 0.5

    gud_success = get_metric(df, 'gud_syndromic_success', start_year, end_year).mean()
    gud_missed = get_metric(df, 'gud_syndromic_missed', start_year, end_year).mean()
    gud_failure = get_metric(df, 'gud_syndromic_failure', start_year, end_year).mean()
    gud_unnecessary = get_metric(df, 'gud_syndromic_unnecessary', start_year, end_year).mean()

    # Cascade per 100 primary infections.
    # new_inf (from calib_stats "new_" prefix) uses SUM aggregation → annual total.
    # treatment_outcomes results (no "new_"/"n_" prefix) use MEAN aggregation → per-timestep avg.
    # Multiply by 12 to convert monthly averages to annual counts.
    # TODO: fix properly by setting summarize_by='sum' in treatment_outcomes Result definitions.
    visible_total = new_inf * p_visible
    sought_care_total = (gud_success + gud_missed + gud_failure) * 12
    seek_rate = sought_care_total / visible_total if visible_total > 0 else 0

    # Annotation values: non-syphilis GUD treated per 100 primary cases, PPV
    false_pos_per_100 = (gud_unnecessary * 12) / new_inf * 100 if new_inf > 0 else 0
    ppv_denom = (gud_success + gud_unnecessary + gud_failure) * 12
    ppv = (gud_success * 12) / ppv_denom * 100 if ppv_denom > 0 else 0

    # Steps 4+5 merged: syndromic management doesn't involve a discrete test —
    # the clinician assesses and treats presumptively. Combine sensitivity (80%) and
    # drug efficacy (98%) into a single "presumptively treated" step.
    steps = [
        100,
        100 * p_visible,
        100 * p_visible * seek_rate,
        100 * p_visible * seek_rate * 0.80 * 0.98,
    ]
    labels = [
        'Primary\nsyphilis',
        'Visible\nchancre',
        'Seek\ncare',
        'Presumptively\ntreated\nfor syphilis',
    ]
    loss_labels = [
        '',
        'Painless/internal',
        'No care-seeking',
        'Missed or treatment failure',
    ]

    # Colors: green for retained, with fading
    bar_color = '#4a90d9'
    loss_color = '#dddddd'
    y = np.arange(len(steps))[::-1]

    # Main bars
    bars = ax.barh(y, steps, color=bar_color, alpha=0.85, edgecolor='white', linewidth=0.5, height=0.7)

    # Loss shading (show what was lost at each step)
    for i in range(1, len(steps)):
        ax.barh(y[i], steps[i-1] - steps[i], left=steps[i], color=loss_color, alpha=0.5, height=0.7)

    # Labels
    for i, (step, label) in enumerate(zip(steps, labels)):
        if step >= 1:
            ax.text(step + 1, y[i], f'{step:.0f}', ha='left', va='center', fontsize=16, fontweight='bold')
        else:
            ax.text(max(step, 0) + 1, y[i], f'{step:.1f}', ha='left', va='center', fontsize=16, fontweight='bold')

    # Loss annotations
    for i in range(1, len(steps)):
        lost = steps[i-1] - steps[i]
        if lost > 0.5:
            ax.text(steps[i-1] - 0.5, y[i] + 0.35, f'\u2212{lost:.0f}: {loss_labels[i]}',
                    ha='right', va='bottom', fontsize=12, color='#888888', style='italic')

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 115)
    ax.set_title('GUD syndromic cascade\nper 100 primary syphilis cases')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_congenital_cascade(df, cs, ax, start_year=2015, end_year=END_YEAR, anc_attendance=0.90):
    """Congenital prevention cascade: per 100 at-risk pregnancies (mother infected with syphilis).

    Denominator uses pregnant syphilis prevalence from calib_stats.
    Treatment flows from treatment_outcomes are per-timestep averages (monthly dt),
    so they are annualized (×12) before use alongside annual calib_stats totals.
    """
    cols = cs.columns.get_level_values(0)

    # treatment_outcomes results (no "new_"/"n_" prefix) use MEAN aggregation → per-timestep avg.
    # calib_stats "new_mtc_*" uses SUM aggregation → annual total.
    # Multiply by 12 to convert monthly averages to annual counts before mixing.
    # TODO: fix properly by setting summarize_by='sum' in treatment_outcomes Result definitions.
    anc_success = get_metric(df, 'anc_screen_success', start_year, end_year).mean() * 12
    anc_missed = get_metric(df, 'anc_screen_missed', start_year, end_year).mean() * 12
    anc_failure = get_metric(df, 'anc_screen_failure', start_year, end_year).mean() * 12
    nb_success = get_metric(df, 'newborn_success', start_year, end_year).mean() * 12

    # Total MTC transmissions per year (annual totals from calib_stats)
    mtc_total = sum(
        cs.loc[cs.index >= start_year, (f'transmission_by_stage.new_mtc_{s}', '50%')].mean()
        for s in ['primary', 'secondary', 'early', 'late', 'tertiary']
        if f'transmission_by_stage.new_mtc_{s}' in cols
    )

    # Convert everything to "pregnancies" (not MTC events).
    # beta_m2c=0.075/timestep over ~9 months → ~50% per-pregnancy MTC probability
    p_mtc = 1 - (1 - 0.075)**9  # ≈ 0.50

    # Denominator: total infected pregnancies (in pregnancy units)
    # = (MTC actual + MTC prevented) / p_mtc
    fetal_prevented_mtc = anc_success * p_mtc * 0.90  # MTC events prevented by treatment
    total_at_risk = (mtc_total + fetal_prevented_mtc) / p_mtc  # Convert to pregnancies

    # Mothers screened: successfully treated + treatment failures + false negatives
    mothers_screened_active = anc_success + anc_missed + anc_failure

    # Normalize to per 100 at-risk pregnancies
    f = 100 / total_at_risk if total_at_risk > 0 else 1
    attend_anc = 100 * anc_attendance
    screened = mothers_screened_active * f
    detected = (anc_success + anc_failure) * f  # tested positive (infected who tested pos)
    treated = anc_success * f  # mother treated successfully
    fetus_cured = anc_success * 0.90 * f  # 90% avg fetal treatment efficacy
    remaining = 100 - fetus_cured
    nb_saved = nb_success / mtc_total * remaining if mtc_total > 0 else 0

    # Two-part cascade: ANC pathway + newborn pathway
    steps = [100, attend_anc, screened, detected, treated, fetus_cured]
    labels = [
        'At-risk\npregnancies',
        'Attend\nANC',
        'Screened for\nsyphilis',
        'Test\npositive',
        'Mother\ntreated',
        'Fetus cured\nin utero',
    ]
    loss_labels = ['', 'No ANC', 'Not screened', 'False negative', 'Not treated', 'Treatment failure']

    # Newborn section: known = mother was ANC-positive but baby still at risk
    known_at_risk = screened - fetus_cured  # ANC false neg + fetal treatment failure

    # Colors
    bar_color = '#377eb8'
    loss_color = '#dddddd'
    y = np.arange(len(steps) + 4)[::-1]  # Extra space for newborn section (3 bars + gap)

    # ANC bars
    ax.barh(y[0:len(steps)], steps, color=bar_color, alpha=0.85,
            edgecolor='white', linewidth=0.5, height=0.7)
    for i in range(1, len(steps)):
        ax.barh(y[i], steps[i-1] - steps[i], left=steps[i], color=loss_color, alpha=0.5, height=0.7)

    # Newborn section: 3 bars
    nb_steps = [remaining, known_at_risk, nb_saved]
    nb_labels = ['Still at risk\nat birth', 'Flagged for\nnewborn test', 'Newborn\ntreated']
    nb_y = y[len(steps)+1:len(steps)+4]
    # "Still at risk" bar
    ax.barh(nb_y[0], remaining, color='#e41a1c', alpha=0.7, edgecolor='white', linewidth=0.5, height=0.7)
    ax.barh(nb_y[0], 100 - remaining, left=remaining, color=bar_color, alpha=0.3, height=0.7)
    # "Known" and "treated" bars
    ax.barh(nb_y[1], known_at_risk, color='#e41a1c', alpha=0.7, edgecolor='white', linewidth=0.5, height=0.7)
    ax.barh(nb_y[1], remaining - known_at_risk, left=known_at_risk, color=loss_color, alpha=0.5, height=0.7)
    ax.barh(nb_y[2], nb_saved, color='#e41a1c', alpha=0.7, edgecolor='white', linewidth=0.5, height=0.7)

    # Labels on bars
    all_steps = list(steps) + [None] + list(nb_steps)
    all_labels = list(labels) + [''] + list(nb_labels)
    for i, (step, _) in enumerate(zip(all_steps, all_labels)):
        if step is None:
            continue
        if step >= 1:
            ax.text(step + 1, y[i], f'{step:.0f}', ha='left', va='center', fontsize=16, fontweight='bold')
        elif step > 0:
            ax.text(max(step, 0) + 1, y[i], f'{step:.1f}', ha='left', va='center', fontsize=16, fontweight='bold')

    # Loss annotations for ANC steps
    for i in range(1, len(steps)):
        lost = steps[i-1] - steps[i]
        if lost > 1:
            ax.text(steps[i-1] - 0.5, y[i] + 0.35, f'\u2212{lost:.0f}: {loss_labels[i]}',
                    ha='right', va='bottom', fontsize=12, color='#888888', style='italic')

    # Newborn loss annotations
    ax.text(remaining - 0.5, nb_y[1] + 0.35, f'\u2212{remaining - known_at_risk:.0f}: Mother not screened',
            ha='right', va='bottom', fontsize=12, color='#888888', style='italic')

    # Section divider
    sep_y = (y[len(steps)-1] + y[len(steps)]) / 2
    ax.axhline(y=sep_y, color='grey', linewidth=0.8, linestyle='--', alpha=0.4)
    ax.text(50, sep_y - 0.15, 'Remaining: newborn pathway', ha='center', va='top',
            fontsize=13, color='grey', style='italic')

    all_y_labels = list(labels) + [''] + list(nb_labels)
    ax.set_yticks(y)
    ax.set_yticklabels(all_y_labels)
    ax.set_xlim(0, 115)
    ax.set_title('Congenital prevention cascade\nper 100 at-risk pregnancies')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


if __name__ == '__main__':

    scenario = 'soc'
    df = load_data(scenario)
    cs = sc.loadobj(f'{RESULTS_DIR}/zimbabwe_calib_stats_all.df')

    set_font(size=20)
    fig = pl.figure(figsize=(22, 8))
    gs = GridSpec(1, 3, left=0.05, right=0.98, bottom=0.10, top=0.88,
                  wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    plot_stacked_outcomes(df, ax)

    ax = fig.add_subplot(gs[0, 1])
    plot_overtreatment_rate_bars(df, ax)

    ax = fig.add_subplot(gs[0, 2])
    plot_stage_at_detection(df, ax)

    pl.savefig(f'{FIGURES_DIR}/fig3_treatment_outcomes_soc.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig3_treatment_outcomes_soc.png')

    # --- Cascade figures ---
    fig2, ax2 = pl.subplots(1, 1, figsize=(12, 10))
    plot_gud_cascade(df, cs, ax2)
    pl.tight_layout()
    pl.savefig(f'{FIGURES_DIR}/fig2_cascades_soc.png', dpi=200, bbox_inches='tight')
    print(f'Saved {FIGURES_DIR}/fig2_cascades_soc.png')
