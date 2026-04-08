"""
Slide panel 3: "The Value" — unnecessary treatments avoided per test, by use case.

Single figure (3.83 × 3.83 in). Bar chart showing that each test avoids 0.4–0.7
unnecessary treatments, implying cost-savings when the test costs less than treatment.
"""

import sciris as sc
import numpy as np
import matplotlib.pyplot as pl
import sys; sys.path.insert(0, '..')
from utils import set_font

RESULTS_DIR = '../results'
FIGURES_DIR = '.'

PERIOD = (2026, 2035)

USE_CASES = [
    ('gud_syndromic', 'GUD\nPOC test',   '#E8453C'),
    ('plhiv_screen',  'PLHIV\nPOC NT',   '#ff7f00'),
    ('kp_screen',     'KP\nPOC NT',      '#984ea3'),
    ('anc_screen',    'ANC\nPOC NT',     '#377eb8'),
]


def load_scenario(scenario):
    df = sc.loadobj(f'{RESULTS_DIR}/treatment_outcomes_{scenario}.df').copy()
    sub = df[(df['year'] >= PERIOD[0]) & (df['year'] <= PERIOD[1])]
    return sub


def get_avoided_per_test():
    soc_sub  = load_scenario('soc')
    both_sub = load_scenario('both')
    results = []
    for pw, label, color in USE_CASES:
        soc_u = soc_sub[f'{pw}_unnecessary'].mean()
        both_u = both_sub[f'{pw}_unnecessary'].mean()
        soc_t = soc_sub[f'{pw}_treated'].mean()
        apt = (soc_u - both_u) / soc_t if soc_t > 0 else 0
        results.append((pw, label, color, apt))
    return results


def plot():
    set_font(size=12)
    fig, ax = pl.subplots(figsize=(3.83, 3.83), facecolor='white')
    fig.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.18)

    data = get_avoided_per_test()

    x = np.arange(len(data))
    width = 0.6
    labels = [d[1] for d in data]
    colors = [d[2] for d in data]
    values = [d[3] for d in data]

    bars = ax.bar(x, values, width, color=colors, alpha=0.85, edgecolor='white', linewidth=0.5)

    # Value labels above bars
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.02,
                f'{val:.2f}', ha='center', va='bottom', fontsize=11,
                fontweight='bold', color='#333333')

    # Reference line at 1.0
    ax.axhline(1.0, color='#999999', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.text(len(data) - 0.5, 1.02, '1 treatment', ha='right', va='bottom',
            fontsize=8, color='#999999')

    # Annotation
    ax.text(0.5, 0.95, transform=ax.transAxes,
            s='Cost-saving if test\ncosts less than treatment',
            fontsize=9, color='#2e7d32', fontweight='bold',
            ha='center', va='top', linespacing=1.3,
            bbox=dict(boxstyle='round,pad=0.4', facecolor='#e8f5e9', edgecolor='#c8e6c9',
                      linewidth=0.8, alpha=0.9))

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel('Unnecessary treatments\navoided per test', fontsize=10)
    ax.set_ylim(0, 1.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='x', length=0)

    pl.savefig(f'{FIGURES_DIR}/slide_cost.png', dpi=300, facecolor='white', bbox_inches=None)
    print(f'Saved {FIGURES_DIR}/slide_cost.png')
    pl.close(fig)


if __name__ == '__main__':
    plot()
