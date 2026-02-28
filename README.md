# Estimating the value of novel syphilis diagnostics in Zimbabwe

Agent-based model analysis of syphilis diagnostic scenarios using STIsim.

## Quick start

```bash
# 1. Calibrate (2,000 Optuna trials, ~hours on HPC)
python run_calibrations.py

# 2. Run multi-sim with top 200 calibrated parameter sets (~20 min)
python run_msim.py --n_pars 200

# 3. Run diagnostic scenarios (5 scenarios x 10 pars, ~15 min)
python run_scenarios.py

# 4. Generate all figures
python plot_fig1_epi.py           # Fig 1: Syphilis & HIV epidemiology
python plot_fig2_treatment.py     # Fig 2: Care-seeking cascades + Fig 3: Treatment outcomes
python plot_fig4_scenarios.py     # Fig 4: Scenario comparison
python plot_figs2_network.py      # Fig S2: Network structure (supplementary)
python plot_figs3_calibration.py  # Fig S3: HIV calibration (supplementary)
```

## Pipeline

| Step | Script | Output | Time |
|------|--------|--------|------|
| Calibrate | `run_calibrations.py` | `results/zimbabwe_pars_all.df` | Hours (HPC) |
| Multi-sim | `run_msim.py` | `results/zimbabwe_calib_stats_all.df`, `results/sw_prev_df.df` | ~20 min |
| Scenarios | `run_scenarios.py` | `results/treatment_outcomes_{scenario}.df` | ~15 min |

- `run_msim.py` captures all results via `sim.to_df()` (~744 columns, pruned to ~100).
  To add new analyzer results, just rerun `run_msim.py` — no recalibration needed.
- `run_scenarios.py` runs 5 scenarios: `soc`, `gud`, `conf`, `both`, `cs`.

## Figures

| Figure | Script | Description |
|--------|--------|-------------|
| Fig 1 | `plot_fig1_epi.py` | Syphilis & HIV epidemiology (5 panels) |
| Fig 2 | `plot_fig2_treatment.py` | Care-seeking cascades (GUD + congenital) |
| Fig 3 | `plot_fig2_treatment.py` | Treatment outcomes under SOC (3 panels) |
| Fig 4 | `plot_fig4_scenarios.py` | Scenario comparison (3 panels) |
| Fig S2 | `plot_figs2_network.py` | Network structure (supplementary) |
| Fig S3 | `plot_figs3_calibration.py` | HIV calibration to UNAIDS data (supplementary) |

## Key files

| File | Description |
|------|-------------|
| `diseases.py` | HIV + syphilis disease configuration |
| `interventions.py` | Diagnostic testing algorithms and treatment pathways |
| `analyzers.py` | Treatment outcomes, transmission by stage, epi time series |
| `run_sims.py` | Core sim-building functions |
| `data/syph_dx.csv` | Diagnostic test sensitivities by syphilis state |
| `data/syph_products.csv` | Algorithm routing by scenario and year |
| `syph_dx_zim.md` | Main manuscript |
| `sm_syph_dx_zim.md` | Supplementary materials |
| `archive/` | Old/exploratory plotting scripts (not in manuscript) |

## Dependencies

- Python 3.11+
- [STIsim](https://github.com/starsimhub/stisim) v1.5+
- [Starsim](https://github.com/starsimhub/starsim) v3.1+
- sciris, numpy, pandas, matplotlib
