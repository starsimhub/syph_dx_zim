# Estimating the value of novel syphilis diagnostics in Zimbabwe

Agent-based model analysis of syphilis diagnostic scenarios using STIsim.

## Quick start

```bash
# 1. Calibrate (2,000 Optuna trials, ~hours on HPC)
python run_calibrations.py

# 2. Run multi-sim with all surviving calibrated parameter sets (~20 min)
python run_msim.py

# 3. Run diagnostic scenarios (6 scenarios × all pars, ~15 min)
python run_scenarios.py

# 4. Generate all figures
python plot_fig1_epi.py           # Fig 1: Syphilis & HIV epidemiology
python plot_fig2_treatment.py     # Fig 2: Care-seeking cascades + Fig 3: Treatment outcomes
python plot_fig4_scenarios.py     # Fig 4: Scenario comparison + cost threshold analysis
python plot_figs2_network.py      # Fig S2: Network structure (supplementary)
python plot_figs3_calibration.py  # Fig S3: HIV calibration (supplementary)
```

## Pipeline

| Step | Script | Output | Time |
|------|--------|--------|------|
| Calibrate | `run_calibrations.py` | `results/zimbabwe_pars_all.df` | Hours (HPC) |
| Multi-sim | `run_msim.py` | `results/zimbabwe_calib_stats_all.df`, `results/sw_prev_df.df` | ~20 min |
| Scenarios | `run_scenarios.py` | `results/treatment_outcomes_{scenario}.df` | ~15 min |

Notes:
- `run_msim.py` replays each calibrated parameter set with its stored `rand_seed`, ensuring exact reproducibility across the pipeline.
- `run_scenarios.py` runs 6 scenarios: `soc`, `gud`, `anc`, `kp`, `plhiv`, `both`.
- To add new analyzer results, rerun `run_msim.py` — no recalibration needed.

## Scenarios

| Scenario | Description |
|----------|-------------|
| `soc` | Standard of care (syndromic management throughout) |
| `gud` | GUD POC NT active infection diagnostic |
| `anc` | ANC POC NT active infection diagnostic (confirmatory) |
| `kp` | KP dual RDT + POC NT active infection diagnostic (confirmatory) |
| `plhiv` | PLHIV dual RDT + POC NT active infection diagnostic (confirmatory) |
| `both` | All four diagnostic use cases active simultaneously |

## Figures

| Figure | Script | Description |
|--------|--------|-------------|
| Fig 1 | `plot_fig1_epi.py` | Syphilis & HIV epidemiology (5 panels) |
| Fig 2 | `plot_fig2_treatment.py` | Care-seeking cascades (GUD + congenital) |
| Fig 3 | `plot_fig2_treatment.py` | Treatment outcomes under SOC (3 panels) |
| Fig 4 | `plot_fig4_scenarios.py` | Scenario comparison (panels A–C) + cost threshold analysis (panel D) |
| Fig S2 | `plot_figs2_network.py` | Network structure (supplementary) |
| Fig S3 | `plot_figs3_calibration.py` | HIV calibration to UNAIDS data (supplementary) |

## Key files

| File | Description |
|------|-------------|
| `diseases.py` | HIV + syphilis disease configuration (natural history parameters) |
| `interventions.py` | Diagnostic testing algorithms and treatment pathways |
| `analyzers.py` | Treatment outcomes, transmission by stage, epi time series |
| `run_sims.py` | Core sim-building functions (`make_sim`, `load_calib_pars`) |
| `utils.py` | Shared plotting utilities (`set_font`, `get_metric`) |
| `data/syph_dx.csv` | Diagnostic test sensitivities by syphilis state |
| `data/syph_products.csv` | Algorithm routing by scenario and year |
| `syph_dx_zim_v2.md` | Main manuscript (working file) |
| `sm_syph_dx_zim.md` | Supplementary materials |
| `archive/` | Old and exploratory scripts (not part of the manuscript pipeline) |

## Dependencies

- Python 3.11+
- [STIsim](https://github.com/starsimhub/stisim) v1.5.2
- [Starsim](https://github.com/starsimhub/starsim) v3.3.0
- sciris, numpy, pandas, matplotlib
