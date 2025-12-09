# syph_dx_zim

**Estimating the value of novel syphilis diagnostics in Zimbabwe: how much overtreatment can be avoided?**

This repository contains an agent-based model of co-transmitting HIV and syphilis in Zimbabwe, used to evaluate the potential impact of novel point-of-care diagnostic tests on syphilis overtreatment.

## Background

Current approaches to syphilis diagnosis in Zimbabwe result in substantial overtreatment:

1. **Symptomatic detection (genital ulcer disease)**: WHO guidelines recommend syndromic management, which treats all patients presenting with genital ulcers for multiple pathogens including syphilis. This achieves 100% sensitivity but 0% specificity, as syphilis causes only 10-30% of GUD cases in sub-Saharan Africa.

2. **Asymptomatic detection (antenatal screening)**: Dual HIV-syphilis rapid diagnostic tests use treponemal antibody detection, which cannot distinguish active infection from past/treated infection. Studies suggest 40-60% of treponemal-positive pregnant women may have past rather than active infections.

This study uses STIsim, an agent-based model of co-transmitting sexually transmitted infections, to quantify current levels of overtreatment and estimate the reduction achievable through two novel point-of-care diagnostics:
- A rapid test for detecting treponemes directly from genital ulcers
- A confirmatory test for active syphilis following positive treponemal screening in pregnancy

## Installation

1. Create a new virtual environment:
   ```bash
   conda create -n syph_dx python=3.11 -y
   conda activate syph_dx
   ```

2. Install requirements:
   ```bash
   pip install -r requirements.txt
   ```

3. Test the installation:
   ```bash
   python run_sims.py
   ```

## Repository Structure

### Data
- `data/` - Input data for calibration
  - `zimbabwe_hiv_data.csv` - HIV prevalence, incidence, and mortality targets
  - `zimbabwe_syph_data.csv` - Syphilis prevalence and symptomatic presentation rates
  - Other demographic and behavioral data files

### Model Files
- `diseases.py` - Transmission models for HIV, syphilis, and background GUD
- `run_sims.py` - Construct main HIV-syphilis coinfection model and define run configs
- `interventions.py` - Diagnostic intervention implementations
- `analyzers.py` - Custom analyzers for tracking overtreatment and health outcomes
- `utils.py` - Utility functions for data processing and plotting

### Analysis Workflow

The analysis proceeds in four sequential steps:

#### Step 1: Calibrate HIV Model
```bash
python run_calibrations.py
```
Calibrates HIV prevalence, incidence, and mortality to match Zimbabwe data.

**Outputs:**
- `results/zim_calib_stats_hiv.df`
- `results/zim_par_stats_hiv.df`

#### Step 2: Calibrate Syphilis Model
```bash
python run_calibrations.py
```
Calibrates syphilis prevalence and symptomatic presentation rates.

**Outputs:**
- `results/zim_calib_stats_syph.df`
- `results/zim_par_stats_syph.df`

#### Step 3: Run Diagnostic Scenarios
```bash
python run_scenarios.py
```
Simulates baseline and intervention scenarios to quantify overtreatment reduction.

**Scenarios:**
- **Baseline**: Current syndromic management + treponemal RDT screening
- **Scenario 1**: Baseline + POC test for GUD (detects treponemes from ulcers)
- **Scenario 2**: Baseline + POC confirmatory test for ANC (confirms active infection)
- **Scenario 3**: Both POC tests

**Outputs:**
- `results/overtreatment.obj`
- `results/diagnostic_scenarios.obj`
- `results/health_outcomes.obj`

### Plotting Results

Generate figures for the manuscript:

```bash
# Calibration results (supplementary materials)
python plot_calibrations.py

# Epidemiological trends (Figure 2)
python plot_fig2_epi.py

# Overtreatment analysis (Figure 3)
python plot_fig3_overtreatment.py

# Diagnostic scenario comparison (Figure 4)
python plot_fig4_scenarios.py
```

## Key Model Features

### HIV Transmission
- Sexual network structure with casual, main, and one-time partnerships
- HIV transmission with stage-specific infectiousness
- ART coverage and viral suppression
- Interaction with syphilis (increased HIV acquisition risk)

### Syphilis Transmission
- Primary, secondary, early latent, and late latent stages
- Symptomatic presentation (genital ulcers) in ~40% of primary/secondary cases
- Congenital syphilis transmission during pregnancy
- Interaction with HIV (increased syphilis acquisition and progression)

### Diagnostic Pathways
1. **Symptomatic (GUD) pathway:**
   - Probability of care-seeking with genital ulcers
   - Syndromic management (baseline): treat all with antibiotics
   - POC treponeme test (intervention): test-guided treatment

2. **Asymptomatic (ANC) pathway:**
   - Antenatal care attendance rates
   - Dual HIV-syphilis screening
   - Treponemal RDT (baseline): treat all positives
   - POC confirmatory test (intervention): test-guided treatment

### Outcome Measures
- Syphilis overtreatment (unnecessary treatments)
- Syphilis undertreatment (missed cases)
- Congenital syphilis cases
- Disability-adjusted life years (DALYs)
- Treatment costs
- Benzathine penicillin doses used

## Model Parameters

Key parameters are calibrated to match Zimbabwe-specific data:
- HIV prevalence: ~12-13% in adults (ZIMPHIA 2015-2016)
- Active syphilis prevalence: ~0.8-0.9% in general population
- Syphilis seroprevalence: ~2.7% (ever infected)
- ANC attendance: ~90%
- Proportion of GUD caused by syphilis: ~15-25%

## Requirements

See `requirements.txt` for full dependencies. Key packages:
- `starsim` (v2.0+)
- `stisim` (v1.0+)
- `numpy`
- `pandas`
- `matplotlib`
- `scienceplots`
- `scipy`

## Citation

If you use this model or code, please cite:

[Citation will be added upon publication]

## License

MIT License - see LICENSE file for details

## Contact

For questions or collaboration inquiries, please open an issue on this repository.

## Acknowledgments

This work uses the Starsim modeling framework developed by the Institute for Disease Modeling. 
