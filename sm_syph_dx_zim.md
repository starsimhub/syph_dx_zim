# Estimating the value of novel syphilis diagnostics in Zimbabwe: how much overtreatment can be avoided?

## Supplementary Materials

[**Syphilis model	2**](#syphilis-model)

[Natural history	2](#natural-history)

[Untreated infection	2](#untreated-infection)

[Treated infection	3](#treated-infection)

[Transmission	3](#transmission)

[Sexual transmission	3](#sexual-transmission)

[Maternal transmission	4](#maternal-transmission)

[Model Initialization	5](#model-initialization)

[Diagnostic Testing	6](#diagnostic-testing)

[**References	6**](#references)

# 

# Syphilis model  {#syphilis-model}

## Natural history {#natural-history}

### Untreated infection {#untreated-infection}

A depiction of the natural history of syphilis captured in our model is given in Figure S1.

![][image1]  
Figure S1: schematic of the modeled phases of adult syphilis in the absence of treatment.

After an agent is exposed to & acquires syphilis, there is an initial latent period of 10-90 days (mean 21 days) after which the primary phase of infection begins [(1,2)](https://www.zotero.org/google-docs/?JT7JOJ). This primary phase is characterized by the appearance of smooth chancres, which are often painless and/or occluded from sight, and therefore reportedly go unnoticed in around 60% of cases [(2)](https://www.zotero.org/google-docs/?Waa7Tb). Studies indicate that these chancres self-resolve within 3-10 weeks [(2,3)](https://www.zotero.org/google-docs/?w8eQGJ). The secondary phase of infection commences 2-26 weeks after initial infection, and is accompanied by a disseminated rash ranging in severity [(2,3)](https://www.zotero.org/google-docs/?Q0cU29). Secondary symptoms may persist for 3-12 weeks [(3,4)](https://www.zotero.org/google-docs/?gsTDiJ), after which the infection progresses to latency. Within our syphilis model, we determine all these durations (initial latency, primary infection, and secondary infection) for each agent at the point of infection. The durations are drawn from distributions shown in Table S1, which are parameterized to match the studies referenced above.

Especially within the early latent phase (variously defined as within the first 1-2 years), there is evidence that around 25-35% of infections can reactivate, returning to secondary syphilis [(2)](https://www.zotero.org/google-docs/?yos7Me). Within the model, when agents enter the latent phase we evaluate which of them will reactivate (using a probability of 0.35) and when, following the distribution given in Table S1. Reactivated infections will return to latency if untreated.

Latent infections can persist indefinitely or can progress to tertiary syphilis, a phase with diverse sequelae including neurologic, cardiovascular, and gummatous manifestations [(2,4)](https://www.zotero.org/google-docs/?EUO8Vl). Evidence indicates that tertiary syphilis commences 1-46 years after a person is initially infected with syphilis. In our model we assume that 35% of all untreated cases will progress to tertiary syphilis after dur\_latent years, where dur\_latent is assumed to be log-normally distributed with a mean of 20 years and standard deviation of 2 years (Table S1). Among those progressing to tertiary syphilis, we assume 5% will die from tertiary complications after an additional period drawn from a log-normal distribution with a mean of 5 years (Table S1).

| Parameter | Description | Value | Source |
| :---- | :---- | :---- | :---- |
| dur\_initial | Duration of initial latent period (weeks) | U(2, 10\) | [(2)](https://www.zotero.org/google-docs/?AP4Pp6) |
| dur\_primary | Duration of primary syphilis (weeks) | U(3, 10\) | [(2)](https://www.zotero.org/google-docs/?2VRFsw) |
| dur\_secondary | Duration of secondary syphilis (weeks) | LogN(mean=15, stdev=1.5) | [(2)](https://www.zotero.org/google-docs/?DoxisM) |
| p\_reactivation | Proportion of early latent cases that reactivate to secondary | 0.35 | [(2)](https://www.zotero.org/google-docs/?iYoUqO) |
| time\_to\_reactivate | Duration of latent infection before reactivation (years) | LogN(mean=1, stdev=0.25) | [(2)](https://www.zotero.org/google-docs/?UMPBG7) |
| p\_tertiary | Probability of a latent infection progressing to tertiary syphilis | 0.35 | [(2)](https://www.zotero.org/google-docs/?hM0InE) |
| time\_to\_tertiary | Duration of latent infection before tertiary syphilis (years) | LogN(mean=20, stdev=8) | [(2)](https://www.zotero.org/google-docs/?LgpuNc) |

Table S1: durations and probabilities of progressing through stages of adult syphilis

### Treated infection {#treated-infection}

The recommended treatment for syphilis is benzathine penicillin G (BPG), administered as a single intramuscular injection in cases of primary, secondary, or early latent syphilis, and a course of three injections for late latent syphilis or for pregnant women [(5)](https://www.zotero.org/google-docs/?wapNc8). A recent systematic review found that the recommended dosage of BPG is 90-100% successful in treating primary, secondary, and early latent infection [(6)](https://www.zotero.org/google-docs/?dxrRF2). The evidence on success rates for individuals with late latent or tertiary syphilis is more limited; one study reported lower success rates of around 60% [(6)](https://www.zotero.org/google-docs/?YsDudI), but this was averaged across a combination of different treatment regimens and did not include the current WHO-recommended three-dose schedule. Within the model, we assume a treatment success rate of 95%.

Following successful treatment, agents return to a susceptible state. In our model, we track treatment history through a persistent ever\_exposed state, which allows for modeling of serological test performance (see Diagnostic Testing section). We assume that successfully treated individuals have temporary immunity to reinfection, with immunity waning over time (parameters in Table S2). 

| Parameter | Description | Value | Source |
| ----- | ----- | ----- | ----- |
| treatment\_efficacy | Probability of successful treatment with BPG | 0.95 | (5,6) |
| dur\_immunity | Duration of post-treatment immunity | \[To be calibrated\] | Assumption |

Table S2: Treatment parameters

## Transmission  {#transmission}

### Sexual transmission {#sexual-transmission}

Syphilis can be transmitted when an infected person has sexual contact with a susceptible person. The probability of transmission varies according to the phase of infection, and is an active area of current research. During primary infection, *T pallidum* bacteria are highly concentrated at the site of the chancre and can easily be transmitted via contact, while during secondary infection the bacterial load is disseminated and systemic [(2)](https://www.zotero.org/google-docs/?X5HwJT). Various studies have estimated the probability of transmission during primary/secondary infection to be between 51-64% per partnership [(7,8)](https://www.zotero.org/google-docs/?Qd8Qoz). An accurate representation of the true likelihood of transmission would also need to account for behavioral changes during the early phases of infection.

Our model includes baseline parameters for the probability of infection per sexual act from males to females (beta\_m2f) and from females to males (beta\_f2m), which are estimated as part of the model calibration process. These baseline parameters are then scaled according to the phase of infection. After accounting for the effects of behavioral changes, bacterial load and concentration, we assume that the relative transmissibility during primary and secondary infection is the same, both set to a baseline value of 1.0 (Table S3).  
There is significant uncertainty about the extent of transmission during latent infection. Recent studies have identified *T pallidum* bacteria at mucosal sites of asymptomatically-infected people [(9)](https://www.zotero.org/google-docs/?ArDWWg), and ongoing work is investigating how readily these bacteria are able to transmit and establish infection during latent phases of infection. In our model, we scale the probability of transmission during latent infection using a factor rel\_trans\_latent, which we assume to be an exponentially-decaying function of how long a person has been in the latent phase. The transmissibility during latent infection is calculated as:  
rel\_trans(t) \= rel\_trans\_latent × exp(-ln(2)/half\_life × t)  
where t is the duration in latent phase and the half-life is assumed to be 1 year (Table S3). This formulation reflects the biological understanding that early latent infection is more infectious than late latent infection, with transmissibility declining over time. We assume that tertiary syphilis is not infectious (rel\_trans\_tertiary \= 0).

| Parameter | Description | Value | Source |
| :---- | :---- | :---- | :---- |
| beta\_m2f | Per-act transmission probability from an infected male to a susceptible female | 0.17 (0.15–0.22) | Calibrated |
| beta\_m2c | Per-timestep MTCT probability (stage-independent) | 0.075 | Calibrated to ~1,500 CS cases/yr |
| eff\_condom | Condom efficacy for syphilis prevention | 0.53 (0.41–0.58) | Calibrated |
| rel\_trans\_primary | Scale factor on transmission probability for persons with primary syphilis | 8.0 (6.9–9.7) | Calibrated |
| rel\_trans\_secondary | Scale factor on transmission probability for persons with secondary syphilis | 1.0 (reference) | Assumption |
| rel\_trans\_latent | Baseline scale factor for persons with latent syphilis | 0.1, decaying with half-life of 6 months | Assumption |
| rel\_trans\_tertiary | Scale factor for persons with tertiary syphilis | 0.0 | Assumption |

Table S3: probability of syphilis transmission through sexual contact. Values in parentheses show the 5th–95th percentile range across the top 200 calibrated parameter sets. Sexual transmission uses stage-specific relative infectiousness; maternal-to-child transmission uses stage-independent probability (rel_trans = 1 for all stages on the maternal network), reflecting that spirochetemia during pregnancy affects the fetus regardless of clinical stage.

### Maternal transmission  {#maternal-transmission}

Maternal transmission of syphilis usually occurs transplacentally, although transmission during birth is possible [(7)](https://www.zotero.org/google-docs/?Xp3Tep). Outcomes for the fetus depend on when transmission occurs, as well as the mother’s phases of infection during the course of the pregnancy. Studies of the prognoses are limited, but generally indicate significantly poorer outcomes in cases where the mother has primary, secondary, or early latent infection [(10)](https://www.zotero.org/google-docs/?lydgEi). Early latent infection is usually defined to mean the 1-2 years of latent infection, a period which is also associated with higher likelihood of secondary reactivation. In our model, we distinguish between three maternal infection states that determine birth outcomes:

1. **Active infection** (primary or secondary)  
2. **Early latent infection** (within the first 12-14 months of latent infection)  
3. **Late latent or tertiary infection**

We model five mutually exclusive birth outcomes:

* **Miscarriage** (pregnancy loss before viability)  
* **Neonatal death** (death within 28 days of birth)  
* **Stillbirth** (fetal death after 20 weeks gestation)  
* **Congenital syphilis** (live birth with congenital infection)  
* **Live birth without syphilis-related complications** (may include preterm birth or low birth weight but survives without active infection)

The probabilities of each outcome vary by maternal infection stage, as shown in Table S4. These probabilities are applied at the time of birth, with outcomes scheduled according to the timing of delivery in the model's pregnancy module

| Outcome | Active (Primary/Secondary) | Early Latent | Late Latent/Tertiary |
| ----- | ----- | ----- | ----- |
| Miscarriage | 0.00 | 0.00 | 0.00 |
| Neonatal death | 0.10 | 0.05 | 0.00 |
| Stillbirth | 0.20 | 0.10 | 0.10 |
| Congenital syphilis | 0.45 | 0.40 | 0.10 |
| Live birth without complications | 0.25 | 0.45 | 0.80 |
| **Total adverse outcomes** | **75%** | **55%** | **20%** |

Table S4: Probabilities of adverse birth and pregnancy outcomes associated with maternal syphilis

The model tracks all congenital outcomes separately in results, including cumulative burden over time. Neonatal deaths and stillbirths are processed through the mortality module, while congenital syphilis cases surviving infancy remain infected and contribute to the overall disease burden.

## Model Initialization {#model-initialization}

To initialize the model with realistic starting conditions reflecting historical syphilis prevalence, we use age- and sex-stratified prevalence data where available. The model supports two types of initial conditions:

1. **Active infection prevalence** (init\_prev\_data): Agents with primary or secondary infection at model start  
2. **Latent infection prevalence** (init\_prev\_latent\_data): Agents with latent infection at model start

For agents initialized with active infection, disease progression follows the standard natural history from their assigned stage. For agents initialized with latent infection, the model assigns:

* Time since infection (ti\_init\_infected): Determines position in latent phase  
* Probability of reactivation and tertiary progression: Same as for newly infected agents  
* Time to potential tertiary progression: Scheduled based on parameters in Table S1

This initialization approach ensures the model begins with a realistic distribution of disease stages and accounts for individuals at various points in their infection trajectory.

## Diagnostic Testing {#diagnostic-testing}

See Tables 1 and 2 in the main text for diagnostic test characteristics and intervention scenario descriptions. Test sensitivities by syphilis state are detailed in Table S5.

| Test | Naive | Previously exposed | Exposed (incubating) | Primary | Secondary | Latent | Tertiary |
|---|---|---|---|---|---|---|---|
| Syndromic GUD (SOC) | 80% | 80% | 80% | 80% | 80% | 80% | 80% |
| Dual treponemal RDT | 1% | 95% | 20% | 20% | 95% | 95% | 95% |
| GUD POC test | 5% | 5% | 5% | 95% | 20% | 5% | 5% |
| Confirmatory test | 0% | 0% | 10% | 95% | 95% | 95% | 95% |

Table S5: Diagnostic test sensitivity (probability of a positive result) by syphilis disease state. The syndromic GUD test represents presumptive treatment of all GUD presentations; 80% reflects that not all presentations result in treatment. The dual treponemal RDT has high sensitivity for previously exposed individuals (including those with resolved infection), which drives overtreatment. The GUD POC test is specific to primary chancres. The confirmatory test distinguishes active infection from serological scarring.

## Calibrated parameters {#calibrated-parameters}

We jointly calibrated 15 model parameters to Zimbabwe-specific HIV and syphilis epidemiological data using 2,000 Optuna trials. The top 200 parameter sets (ranked by weighted mismatch score) were retained for analysis. Table S6 lists all calibrated parameters with their prior ranges and posterior estimates.

| Parameter | Description | Prior range | Posterior median (5th–95th) |
|---|---|---|---|
| **HIV transmission** | | | |
| hiv\_beta\_m2f | Male-to-female HIV transmission probability per act | 0.002–0.014 | 0.012 (0.008–0.013) |
| hiv\_eff\_condom | Condom efficacy for HIV prevention | 0.50–0.90 | 0.86 (0.63–0.90) |
| hiv\_rel\_init\_prev | Relative scaling of initial HIV prevalence in 1985 | 2.0–15.0 | 6.7 (3.1–14.6) |
| hiv\_rel\_dur\_on\_art | Effective ART duration multiplier (base 3 years) | 1.0–20.0 | 13.0 (7.0–16.2) |
| **Syphilis transmission** | | | |
| syph\_beta\_m2f | Male-to-female syphilis transmission probability per act | 0.15–0.35 | 0.17 (0.15–0.22) |
| syph\_eff\_condom | Condom efficacy for syphilis prevention | 0.30–0.70 | 0.53 (0.41–0.58) |
| syph\_rel\_trans\_primary | Relative infectiousness of primary vs secondary syphilis | 3.0–10.0 | 8.0 (6.9–9.7) |
| **Sexual network** | | | |
| nw\_prop\_f0 | Proportion of women in lowest risk group | 0.55–0.90 | 0.67 (0.60–0.82) |
| nw\_prop\_m0 | Proportion of men in lowest risk group | 0.55–0.80 | 0.61 (0.55–0.77) |
| nw\_m1\_conc | Concurrency rate for mid-risk men | 0.05–0.30 | 0.08 (0.05–0.24) |
| **HIV-syphilis coinfection** | | | |
| conn\_rel\_sus\_syph\_hiv | Relative susceptibility to syphilis given HIV | 1.0–3.0 | 1.9 (1.2–2.6) |
| conn\_rel\_sus\_hiv\_syph | Relative susceptibility to HIV given active syphilis | 1.5–4.0 | 2.5 (1.7–3.9) |
| **Diagnostic testing rates** | | | |
| rel\_symp\_test | Symptomatic care-seeking rate multiplier | 0.50–1.50 | 0.84 (0.56–1.17) |
| rel\_anc\_test | ANC syphilis screening rate multiplier | 0.70–1.80 | 1.42 (1.08–1.66) |
| rel\_kp\_test | KP syphilis screening rate multiplier | 0.50–1.70 | 0.88 (0.53–1.61) |

Table S6: Calibrated model parameters. Prior ranges define the bounds of a uniform distribution explored by the Optuna sampler. Posterior estimates are computed from the top 200 parameter sets (out of 2,000 trials), ranked by weighted mean squared error mismatch to calibration targets.

# References {#references}

[1\.](https://www.zotero.org/google-docs/?mcvItb)	[Peeling RW, Hook III EW. The pathogenesis of syphilis: the Great Mimicker, revisited. The Journal of Pathology. 2006;208(2):224–32.](https://www.zotero.org/google-docs/?mcvItb) 

[2\.](https://www.zotero.org/google-docs/?mcvItb)	[O’Byrne P, MacPherson P. Syphilis. BMJ. 2019 Jun 28;365:l4159.](https://www.zotero.org/google-docs/?mcvItb) 

[3\.](https://www.zotero.org/google-docs/?mcvItb)	[Clark EG, Danbolt N. The Oslo study of the natural history of untreated syphilis: An epidemiologic investigation based on a restudy of the Boeck-Bruusgaard material a review and appraisal. Journal of Chronic Diseases. 1955 Sep 1;2(3):311–44.](https://www.zotero.org/google-docs/?mcvItb) 

[4\.](https://www.zotero.org/google-docs/?mcvItb)	[Singh AE, Romanowski B. Syphilis: Review with Emphasis on Clinical, Epidemiologic,  and Some Biologic Features. Clin Microbiol Rev. 1999 Apr;12(2):187–209.](https://www.zotero.org/google-docs/?mcvItb) 

[5\.](https://www.zotero.org/google-docs/?mcvItb)	[RECOMMENDATIONS FOR TREATMENT OF SYPHILIS. In: WHO Guidelines for the Treatment of Treponema pallidum (Syphilis) \[Internet\]. World Health Organization; 2016 \[cited 2024 Sep 17\]. Available from: https://www.ncbi.nlm.nih.gov/books/NBK384905/](https://www.zotero.org/google-docs/?mcvItb) 

[6\.](https://www.zotero.org/google-docs/?mcvItb)	[Clement ME, Okeke NL, Hicks CB. Treatment of Syphilis A Systematic Review. JAMA. 2014 Nov 12;312(18):1905–17.](https://www.zotero.org/google-docs/?mcvItb) 

[7\.](https://www.zotero.org/google-docs/?mcvItb)	[Stoltey JE, Cohen SE. Syphilis transmission: a review of the current evidence. Sex Health. 2015 Apr;12(2):103–9.](https://www.zotero.org/google-docs/?mcvItb) 

[8\.](https://www.zotero.org/google-docs/?mcvItb)	[Garnett GP, Aral SO, Hoyle DV, Cates W, Anderson RM. The natural history of syphilis. Implications for the transmission dynamics and control of infection. Sex Transm Dis. 1997 Apr;24(4):185–200.](https://www.zotero.org/google-docs/?mcvItb) 

[9\.](https://www.zotero.org/google-docs/?mcvItb)	[Aung ET, Fairley CK, Williamson DA, Azzato F, Towns JM, Wigan R, et al. Treponema pallidum Detection at Asymptomatic Oral, Anal, and Vaginal Sites in Adults Reporting Sexual Contact with Persons with Syphilis \- Volume 29, Number 10—October 2023 \- Emerging Infectious Diseases journal \- CDC. \[cited 2024 Sep 17\]; Available from: https://wwwnc.cdc.gov/eid/article/29/10/23-0660\_article](https://www.zotero.org/google-docs/?mcvItb) 

[10\.](https://www.zotero.org/google-docs/?mcvItb)	[Fiumara NJ, Fleming WL, Downing JG, Good FL. The Incidence of Prenatal Syphilis at the Boston City Hospital. New England Journal of Medicine. 1952 Jul 10;247(2):48–52.](https://www.zotero.org/google-docs/?mcvItb) 

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAloAAAAgCAYAAADKfdZkAAAUn0lEQVR4Xu2ceZQV1Z3Hm33f901QQQUFHcQoYgREZEBAIktQgrRKN7gighI30mqCbNp0u4AIyCIYXJFRiAguSCSaKKKZySTnzGQyOTNJ5p+Zc+acmT9m5tx531t16977+9V7za2uFvD8vud8zrv1u1W/qvr2fa9+71a9LisTiUQikUgkEjWolFCSv5TVTzSfwKmPaC6BUx9h/NN8gk99PiPE37oRfxuW+vgL0XyCj5b67b99LRTBNSqj1PYjnwpFgD/UsECpv911SCgC/KGGBYrlFCzwJyarWE7BBx5R0wLE8gk+8IiaFiiWU4iAN4lJtLgQLK5RGcWKC8ECf6hhgWKDW7DAH2pYoFhOwQJ/YrKK5RR84BE1LUAsn+ADj6hpgWI5hQh4k5hEiwvB4hqVUay4ECzwhxoWKDa4BQv8oYYFiuUULPAnJqtYTsEHHlHTAsTyCT7wiJoWKJZTiIA3iUm0uBAsrlEZxYoLwQJ/qGGBYoNbsMAfaligWE7BAn9isorlFHzgETUtQCyf4AOPqGmBYjmFCHiTmESLC8HiGpVRrLgQLPCHGhYoNrgFC/yhhgWK5RQs8Ccmq1hOwQceUdMCxPIJPvCImhYollOIgDeJSbS4ECyuURnFigvBAn+oYYFig1uwwB9qWKBYTsECf2KyiuUUfOARNS1ALJ/gA4+oaYFiOYUIeJOYRIsLweIalVGsuBAs8IcaFig2uAUL/KGGBYrlFCzwJyarWE7BBx5R0wLE8gk+8IiaFiiWU4iAN4lJtLgQLK5RGcWKC8ECf6hhgWKDW7DAH2pYoFhOwQJ/YrKK5RR84BE1LUAsn+ADj6hpgWI5hQh4k5hEiwvB4hqVUay4ECzwhxoWKDa4BQv8oYYFiuUULPAnJqtYTsEHHlHTAsTyCT7wiJoWKJZTiIA3iUm0uBAsrlEZxYoLwQJ/qGGBYoNbsMAfaligWE7BAn9isorlFHzgETUtQCyf4AOPqGmBYjmFCHiTmESLi5PNwc/2sdjJwjUqo1hxcTqxZvfr+vWsIeerTl27sX5QLH4iwB9qWKDY4P62cEaP3iwWCvyhhgWK5Wxofla9g8VOVeBPTFaxnKc6Z/fpr3698yCLNxTwiJoWIJZP8IFH1LRAsZynAvtPgc8ReJOYRIuLvIh3opY+stiLffGPn7J16XY0drJwjcooVlxQsA5o2aq1fqX93yQXjhyVtNfvP5gcT7uOHYseW7H4iRCfe92aWdWchmKxwR3CKz9en/gPGjdqrD7fuo+tdzKo77mZHJ5bxdRA/hqozyMvuJitY8hrn98EzjmV1tw13WkoFstZCqzfqkVL/Vox9QbW/02AfX/90nss3lBgf75lscqrOybtW9dd6fS4YvnyICTv5FHj1M83vsniIVw1YhSL5QXOhZrGVNxfiOUMBTko1333GrZeKa686NKkffSFt3SOmsVVbL1vkvhcIpNocZEXyN2pSyf96sZoobX/k71s2xPh098eYbG8cY0qqcpaVVZZ8wQNF8SKC5d5S+5XPfv1Y/G6eP7A+2rjgQ9Y/MlX3mCxtbt5TOd495DaevioFyt2vG6hRfPRbTYdOqy2fvQJy5EGtvXtKqLI3z+UVa5bQHrY4A7BFAA0Xhfv1e5iMROn3/ZxUTr09E/Zuh88u5vFvtz+brJ9luOiIAfxK13GX15wsZxZoD6n5S31DdT1u9hF6+21W1kMHNnwhvrqpQMsngc4j5jSqqg9rD2+tXoU6WE5i/HIzYtU/559Wdzwxdb9qef57rqX1LFt+70YxuQvt7zN1gXHtv+MjWHw0fpX9SuO2S208HdLK7yO7zhQ9G8VAvbnOYYCq7L2P7SfRpU1v9PLFTXHC0uN7Mon7m8IxfJ++NyrbBw2bdJE7XtqO1s37f1v4p9t/hsvVmx/eYDcjl/pMv7Or76edpXleGzItfvHz7F4sbHtjrG6jqOhx2kaOKbEJFpc5AVyd+/ZXfXt31eNv3ZcEjOFFtqt27RSV08ap9uf/8PRJG5ee/Tq4eX7+798pYYMHaw6duqgRl99pWrUqJGaf+fNbN954RpVUvpCZVi3xOlhxYXLiu0v6330P+dc1udu+9RrUYWOdvOWLdXg4SPUGQMHFbbfpWMo2NB/zrALk/XufGyFbo+aMDHyum1bHUfR1LFL10KBd4aOjxg9Rq3cuVv16NtPL3fv01fNvuNu7xhMoQW69urtHZvbbtS4sTr7/AtUt969k/2VIs5Ztzx/tcem4GKDOwRaABi6duiUfKNC/69efEc9XH63atuqjV4eevZ5+vWnjz+r18EyZsOGDRysmjVtptq0aq3jN4y/TnVs115d850r1fjvfFfHBg8YpLedcOloPX5vmfx9He/RuatqUvAPH8zoTzuuUOI8dcvz1iu4WM4sUJ9pG7Ru2crrg9+mb2DfAfq1d9ceqlvHLrq9+aG1er0510xTrVu0VONGXKHjn27am+Rp1qSpBjNoWDYf1stvuUf16tqdHWco5vhgVElV1BzwPa6eFvewnMV4a/VmvS+8uvGn7l6ejL2+3Xupdq3b6DgKM6w//NwLtD9mnGF8dW7fUfXs3M37O3Rq10F7grGP+LK5t+n4fXMW6uWuHTsn54sL1s+ff0O3MbZN3OQyywDHg7+l6cPf7/HKpez8ihHnwcW+qswUWAYjUwgYbMHF8uVBWl7ERg0b4Y3D2eOn6jbe27jlivVwYUeMvv8RA4P6nak/Q7q076TjeIQAcbwWK47rQ7zf0qL+RgWXKWhZzqwgFy200sY2xio+V7E+OKNnn8SjpXMWeLnoOF08e36Sw2wP3HGK8Z3HecW5tVhxkRfI3b1HN7Vu01rdNjFTaC2rWqKO/T5qt23XVj3zYnWyjnl1j8+NH/7qkG4/uWGVata8Gdt3XpT1PKswqGovK6uovqiscu15hW9TA8rKV/fU36rKq1oaE/0PUf1GXxj3sOKCYs5z5PgJLG7abqHlxrd8cESt3/eealMYdBNn36hjoydP1a+du3dXMytv020UZWY7FE1NCh+2aM9ccLuXjx6vuw3aFQ8uT+Krdr3qrYOZLNv+iOVKQ5/7/OoxdUL91ehiiw3uEEwBMOWKqzXvPxPNPK2962Ed3/fkNv2KmLnwXzRoiF7GRWniyLG6rc91ebVu4xkjsw0uWjuranX7jZUbk3X7dOup26vvfEh/qJr4nlWb9LcrfGDX99xMTuZlGszbZPaQ5cyCW2jh/Ny8aOMiZIog6jc+ZE0c4EMQxeq1l1+l42sKHpoLEAphc7vArI+/6eH1r+mL2gsPrNJ9lwy+MJdbb2YfZSVV1bjg57vM44raOWWB/mJ9FON712xJYm7BhNvept2hTTtVPmmmbsNzFF6YGUA/Zq1MvjdXvqDbuPDAI/iLi9KZvftFedq2V9PHTtKzXOvuqUr+BshlvmigWHDPBW0Uah9veF3NGjc5ec+YPlMMnwhYPy6y+GeAGb8VNX9kfVGxxfLlQVpec4EH7jikM1oPld9V9P0/5MxBuv3qT4p/Mckb7a/7WUBV2t9ci1nkcgstM17Rdse2KZLMGDPbpuWi4/Scfmel5nDH6daHn1LNm0V/l/qA/MZGVlzkBXKj0DJtU3C5tw73HXlL/fCx+9Sg8waqJQ/fk6yL16rV0cXu7/70pVq+8iF17vnnJP03lM9KaNBzmLoIA+pogWOFwfWbwhv+94XXPxWW/73wQfnfbPD5A/Eojo0WF2ngNl7b9lGF7hYhpu0WWi1atVLtO3VWT+/dp5dnVEbfOGlOxK6adr1myMXRhRtxFE1Vm7am7ofmcbdx+0aMHqvGz5jlrTP79rv0cZl90lxpYJ2CTx/WCfVWUygGyur3JjcFwGsrNmjwbdP0tWzeQvfhdgCWceF3nwNAUYT+tG8/WEaRduM103T7+WUrvb7Rf3WZ+v7VUzRm27Qc7nIW6u1v4yYsZxboM1ovxcUnoOdplk2hZeL9Ct9oN/4wKpSeWfK4GtDL3kZD4XH/D27TM1/mGyu2df9euKjRWbP6Evub4p3n4/8VPiv+l8Urav6c5ThwftgOXwzMMWD2o67xBO6eebMu4s0y1jcFKy485iKHMV8sD5bNLRjM9D46/149q+iul7YN9o1imvbVhfaYzmRZD+PxW/NfrA8E7utEKZYXRT8dh7TQwuxi6Puf9uWJ9tf9LKAq5W/FuullOR4bcrmF1tIbF6SObXesutvSZbOOO05NAWUKLbqN6cNnjNuXBeQzNrLiIi+Q2xRaB37xjtmpd+vwzLMHqHc+3qOumjCGFVoA8ZlzrtfftA4fP6h+8+fjuv/4H36Z8PUfP2f7zgvXqJLyBmDNXqeHFRelwPqLnlidtE3cLbTAkjXVerniwUfUtXPmen1uLswsuSCeV6GF/brrTK9YqM4fcQnbXymwbWRTHfL9ze3WFr2l5YKZA/S9/kQ0E0ULLcxQoR/Pv9AcWD749Mu6/eGzr6j2bdrqqW9TlOEbE76hAWyP2YK0HPSYQkGO2KfS8j48k9lYiOXMQimfadws00ILhZUptDbcvyIptC4ferFu45ba2OEjixZaJobnvXA7gh5HFpAvprQqa3Y6Hv9nWeXKDnEPy3ki0ELosguGJ+PJ/JgjLTdm8cYUPDLLs8ZNUd8b/de67V68zMPEaXmwjHG840c1uo0ZW4xndz26DT6/Ebt+zER9O93tqwtsp53Sz2aRmS0jemsrKgAgli8P0vIihkKLjkO30IJvuJVN3/9pOd1l2pcnyB17VVyuvw34HBxyuQUUivO0sR1SaNFxWqrQwjjFF24azwryJCbR4iIvkNsUWqBd+3Y6hkLry3/6TLfXrl+p+9q0aZ1aaH31z78yB+vlReFmlre88jzbd164RpVUVABsoeGCWHHhMvTSywpF1B7dNrfe8KC7KULwUPnStetU8xbR7Arij23epl/HTLlOzxzhXzA0LryZp95UruPm1iFmvqbMjWLgiR0v61f31iH6mzRtmqyDfWx+/2Nv2WyDNp77MvGaPW9762CGDe1tH//C218psL4xqqS0v+sW5P2wdrECALcB8EbEbRXTby78uO2EZTyzghkrtBHHDBfaeOjSbLP8lkX6FQ/DmxhecbvH7MvczkIcD8LiwwTPb6QdVyjIYYwqKeMvF8uZhWI+Axo3yydaaKFoWnXHg7qNGatShZZ5lmPl7Q948awgV0xpVdRsLnw+/E/Z/Nq+pIflLMYVwy5R79VEPwr43ugJ+oKNNmZIkMc8wG6e4cIXhcppc3QbhdOtU2Yn49D9wYV5KLtYoYVC4Y7p83Qbz8UhjoJh/tTZycP5+FWcey70vOA3YjR+IsTb+TIFl13+nS4AbIFlxPLlAc2LGRM35o5DfE48MO/OpK/82plF3//F9kH78gS5XcNSZfz1CywjljMryOUWUPgBBmJ0bBcrtNwffZhcdJyWKrTMOM3jX+sA5EpMosVFXiB3n369k+Vf/8sX0Yn96zG9fMnI6HYWLvoospY9ujTZzs2DKnPqjMnJ8pGvP9DPZcUnocoXzGX7zgvXqIxixYXLT7btLJxL8+RcTCEDps67WcfwgDt+TYgHzRHv2rOXjrv/vwoFGQorxN1nvbr06KFj2BZ5EEPRdNO99+l4y9ZtksIITJ47T8cxO4Vlc/wdunTRBZl5gH7yD25KtnHPceHyx5JzMfsrRbxufcQGdwiYbUIOg545Xf+abps3Nx5kxUUKF37MnsybNEP3n9d/YJIHFzA8a4E4HoT/ZOMeHXcfFDa/2sJshBu/aeJ0Hd/16NN6GVPli2bdotv0eEOJ91EfsZxZwOyfKQwodB9muWr+vV4fbsdsiR+A3/zgGr2MNooPrIfCAhc3PLxt8tCfxWOWke6vPiBXTGlVVrWmoVgsZzFQ9OMigW3M81MGFJqIAzz8jhhu05kiDGMT4xpx3EIx69YufjTJgYfUzewtfvWG9wLauJ1u/qXE8HOH6veDeW9grCO+YuEy78cFaec1Y+yk1HhdxMeaVSxfHhhfDYiNGDxMt+k4xG1txM2D3CDt/U+P1V02XxBO2sPwpcVyZgXj1DzLakgb2+5YNWC2FutgFgzLaJt13HGKv53JkXbsiJX6BXQI8XFHJtHiQrC4RmUUKy5ONvTW4ckE/lDDAsUGd0NBbx2eDsAfaligWM7TGTw/gxkGGs8K/InJKpbz2wgevscMAh7Ip311AY+oaQFi+QQfeERNCxTLebqCcZrn+SBXYhItLgSLa1RGseLiZIPZqcdf3MHiJwP4Qw0LFBvcDcWPbl2sxl58OYufysAfaligWM7TGZyP+ZVdHiBfTFaxnN9G8E87s54rtqOmBYjlE3zgETUtUCzn6QrGKR7boPGswJvEJFpcCBbXqIxixYVggT/UsECxwS1Y4A81LFAsp2CBPzFZxXIKPvCImhYglk/wgUfUtECxnEIEvElMosWFYHGNyihWXAgW+EMNCxQb3IIF/lDDAsVyChb4E5NVLKfgA4+oaQFi+QQfeERNCxTLKUTAm8QkWlwIFteojGLFhWCBP9SwQLHBLVjgDzUsUCynYIE/MVnFcgo+8IiaFiCWT/CBR9S0QLGcQgS8SUyixYVgcY3KKFZcCBb4Qw0LFBvcggX+UMMCxXIKFvgTk1Usp+ADj6hpAWL5BB94RE0LFMspRMCbxCRaXAgW16iMYsWFYIE/1LBAscEtWOAPNSxQLKdggT8xWcVyCj7wiJoWIJZP8IFH1LRAsZxCBLxJTKLFhWBxjcooVlwIFvhDDQsUG9yCBf5QwwLFcgoW+BOTVSyn4AOPqGkBYvkEH3hETQsUyylEwJvEJFpcCBbXqIxixYVggT/UsECxwS1Y4A81LFAsp2CBPzFZxXIKPvCImhYglk/wgUfUtECxnEIEvElMosWFYHGNyihWXAgW+EMNCxQb3IIF/lDDAsVyChb4E5NVLKfgA4+oaQFi+QQfeERNCxTLKUTAm8Sk7W9uEYrgGpVR6oHa54QiwB9qWKDUiw8/KRQB/lDDAsVyChb4E5NVLKfgA4+oaQFi+QQfeERNCxTLKUTAm/8HydO0UgqplwUAAAAASUVORK5CYII=>