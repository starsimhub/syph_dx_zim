# Estimating the value of novel syphilis diagnostics in Zimbabwe: how much overtreatment can be avoided?

Robyn M. Stuart1,⤉,, Michael Marks2,3,4, Chido Dziva Chikwari5,6, Alina Muellenmeister7, Romesh Abeysuriya7, Griffins Manguro1, Remco P.H. Peters8, Darcy W. Rao1, Lori M. Newman1

1\. Gates Foundation, Seattle, Washington, USA  
2\. Clinical Research Department, Faculty of Infectious and Tropical Diseases, London School of Hygiene & Tropical Medicine  
3\. Hospital for Tropical Diseases, University College London Hospital  
4\. Division of Infection and Immunity, University College London  
5\. MRC International Statistics & Epidemiology Group, Department of Infectious Disease Epidemiology, London School of Hygiene & Tropical Medicine, London, UK  
6\. Biomedical Research and Training Institute, Harare, Zimbabwe  
7\. Burnet Institute, Melbourne, Australia  
8\. Department of Global HIV, Hepatitis and STI Programmes, World Health Organization, Geneva, Switzerland  
⤉ Contractor on assignment

# **Abstract**

# 

# **Introduction**

The global burden of syphilis, a bacterial sexually transmitted infection (STI), is both sizable and increasing. The World Health Organization (WHO) estimates around 8 million adults were newly infected with syphilis in 2022 [(1)](https://www.zotero.org/google-docs/?O8yC4u). In addition, nearly 700,000 cases are transmitted from pregnant women to their babies [(1)](https://www.zotero.org/google-docs/?DNGHWo); these congenital syphilis cases make up essentially all of the estimated mortality and morbidity burden of syphilis and approximately 85% of the overall global disability-adjusted life year (DALY) burden of all STIs [(2)](https://www.zotero.org/google-docs/?PzQKoy). Much of this burden results from persistent challenges in the accurate and timely diagnosis and management of syphilis [(3)](https://www.zotero.org/google-docs/?KMd4Ow). Current diagnostic approaches in resource-limited settings, while improving access to treatment, frequently result in substantial overtreatment. Syndromic management of genital ulcer disease (GUD), recommended by WHO for settings without laboratory capacity, treats all patients presenting with genital ulcers presumptively for multiple pathogens including T pallidum. When implemented as per the guidelines, this approach has sensitivity of 100% but specificity of 0%, meaning that many patients without syphilis receive unnecessary treatment. Studies across sub-Saharan Africa have documented that syphilis accounts for only 10-30% [(4–7)](https://www.zotero.org/google-docs/?YrM7sr) of GUD cases in many settings, with herpes simplex virus being the predominant cause, yet all patients receive anti-syphilitic therapy under syndromic management protocols.

Similarly, screening programs for pregnant women and key populations face diagnostic limitations. The widespread adoption of dual HIV-syphilis rapid diagnostic tests (RDTs) in antenatal care (ANC) has dramatically improved screening coverage across low- and middle-income countries. However, most point-of-care syphilis tests are treponemal antibody-based assays that cannot distinguish between active infection requiring treatment and previously treated and resolved infections. For sex workers who undergo repeated screening over time, the combination of more frequent screening and higher prevalence may lead to repeated unnecessary treatments. This overtreatment has associated opportunity costs, potential adverse reactions, and possible contribution to antimicrobial resistance concerns. Use of benzathine penicillin for these cases can also contribute to stock-outs, which have been a significant global concern over the past 5-10 years [(8–10)](https://www.zotero.org/google-docs/?cBdErb).

The development of novel point-of-care diagnostics offers potential solutions to reduce this overtreatment burden. For symptomatic presentations, rapid tests capable of detecting treponemes directly from genital ulcers could improve the specificity of GUD management. For antenatal screening, point-of-care tests that confirm active infection following a positive treponemal screening test could substantially reduce unnecessary treatment while maintaining the benefits of universal screening coverage.

While the potential for novel point-of-care diagnostics to reduce syphilis overtreatment has been explored in several studies, findings have been mixed. A recent modeling study by Zhang et al. compared the performance of different rapid diagnostic test algorithms among pregnant individuals and men who have sex with men, finding that adding a non-treponemal confirmatory test to standard treponemal screening (T/NT-RDT) could reduce overtreatment by 90% compared to treponemal-only testing, though at the cost of missing more active cases (125 versus 68 per 100,000) [(11)](https://www.zotero.org/google-docs/?Y91T4G). Earlier modeling work suggested that dual point-of-care tests detecting both treponemal and non-treponemal antibodies could reduce overtreatment while saving costs in high-prevalence settings, though these benefits were contingent on disease prevalence and loss to follow-up rates [(12)](https://www.zotero.org/google-docs/?1E7F1K). However, real-world implementation has not always matched these theoretical projections. In a field evaluation in Burkina Faso, Langendorf et al. found that the Dual Path Platform (DPP) screen-and-confirm assay, which detects both treponemal and non-treponemal antibodies, did not reduce overtreatment compared to standard treponemal testing (0.0% vs 2.5%; p=0.218) and had a significantly higher rate of underdiagnosis (48.4% vs 2.2%; p\<0.001) [(13)](https://www.zotero.org/google-docs/?fkdFwD). These findings underscore the context-dependency of diagnostic test performance and the importance of field evaluations in diverse epidemiological settings before widespread implementation of novel testing strategies.

Zimbabwe, with estimated active syphilis prevalence of 0.9% among all adults and 2.9% among HIV+ adults [(14–16)](https://www.zotero.org/google-docs/?2wM2JR), represents an important setting for evaluating the potential impact of improved diagnostic approaches. In this study, we use STIsim, an agent-based model of co-transmitting STIs, to simulate both syphilis and HIV transmission to estimate the extent of overtreatment under current diagnostic algorithms and quantify the potential reduction achievable through implementation of novel point-of-care confirmatory tests for both symptomatic testing and key population screening.

# **Methods**

## **Model overview**

We used STIsim, an agent-based model of co-transmitting sexually transmitted infections (STIs) built on the Starsim disease modeling platform [(17,18)](https://www.zotero.org/google-docs/?jqD32c). We initialized a population of 10,000 agents representative of Zimbabwe in 1985, scaled to the national population using UN estimates of age-specific fertility, age- and sex-disaggregated mortality, and overall population size [(19)](https://www.zotero.org/google-docs/?49JC7U). We created an age- and risk-structured sexual network with three risk groups (low, medium, and high) governing partnership formation rates and concurrency. In addition, we modeled participation in transactional sex: approximately 10% of women in the highest risk group engage in sex work at any given time, with 20% of men seeking transactional sex partners. Pregnant women reduce sexual risk behavior, reverting to their pre-pregnancy risk profiles postpartum.

## **Disease models**

We modeled HIV transmission, progression, diagnosis, and treatment as previously described [(17)](https://www.zotero.org/google-docs/?rqEPrN). Briefly, HIV-infected individuals progress through acute, chronic, and AIDS stages with stage-specific infectiousness. We modeled the scale-up of ART in Zimbabwe from 2004 onwards, with ART reducing transmission probability by 96% among virally suppressed individuals. HIV testing, diagnosis, and ART initiation were modeled separately for sex workers (higher testing rates), individuals with low CD4 counts, and the general population.

We simulated syphilis (*Treponema pallidum*) transmission and natural history. After acquisition, infections progress through an incubation period (mean 21 days), primary infection characterized by chancres (mean duration 6 weeks; visible in approximately 40% of cases), secondary infection with disseminated rash (mean duration 3.6 months), and subsequent latency. The first 12–14 months of latency are classified as early latent; thereafter, infections transition to late latent. Upon entering latency, 35% of individuals are assigned to reactivate to secondary syphilis, with the timing of reactivation drawn from a lognormal distribution (mean 1 year, SD 1 year from entry to latency). The remaining 65% do not reactivate; of these, 35% eventually progress to tertiary syphilis (mean 20 years).

Sexual transmission probability is stage-dependent: primary syphilis is modeled as 5-fold more infectious than secondary, while latent infections have 10% of secondary-stage infectiousness that decays exponentially with a half-life of 6 months. This parameterization produces a transmission profile where approximately 60% of sexual transmissions originate from primary-stage infections, 30% from secondary, and the remainder from early latent — consistent with the epidemiology of a disease sustained by rapid partner turnover during the short but highly infectious primary window.

For maternal-to-child transmission (MTCT), we used stage-independent transmission probability, reflecting that spirochetemia during pregnancy affects the fetus regardless of clinical stage. Birth outcomes among MTCT-affected pregnancies depend on maternal stage at delivery: among mothers with primary or secondary syphilis, 10% of pregnancies result in neonatal death, 20% in stillbirth, 45% in congenital syphilis, and 25% in normal outcomes. Early latent mothers produce somewhat better outcomes (5% neonatal death, 10% stillbirth, 40% congenital, 45% normal), while late latent mothers have predominantly normal outcomes (10% stillbirth, 10% congenital, 80% normal).

We modeled bidirectional interactions between HIV and syphilis via a connector module. HIV infection increases susceptibility to syphilis acquisition (calibrated median relative susceptibility 2.5; Table S6), while active syphilis increases susceptibility to HIV acquisition (calibrated median relative susceptibility 1.9; Table S6).

## **Diagnostic and treatment pathways**

We modeled five pathways through which syphilis infections are detected and treated under the standard of care (Table 1). Because the treponemal-based screening tests used in the ANC, KP, and PLHIV pathways detect antibodies that persist after successful treatment, previously infected individuals who have cleared their infection will continue to test positive. This drives substantial overtreatment, particularly among key populations with high cumulative exposure. Care-seeking probability for the GUD syndromic pathway varies by sex, risk group, and year. For ANC screening, positive results lead to treatment of both the mother (with benzathine penicillin) and — through in-utero treatment — the fetus; newborns of ANC-positive mothers are scheduled for clinical examination at delivery.

**Table 1. Diagnostic and treatment pathways under the standard of care.**

| Pathway | Eligible population | Test | Sensitivity: active syphilis | Sensitivity: past/no infection | Action if positive |
|---|---|---|---|---|---|
| GUD syndromic | Individuals with visible genital ulcers (primary chancre or other GUD) | Syndromic management [(20)](https://www.zotero.org/google-docs/?JsKrfG) | 80% (all stages) | 80% (all states) | Presumptive treatment |
| ANC screening | Pregnant women (coverage scaling to ~50%) | Dual HIV-syphilis RDT (treponemal) | 95% (secondary, latent); 20% (primary, exposed) | 95% (past infection); 1% (naive) | 70% treated; 30% no action |
| KP screening | Female sex workers accessing HIV testing | Dual HIV-syphilis RDT (treponemal) | As above | As above | Treatment |
| PLHIV screening | PLHIV on ART (scaling from 2019) | Dual HIV-syphilis RDT (treponemal) | As above | As above | Treatment |
| Newborn exam | Newborns of ANC-positive mothers | Clinical examination | 30% (congenital) | 10% (susceptible) | 10% direct treatment; 10% exam; 80% no intervention |

## **Intervention scenarios**

We evaluated four novel diagnostic scenarios introduced from 2028 onwards, compared with the current standard of care (Table 2).

**Table 2. Intervention scenarios.**

| Scenario | Pathway affected | Change from SOC (2028 onwards) | Test sensitivity: active syphilis | Test sensitivity: past/no infection |
|---|---|---|---|---|
| SOC (baseline) | — | — | — | — |
| GUD POC test | GUD syndromic | POC *T. pallidum* detection replaces syndromic management | 95% (primary); 20% (secondary); 5% (other) | 5% (all states) |
| Confirmatory test | ANC, KP, PLHIV | POC confirmatory test after treponemal screen-positive | 95% (primary, secondary, latent) | 0% (naive); 0% (past infection); 10% (exposed) |
| Both | GUD + ANC/KP/PLHIV | Combined GUD POC test and confirmatory test | As above per pathway | As above per pathway |
| Newborn POC (CS) | Newborn | POC clinical screening replaces standard exam | 90% (congenital) | 0% (susceptible) |

## **Model calibration**

We jointly calibrated HIV and syphilis parameters to Zimbabwe-specific epidemiological data using Bayesian optimization with the Optuna framework [(21)](https://www.zotero.org/google-docs/?optuna), implemented via the STIsim Calibration class. We fit 15 parameters simultaneously: 4 governing HIV transmission and treatment (transmission probability, condom efficacy, initial prevalence scaling, effective ART duration), 3 governing syphilis transmission (transmission probability, condom efficacy, relative primary-stage infectiousness), 3 governing sexual network structure (proportion of women and men in the lowest risk group, male mid-risk concurrency), 2 governing HIV-syphilis coinfection interactions (relative susceptibility in each direction), and 3 governing diagnostic testing rates (symptomatic care-seeking, ANC screening, and KP screening rate multipliers). Prior ranges and posterior estimates for all parameters are provided in Table S6.

Calibration targets included HIV prevalence among adults aged 15–49, people living with HIV, new HIV infections, HIV-related deaths, and ART coverage from UNAIDS estimates (2000–2024); active syphilis prevalence, active syphilis prevalence by sex, and active syphilis prevalence by HIV status from ZIMPHIA 2015-2016; and syphilis new infections and active case counts from Spectrum-STI estimates. We ran 2,000 Optuna trials. We applied higher calibration weights to syphilis-specific ZIMPHIA survey targets (weight 10–20) relative to HIV time series targets (weight 2–5), reflecting the greater precision and direct relevance of the cross-sectional survey data. To prevent bimodal calibration outcomes in which syphilis fit well but HIV prevalence collapsed, we implemented a post-simulation quality filter rejecting parameter sets where median HIV prevalence in the final 5 years fell below 5% or where syphilis new infections reached zero (epidemic extinction). We also applied a pre-simulation filter rejecting parameter combinations where the effective syphilis transmission force (product of beta, relative primary infectiousness, and 1 minus condom efficacy) was too low to sustain the epidemic.

## **Model analysis**

We selected the top 200 best-fitting parameter sets from the calibration and re-ran each as a single simulation to produce comprehensive epidemiological and treatment outcomes across all model outputs. We computed percentile statistics (median, 10th–90th percentiles) across the parameter sets, used for figures and uncertainty characterization. For the scenario analysis, we ran the top 10 parameter sets (sorted by effective syphilis transmission force to prioritize sets with sustained epidemics) forward from 1985 to 2041 under each of the five diagnostic scenarios, using a single random seed per parameter set. Parameter sets where syphilis transmission died out were excluded.

For each scenario, we computed: (1) the overtreatment rate, defined as the proportion of treatments given to individuals without active syphilis (including those with serological scarring from past infection and true false positives); (2) the number of correct treatments (individuals with active syphilis who were treated); and (3) the number of missed active infections per care-seeking pathway. We disaggregated outcomes by diagnostic pathway (GUD syndromic, ANC, KP, PLHIV, newborn) and, for the newborn pathway, by whether the treatment correctly identified congenital syphilis. For congenital syphilis, we tracked maternal transmission events and birth outcomes (death, congenital syphilis, normal) disaggregated by the mother's disease stage at delivery. We report results as medians across the calibrated parameter sets.

To characterize bottlenecks in the current standard of care, we constructed two care-seeking cascades. For the GUD syndromic pathway, we estimated the proportion of new primary syphilis infections that progress through each step of the care-seeking cascade: developing a visible chancre (determined by the sex-weighted probability of symptomatic primary infection, approximately 38%), seeking care during the primary stage, testing positive under syndromic management (80%), and being correctly treated (98% treatment efficacy). For congenital syphilis prevention, we estimated the cascade of opportunities to prevent adverse birth outcomes among pregnancies where the fetus is perinatally exposed to syphilis: ANC attendance (approximately 90% in Zimbabwe), syphilis screening at ANC, testing positive on the dual treponemal test, maternal treatment, successful in-utero fetal treatment (98% efficacy in the first two trimesters, 75% in the third trimester), and — for cases not prevented prenatally — newborn detection and treatment. We estimated the total number of at-risk pregnancies as the sum of observed maternal-to-child transmission events (from unscreened or unsuccessfully treated mothers) and estimated prevention events (from successfully treated mothers, assuming approximately 90% average fetal treatment efficacy).

All analyses were conducted using Python 3.11 with STIsim v1.5.0 and Starsim v3.1.1. Code is publicly available at https://github.com/starsimhub/syph_dx_zim.

# **Results**

## Model calibration to syphilis and HIV epidemiology

Figure 1 summarizes the calibrated model's fit to syphilis and HIV epidemiology in Zimbabwe. The model reproduces key epidemiological patterns across both diseases and is broadly consistent with survey data from the Zimbabwe Population-based HIV Impact Assessment (ZIMPHIA) surveys conducted in 2015-2016 and 2020. Women engaged in transactional sex have substantially higher prevalence than the general female population across all age groups, consistent with their elevated exposure through higher partner turnover and concurrency (Figure 1A). Overall active syphilis prevalence among adults aged 15-50 is approximately 0.7%, consistent with the ZIMPHIA 2016 estimate of 0.8% among women and 0.6% among men.

Figure 1B shows modeled HIV prevalence by age and sex across four time points: pre-ART (2005), 2016, 2020, and 2025 (projected). The model reproduces the characteristic age-sex pattern of HIV in Zimbabwe, with peak prevalence among women aged 30-35 and men aged 35-50. The epidemic shows a marked aging pattern over time: while prevalence among younger age groups has declined substantially between 2005 and 2025 (reflecting reduced incidence due to treatment-as-prevention and behavior change), prevalence among older adults has remained high or even increased, reflecting the success of ART in extending survival. The aging HIV epidemic has direct implications for syphilis-HIV coinfection dynamics, shown in Figure 1C. Our modeling analyses estimate that the ratio of syphilis prevalence among HIV-positive versus HIV-negative adults has been declining over time, from an estimated peak of approximately 14x in the early 2000s to around 5x by 2020. Although HIV increases susceptibility to syphilis acquisition, the declining ratio reflects the fact that a growing share of HIV-positive individuals are older adults on ART who are less sexually active and therefore less exposed to syphilis.

For syphilis, women in transactional sex account for approximately 28% of new acquisitions and 13% of onward transmissions, while men engaged in transactional sex account for 24% of acquisitions and 36% of transmissions. This asymmetry reflects the epidemiological importance of male clients as a bridge population, who transmit syphilis acquired during transactional sex to their other sexual partners. A similar though less pronounced pattern is seen for HIV (Figure 1D). 

We estimate that sexual transmission of syphilis is dominated by primary syphilis (62%), followed by secondary (31%) and early latent (8%) stages, reflecting the high infectiousness of early-stage disease (Figure 1E). Maternal transmission shows a different pattern: late latent syphilis contributes the largest share of adverse congenital outcomes (62%), because although per-pregnancy transmission risk is lower for late latent infections, the much larger pool of women with latent syphilis means more maternal transmissions originate from this stage. Early latent (11%) and primary (13%) stages also contribute. Deaths (stillbirths and neonatal deaths) and congenital syphilis cases are shown separately within each stage, with an estimated 1,477 deaths and 1,715 congenital syphilis cases annually.


**Figure 1. Syphilis and HIV epidemiology in Zimbabwe.** (A) Active syphilis prevalence by age group and sex work status, with ZIMPHIA 2016 survey data (diamonds). Dark pink bars show prevalence among women in transactional sex; light pink and light blue bars show overall female and male prevalence, respectively. The rightmost group shows aggregate prevalence across ages 15-50. (B) HIV prevalence by age group and sex (pink = female, blue = male) at four time points, with ZIMPHIA 2016 (diamonds) and 2020 (squares) survey estimates and 95% confidence intervals. Model uncertainty bands show the 10th-90th percentile range across calibrated parameter sets. (C) Ratio of syphilis prevalence among HIV-positive versus HIV-negative adults over time, with ZIMPHIA 2016 estimate (diamond). (D) Distribution of syphilis and HIV infections acquired and transmitted by sex work status (2000-2019 average). (E) Distribution of syphilis transmission events by disease stage. Left: share of sexual transmissions by stage. Right: share of adverse maternal outcomes (deaths and congenital syphilis) by mother's disease stage at delivery.

## Care-seeking cascades under the standard of care

Figure 2 illustrates the care-seeking cascades for the two primary pathways through which syphilis is detected: syndromic GUD management and antenatal screening for congenital syphilis prevention.

The GUD syndromic cascade reveals that the vast majority of primary syphilis infections go undetected (Figure 2, left). Of every 100 new primary syphilis infections, only 38 develop a visible chancre — the remainder occur at internal or painless sites that do not prompt clinical attention. Of those with visible ulcers, fewer than 3% seek care during the approximately 6-week primary window, reflecting low care-seeking rates for genital ulcers in the general population. Of those who do present, syndromic management correctly identifies 80% as requiring treatment. The net result is that fewer than 1 in 100 primary syphilis infections are correctly treated through the GUD pathway.

For congenital syphilis prevention, the cascade shows that ANC screening coverage is the dominant bottleneck (Figure 2, right). Although approximately 90% of pregnant women in Zimbabwe attend at least one ANC visit, only about 25% of at-risk pregnancies result in the mother being screened for syphilis — reflecting the gap between ANC attendance and syphilis test availability or uptake. Among mothers who are screened, the dual treponemal test has high sensitivity for the latent infections that predominate in pregnancy, and nearly all test-positive mothers are treated. In-utero treatment is highly effective: of the 24 mothers treated per 100 at-risk pregnancies, approximately 22 fetuses are successfully cured. However, the 78 at-risk pregnancies where the mother was not screened or treated produce the bulk of congenital cases. Under the current standard of care, newborn detection is minimal: only about 1 in 100 at-risk births results in the newborn being identified and treated, as 80% of newborns receive no testing at all and clinical examination has only 30% sensitivity for detecting congenital syphilis.

**Figure 2. Care-seeking cascades under the standard of care.** Left: GUD syndromic management cascade per 100 new primary syphilis infections. Blue bars show the number retained at each step; grey bars show losses. Right: Congenital syphilis prevention cascade per 100 at-risk pregnancies (pregnancies where the fetus is perinatally exposed to syphilis). Blue bars show the ANC screening and treatment pathway; red bars show the residual burden at birth. Loss annotations indicate the reason and magnitude of attrition at each step. Numbers represent annual averages over 2015–2024.

## Treatment outcomes under the standard of care

GUD syndromic management accounts for the largest volume of syphilis treatments (~15,000 per year) but also the highest overtreatment rate: 92% of treatments through this pathway are given to individuals without active syphilis (Figure 3A–B). This reflects the fundamental limitation of syndromic management, which treats all GUD presentations presumptively regardless of syphilis status. Because background GUD prevalence (1%) far exceeds active syphilis prevalence (~0.6%), most people presenting with genital ulcers do not have syphilis but are treated nonetheless.

ANC screening treats approximately 4,000 women per year with an overtreatment rate of 35%, driven by the treponemal test detecting past cleared infections. KP screening through dual HIV-syphilis RDTs has a 46% overtreatment rate, while PLHIV screening has the second-highest overtreatment rate at 53%, reflecting the high cumulative syphilis exposure in these populations. Newborn treatments are predominantly overtreatment (81%), as most newborns flagged for treatment through the ANC pathway do not have congenital syphilis.

When syphilis is correctly detected and treated, the stage at which infection is identified varies markedly across pathways (Figure 3C). GUD syndromic management predominantly detects primary-stage infections (~80%), consistent with its reliance on visible ulcers as the trigger for care-seeking. In contrast, ANC, KP, and PLHIV screening pathways detect infections predominantly at late latent stage (~80%), reflecting that the treponemal screening test identifies long-standing infections that have progressed beyond the early infectious stages. This has implications for the epidemiological impact of each pathway: GUD syndromic management, while inefficient due to overtreatment, at least targets the most infectious stage of disease, whereas the screening pathways primarily detect infections that contribute relatively little to ongoing transmission.

**Figure 3. Treatment outcomes under the standard of care.** (A) Total syphilis treatments per year by diagnostic pathway, with correctly treated (green) and overtreated (red) stacked. Numbers above bars indicate total treatments. (B) Overtreatment rate by pathway, defined as the proportion of treatments given to individuals without active syphilis. (C) Stage of syphilis infection at the time of correct detection, by pathway (excluding newborn pathway). All values are annual averages over 2000–2024.

# **Discussion**

# **Competing interests**

The authors have no competing interest to declare.   
 

# **Ethics**

This study is a computer simulation study. This study does not involve human participants. No ethical approval is required.

# **References**

[1\.](https://www.zotero.org/google-docs/?xWpmcA)	[WHO. Data on syphilis \[Internet\]. \[cited 2025 Dec 9\]. Available from: https://www.who.int/data/gho/data/themes/topics/data-on-syphilis](https://www.zotero.org/google-docs/?xWpmcA) 

[2\.](https://www.zotero.org/google-docs/?xWpmcA)	[Vos T, Lim SS, Abbafati C, Abbas KM, Abbasi M, Abbasifard M, et al. Global burden of 369 diseases and injuries in 204 countries and territories, 1990–2019: a systematic analysis for the Global Burden of Disease Study 2019\. The Lancet. 2020 Oct 17;396(10258):1204–22.](https://www.zotero.org/google-docs/?xWpmcA) 

[3\.](https://www.zotero.org/google-docs/?xWpmcA)	[Pham MD, Ong JJ, Anderson DA, Drummer HE, Stoové M. Point-of-Care Diagnostics for Diagnosis of Active Syphilis Infection: Needs, Challenges and the Way Forward. Int J Environ Res Public Health. 2022 Jul 4;19(13):8172.](https://www.zotero.org/google-docs/?xWpmcA) 

[4\.](https://www.zotero.org/google-docs/?xWpmcA)	[Kularatne R, Venter JME, Maseko V, Muller E, Kufa T. Etiological Surveillance of Genital Ulcer Syndrome in South Africa: 2019 to 2020\. Sex Transm Dis. 2022 Aug 1;49(8):571–5.](https://www.zotero.org/google-docs/?xWpmcA) 

[5\.](https://www.zotero.org/google-docs/?xWpmcA)	[Mungati M, Machiha A, Mugurungi O, Tshimanga M, Kilmarx PH, Nyakura J, et al. The Etiology of Genital Ulcer Disease and Coinfections With Chlamydia trachomatis and Neisseria gonorrhoeae in Zimbabwe: Results From the Zimbabwe STI Etiology Study. Sex Transm Dis. 2018 Jan;45(1):61–8.](https://www.zotero.org/google-docs/?xWpmcA) 

[6\.](https://www.zotero.org/google-docs/?xWpmcA)	[Tshaka TR, Singh R, Apalata TR, Mbulawa ZZA. Aetiology of genital ulcer disease and associated factors among Mthatha public clinic attendees. S Afr J Infect Dis. 2022 Dec 7;37(1):444.](https://www.zotero.org/google-docs/?xWpmcA) 

[7\.](https://www.zotero.org/google-docs/?xWpmcA)	[Paz-Bailey G, Rahman M, Chen C, Ballard R, Moffat HJ, Kenyon T, et al. Changes in the etiology of sexually transmitted diseases in Botswana between 1993 and 2002: implications for the clinical management of genital ulcer disease. Clin Infect Dis. 2005 Nov 1;41(9):1304–12.](https://www.zotero.org/google-docs/?xWpmcA) 

[8\.](https://www.zotero.org/google-docs/?xWpmcA)	[Nurse-Findlay S, Taylor MM, Savage M, Mello MB, Saliyou S, Lavayen M, et al. Shortages of benzathine penicillin for prevention of mother-to-child transmission of syphilis: An evaluation from multi-country surveys and stakeholder interviews. PLoS Med. 2017 Dec 27;14(12):e1002473.](https://www.zotero.org/google-docs/?xWpmcA) 

[9\.](https://www.zotero.org/google-docs/?xWpmcA)	[WHO. Shortages of benzathine penicillin. How big is the problem? And why it matters \[Internet\]. \[cited 2025 Dec 9\]. Available from: https://www.who.int/news/item/26-12-2017-shortages-of-benzathine-penicillin.-how-big-is-the-problem-and-why-it-matters](https://www.zotero.org/google-docs/?xWpmcA) 

[10\.](https://www.zotero.org/google-docs/?xWpmcA)	[Seghers F, Taylor MM, Storey A, Dong J, Wi TC, Wyber R, et al. Securing the supply of benzathine penicillin: a global perspective on risks and mitigation strategies to prevent future shortages. Int Health. 2024 May 1;16(3):279–82.](https://www.zotero.org/google-docs/?xWpmcA) 

[11\.](https://www.zotero.org/google-docs/?xWpmcA)	[Zhang Y, Chow EPF, Zhang L, Fairley CK, Ong JJ. Comparison of the performance and costs of testing algorithms using rapid diagnostic tests for detection and treatment of syphilis among pregnant individuals and men who have sex with men: a modelling study. The Lancet Infectious Diseases \[Internet\]. 2025 Nov 21 \[cited 2025 Dec 9\]; Available from: https://www.sciencedirect.com/science/article/pii/S1473309925005882](https://www.zotero.org/google-docs/?xWpmcA) 

[12\.](https://www.zotero.org/google-docs/?xWpmcA)	[Owusu-Edusei K, Gift TL, Ballard RC. Cost-effectiveness of a dual non-treponemal/treponemal syphilis point-of-care test to prevent adverse pregnancy outcomes in sub-Saharan Africa. Sex Transm Dis. 2011 Nov;38(11):997–1003.](https://www.zotero.org/google-docs/?xWpmcA) 

[13\.](https://www.zotero.org/google-docs/?xWpmcA)	[Langendorf C, Lastrucci C, Sanou-Bicaba I, Blackburn K, Koudika MH, Crucitti T. Dual screen and confirm rapid test does not reduce overtreatment of syphilis in pregnant women living in a non-venereal treponematoses endemic region: a field evaluation among antenatal care attendees in Burkina Faso. Sex Transm Infect. 2019 Sep;95(6):402–4.](https://www.zotero.org/google-docs/?xWpmcA) 

[14\.](https://www.zotero.org/google-docs/?xWpmcA)	[Ruangtragool L, Silver R, Machiha A, Gwanzura L, Hakim A, Lupoli K, et al. Factors associated with active syphilis among men and women aged 15 years and older in the Zimbabwe Population-based HIV Impact Assessment (2015–2016). PLOS ONE. 2022 Mar 17;17(3):e0261057.](https://www.zotero.org/google-docs/?xWpmcA) 

[15\.](https://www.zotero.org/google-docs/?xWpmcA)	[Farahani M, Killian R, Reid GA, Musuka G, Mugurungi O, Kirungi W, et al. Prevalence of syphilis among adults and adolescents in five sub-Saharan African countries: findings from Population-based HIV Impact Assessment surveys. The Lancet Global Health. 2024 Sep 1;12(9):e1413–23.](https://www.zotero.org/google-docs/?xWpmcA) 

[16\.](https://www.zotero.org/google-docs/?xWpmcA)	[Korenromp EL, Mahiané G, Rowley J, Nagelkerke N, Abu-Raddad L, Ndowa F, et al. Estimating prevalence trends in adult gonorrhoea and syphilis in low- and middle-income countries with the Spectrum-STI model: results for Zimbabwe and Morocco from 1995 to 2016\. Sex Transm Infect. 2017 Dec;93(8):599–606.](https://www.zotero.org/google-docs/?xWpmcA) 

[17\.](https://www.zotero.org/google-docs/?xWpmcA)	[Stuart RM, Newman L, Manguro G, Chikwari CD, Marks M, Peters RPH, et al. Reduction in overtreatment of gonorrhea and chlamydia through point-of-care testing compared with syndromic management for vaginal discharge: a modeling study for Zimbabwe. BMJ STI.](https://www.zotero.org/google-docs/?xWpmcA) 

[18\.](https://www.zotero.org/google-docs/?xWpmcA)	[Institute for Disease Modeling. Starsim. 2024\.](https://www.zotero.org/google-docs/?xWpmcA) 

[19\.](https://www.zotero.org/google-docs/?xWpmcA)	[United Nations. World Population Prospects 2022\. Department of Economic and Social Affairs, Population Division; 2022\.](https://www.zotero.org/google-docs/?xWpmcA) 

[20\.](https://www.zotero.org/google-docs/?xWpmcA)	[World Health Organisation. Guidelines for the management of symptomatic sexually transmitted infections. 2021\.](https://www.zotero.org/google-docs/?xWpmcA) 