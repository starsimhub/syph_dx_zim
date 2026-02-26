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

## **Model details**

We used STIsim, an agent-based model capable of simulating multiple co-infecting sexually transmitted infections (STIs) including HIV transmitting across age- and risk-structured sexual networks [(17)](https://www.zotero.org/google-docs/?jqD32c). STIsim leverages the Starsim disease modeling platform [(18)](https://www.zotero.org/google-docs/?gVC0mR). We defined a population of agents representative of Zimbabwe in 1980 and used UN estimates of age-specific fertility, age- and sex-disaggregated mortality, and overall population size to project births and deaths over 1980–2030 [(19)](https://www.zotero.org/google-docs/?49JC7U). We created an age- and risk-stratified sexual network for agent interactions, with the population divided into three risk groups representing people who marry and do not have extramarital relationships, people who marry and divorce and/or have extramarital relationships, and people who never marry. In addition, we modeled participation in transactional sex. 

We modeled HIV transmission, progression, diagnosis, and treatment as previously described [(17)](https://www.zotero.org/google-docs/?rqEPrN). Briefly, HIV-infected individuals progress through acute, chronic, and AIDS stages with stage-specific infectiousness. We modeled the scale-up of antiretroviral therapy (ART) in Zimbabwe from 2004 onwards, with ART reducing transmission probability by 96% among virally suppressed individuals. 

We simulated syphilis (*T pallidum*) transmission and natural history. After acquisition, infections progress through an initial incubation period (mean 21 days), primary infection characterized by chancres (which go unnoticed in approximately 60% of cases), secondary infection with disseminated rash, and subsequent latency. Approximately 35% of early latent infections reactivate to secondary syphilis. Untreated latent infections can persist indefinitely or progress to tertiary syphilis (assumed in 35% of cases after a mean of 20 years). Detailed natural history parameters and supporting references are provided in Table S1 and in the Supplementary Materials.

We modeled bidirectional interactions between HIV and syphilis. HIV infection increases susceptibility to syphilis acquisition, accelerates syphilis progression, and increases infectiousness. Syphilis increases susceptibility to HIV acquisition, with the greatest effect during primary and secondary stages when genital ulcers are present (relative risk for primary syphilis: \[X\]; secondary: \[X\]; early latent: \[X\]). HIV-syphilis coinfection increases HIV viral load and therefore HIV infectiousness (multiplier: \[X\]) \[XX\].

Vertical transmission risk varies by maternal syphilis stage (primary/secondary: 70%; early latent: 40%; late latent: 10%) \[XX\]. Congenital syphilis outcomes include stillbirth (40%), neonatal death (20%), and long-term morbidity (40%) \[XX\].

## **Diagnostic pathways**

We modeled three primary pathways for syphilis detection and treatment:

* **Symptomatic detection (genital ulcer disease).** Individuals with primary or secondary syphilis develop genital ulcers, and among those with symptomatic disease, we assume that XX% seek care at STI or primary healthcare clinics. Under the current standard of care, clinicians follow WHO guidelines for syndromic management of genital ulcer disease (GUD), which recommend treating all GUD cases presumptively for multiple pathogens including syphilis [(20)](https://www.zotero.org/google-docs/?JsKrfG). This approach achieves 100% sensitivity but 0% specificity for syphilis, as only 15–25% of GUD cases in Zimbabwe are caused by syphilis, with herpes simplex virus being the predominant cause.  
* **Asymptomatic detection via antenatal care screening.** We modeled routine ANC attendance and syphilis screening consistent with Zimbabwe's antenatal care program. Approximately 90% of pregnant women attend at least one ANC visit \[X\], where they receive dual HIV-syphilis rapid diagnostic tests \[X \- TODO check algo\]. The treponemal component has sensitivity of 95% and specificity of 98% \[XX\]. These tests use treponemal antibody detection, which cannot distinguish between active infection requiring treatment and past/treated infection. All women testing positive receive treatment with benzathine penicillin under the current standard of care.   
* **Asymptomatic detection via key population screening.** Sex workers in Zimbabwe have elevated syphilis prevalence (approximately \[X\]% active infection, \[X\]% seropositive) \[refs\]. We assume that \[X\]% of sex workers access regular STI screening services, where they receive dual HIV-syphilis rapid diagnostic tests. As with ANC screening, the treponemal component cannot distinguish active from past infection. Given the high background seroprevalence in this population and repeated testing over time, this leads to substantial overtreatment of women with past rather than active infection.

## **Intervention scenarios**

We evaluated two novel point-of-care (POC) diagnostic interventions, implemented individually and in combination from 2026 onwards. We evaluate these strategies compared with the current standard of care with syndromic management of GUD and treponemal RDT screening without confirmatory testing.

Scenario 1: POC test for treponeme detection from genital ulcers. We modeled a rapid test capable of detecting T. pallidum directly from genital ulcer specimens among individuals seeking care for GUD. We assumed test sensitivity of 95% and specificity of 95%, and we assumed that this test would be integrated into the syndromic management approach for GUD. Only test-positive individuals receive syphilis treatment.

Scenario 2: POC confirmatory test for active syphilis in sex workers. Following a positive treponemal screening test, we modeled a confirmatory rapid test to distinguish active infection from past infection. This could be based on non-treponemal antibody detection (analogous to rapid plasma reagin) or direct antigen detection. We assumed sensitivity of 90% for detecting active infection and specificity of 85% for ruling out past infection. Only women with both a positive screening test and positive confirmatory test receive treatment.

## **Model calibration**

We calibrated the model to Zimbabwe-specific epidemiological data using Bayesian optimization with the Optuna library. Calibration targets included: HIV prevalence (15–49 years): 12.8% in 2015, declining to approximately 11.5% by 2020 (ZIMPHIA, UNAIDS); HIV epidemiological trends over time including people living with HIV, new HIV infections, and HIV-related deaths (UNAIDS); active syphilis prevalence: 0.8–0.9% in general population (ZIMPHIA 2015-2016); syphilis prevalence among pregnant women at ANC: 2.0–4.3% over 2000–2022.

We fit model parameters governing syphilis transmission probability, symptomatic presentation rates, and care-seeking behavior. For each calibration iteration, we computed a mismatch score based on weighted mean squared error between model outputs and calibration targets. We selected the best-fitting 100 parameter sets and ran each for 10 stochastic realizations to characterize uncertainty. Additional calibration details are provided in the Supplementary Materials.

## **Model analysis**

For each scenario, we projected syphilis and HIV dynamics from 2025 to 2040 using the 50 calibrated parameter sets with 10 stochastic realizations each (N=500 simulations per scenario). We focused on two primary outcomes over 2025–2040: (1) overtreatment, defined syphilis treatment given to individuals without active syphilis infection including treatment of individuals with past/cleared infection (seropositive but not currently infected) and individuals without any syphilis infection history (false positives);  (2) undertreatment, defined as care-seeking individuals with active syphilis who were not treated due to false-negative test results. We summarize results as medians and 95% uncertainty intervals (2.5th–97.5th percentiles) across the 500 simulations per scenario.

All analyses were conducted using Python 3.11 with STIsim v1.5 and Starsim v3.1. Data and code are publicly available via Github.

# **Results**

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