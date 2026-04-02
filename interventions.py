"""
Define syphilis diagnostics
"""

import numpy as np
import pandas as pd
import sciris as sc
import starsim as ss
import stisim as sti


def load_syph_dx():
    """
    Create default diagnostic products
    """
    df = sc.dataframe.read_csv('data/syph_dx.csv')
    dxprods = dict()
    for name in df.name.unique():
        dxprods[name] = sti.SyphDx(df[df.name == name], name=f'SyphDx_{name}')
    return dxprods


def load_syph_products():
    """
    Create syphilis screening / testing algorithms
    """
    df = sc.dataframe.read_csv('data/syph_products.csv').fillna(value = 0)
    dxalgos = dict()
    df['algo_scen'] = df['algorithm'] + '_' + df['scenario']
    excl_cols = ['algorithm', 'scenario', 'algo_scen', 'products']
    for algo_scen in df.algo_scen.unique():
        dxalgos[algo_scen] = sti.ProductMix(df[df.algo_scen == algo_scen], excl_cols=excl_cols, name=f'SyphAlgo_{algo_scen}')
    return dxalgos


class DualTest(sti.HIVTest):
    """
    Dual HIV syphilis test
    """
    def __init__(self, product=None, syph_test=None, pars=None, test_prob_data=None, years=None, start=None, eligibility=None, name=None, label=None, syph_years=None, syph_prob=None, **kwargs):
        super().__init__(product=product, pars=pars, test_prob_data=test_prob_data, years=years, start=start, eligibility=eligibility, name=name, label=label, **kwargs)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed
        self.syph_test = syph_test
        self.syph_years = np.array(syph_years) if syph_years is not None else np.array([2022, 2023, 2041])
        self.syph_prob = np.array(syph_prob) if syph_prob is not None else np.array([0.0, 0.1, 0.5])

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        sim.diseases.hiv.diagnosed[pos_uids] = True
        sim.diseases.hiv.ti_diagnosed[pos_uids] = self.ti

        if self.syph_test is not None:
            prob = np.interp(self.t.now('year'), self.syph_years, self.syph_prob)
            if prob > 0:
                testers = (self.ti_tested == self.ti).uids
                if len(testers) > 0:
                    if prob >= 1.0:
                        selected = testers
                    else:
                        mask = np.random.random(len(testers)) < prob
                        selected = testers[mask]
                    if len(selected) > 0:
                        self.syph_test.ti_scheduled[selected] = self.ti

        return outcomes


def get_testing_products(add_dual=False, syph_years=None, syph_prob=None):
    """
    Define HIV products and testing interventions
    """

    scaleup_years = np.arange(1990, 2021)  # Years for testing
    years = np.arange(1990, 2041)  # Years for simulation
    n_years = len(scaleup_years)
    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_years), np.linspace(0.75, 0.85, len(years) - n_years)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_years), np.linspace(0.85, 0.95, len(years) - n_years)])
    gp_prob = np.concatenate([np.linspace(0, 0.5, n_years), np.linspace(0.5, 0.6, len(years) - n_years)])

    if add_dual:
        # Make syphilis test product
        df = sc.dataframe.read_csv('data/syph_dx.csv')
        dual_kp = sti.SyphDx(df[df.name == 'dual'], name=f'SyphDx_dual_kp')

        dual_test = sti.STITest(
            product=dual_kp,
            eligibility=lambda sim: ss.uids(),  # This is set in the DualTest class
            dt_scale=False,
            name='dual_hiv',
            label='dual_hiv',
        )
    else:
        dual_test = None

    # FSW agents who haven't been diagnosed or treated yet
    def fsw_eligibility(sim):
        return sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    fsw_testing = DualTest(
        years=years,
        test_prob_data=fsw_prob,
        name='fsw_testing',
        eligibility=fsw_eligibility,
        label='fsw_testing',
        syph_test=dual_test,
        syph_years=syph_years,
        syph_prob=syph_prob,
    )

    # Non-FSW agents who haven't been diagnosed or treated yet
    def other_eligibility(sim):
        return ~sim.networks.structuredsexual.fsw & ~sim.diseases.hiv.diagnosed & ~sim.diseases.hiv.on_art

    other_testing = sti.HIVTest(
        years=years,
        test_prob_data=gp_prob,
        name='other_testing',
        eligibility=other_eligibility,
        label='other_testing',
        # syph_test=dual_test,
    )

    # Agents whose CD4 count is below 200.
    def low_cd4_eligibility(sim):
        return (sim.diseases.hiv.cd4 < 200) & ~sim.diseases.hiv.diagnosed

    low_cd4_testing = sti.HIVTest(
        years=years,
        test_prob_data=low_cd4_prob,
        name='low_cd4_testing',
        eligibility=low_cd4_eligibility,
        label='low_cd4_testing',
        # syph_test=dual_test,
    )
    tests = [fsw_testing, other_testing, low_cd4_testing]
    if add_dual:
        tests += [dual_test]

    return tests


def make_hiv_intvs(add_dual=False, syph_years=None, syph_prob=None):

    n_art = pd.read_csv(f'data/n_art.csv').set_index('year')
    n_vmmc = pd.read_csv(f'data/n_vmmc.csv').set_index('year')
    tests = get_testing_products(add_dual=add_dual, syph_years=syph_years, syph_prob=syph_prob)
    art = sti.ART(coverage=n_art, future_coverage={'year': 2022, 'prop': 0.90})
    vmmc = sti.VMMC(coverage=n_vmmc)
    prep = sti.Prep()

    interventions = tests + [
        art,
        vmmc,
        prep,
    ]

    return interventions


class pregnancy_risk_reduction(ss.Intervention):
    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.define_pars(
            fsw_redux=ss.bernoulli(1.),  # Proportion of FSW who reduce risk during pregnancy
            rg2_redux=ss.bernoulli(1.),  # Proportion of risk group 2 who reduce risk during pregnancy
            conc_redux=ss.bernoulli(1.),  # Proportion of concurrency who reduce risk during pregnancy
        )
        self.update_pars(pars, **kwargs)
        self.define_states(
            ss.BoolState('ever_fsw'),  # Ever FSW
            ss.BoolState('ever_rg2'),  # Ever risk group 2
            ss.FloatArr('default_concurrency', default=0.0),  # Default concurrency
        )
        return

    def init_post(self):
        super().init_post()
        self.ever_fsw[:] = self.sim.networks.structuredsexual.fsw
        self.ever_rg2[:] = (self.sim.networks.structuredsexual.risk_group == 2)
        self.default_concurrency[:] = self.sim.networks.structuredsexual.concurrency
        return

    def step(self):
        # Update record of all women who have ever engaged in sex work
        self.ever_fsw[:] = self.ever_fsw[:] | self.sim.networks.structuredsexual.fsw
        self.ever_rg2[:] = self.ever_rg2[:] | (self.sim.networks.structuredsexual.risk_group == 2)
        self.default_concurrency[:] = np.maximum(self.sim.networks.structuredsexual.concurrency, self.default_concurrency)

        # Reset sexual preferences for those who are postpartum (breastfeeding and no longer pregnant)
        is_postpartum = self.sim.demographics.pregnancy.breastfeeding & ~self.sim.demographics.pregnancy.pregnant
        postpartum_fsw = self.ever_fsw & is_postpartum
        postpartum_rg2 = self.ever_rg2 & is_postpartum
        self.sim.networks.structuredsexual.fsw[postpartum_fsw] = True
        self.sim.networks.structuredsexual.risk_group[postpartum_rg2] = 2
        self.sim.networks.structuredsexual.concurrency[is_postpartum] = self.default_concurrency[is_postpartum]

        # Set temporary sexual preferences for those who are pregnant
        is_preg = self.sim.demographics.pregnancy.pregnant
        is_rg2 = self.sim.networks.structuredsexual.risk_group == 2
        is_fsw = self.sim.networks.structuredsexual.fsw
        preg_fsw = (is_preg & is_fsw).uids
        preg_rg2 = (is_preg & is_rg2).uids
        fsw_quitters = self.pars.fsw_redux.filter(preg_fsw)
        rg2_quitters = self.pars.rg2_redux.filter(preg_rg2)
        conc_reducers = self.pars.conc_redux.filter(is_preg)

        self.sim.networks.structuredsexual.fsw[fsw_quitters] = False
        self.sim.networks.structuredsexual.risk_group[rg2_quitters] = 0
        self.sim.networks.structuredsexual.concurrency[conc_reducers] = 0

        return


###########################################################################
# Algorithms
###########################################################################

__all__ = ['make_syph_testing']


def parse_scenario(scenario):
    """
    Parse a scenario name into component flags.

    Scenario names follow a 'component_component_...' convention where
    components are drawn from {'gud', 'anc', 'kp', 'plhiv'}.
    Legacy names and special aliases are mapped before parsing.

    Returns:
        use_gud (bool): GUD POC diagnostic active
        conf_channels (frozenset): confirmation test channels, subset of {'anc', 'kp', 'plhiv'}
    """
    _ALIASES = {
        'soc':        '',
        'conf':       'anc_kp_plhiv',
        'both':       'gud_anc_kp_plhiv',
        'conf_anc':   'anc',
        'conf_kp':    'kp_plhiv',   # original conf_kp routes both kp and plhiv
        'conf_fsw':   'kp',         # conf_fsw routes kp/fsw only (not plhiv)
        'conf_plhiv': 'plhiv',
    }
    normalized = _ALIASES.get(scenario, scenario)
    parts = set(normalized.split('_')) if normalized else set()
    use_gud = 'gud' in parts
    conf_channels = frozenset(parts - {'gud'})
    return use_gud, conf_channels


def make_scen_specs(use_gud, conf_channels):
    """Build scenario spec from parsed flags."""
    class ScenSpec:
        pass

    s = ScenSpec()
    s.symp_algo = 'symp_testing_' + ('gud' if use_gud else 'soc')
    # conf_algo product uses 'conf' only when ANC channel gets confirmatory test
    s.conf_algo = 'anc_confirm_' + ('conf' if 'anc' in conf_channels else 'soc')
    s.newborn_algo = 'newborn_testing_soc'
    s.newborn_test = 'exam'
    return s


def make_syph_testing(scenario='soc', rel_symp_test=1.0, rel_anc_test=1.0, plhiv_years=None, plhiv_prob=None):
    """
    Define the full testing algorithm for each scenario.
    Baseline:
        - syndromic testing: syndromic management
        - ANC testing: dual test
    Use case 1:
        - replace all syndromic management with a POC GUD diagnostic
    Use case 2:
        - POC confirmatory test after ANC dual
    """

    orig_scenario = scenario  # Save for special-case checks below

    # Resolve base scenario for flag parsing (handle special-case aliases first)
    _base = scenario
    if scenario == 'gud+secondary':
        _base = 'gud'
    elif scenario == 'soc+scaleup':
        _base = 'soc'

    use_gud, conf_channels = parse_scenario(_base)
    scenspecs = make_scen_specs(use_gud, conf_channels)

    symp_test_data = pd.read_csv('data/symp_test_prob_soc.csv')  # Risk group/sex/year-specific symptomatic testing
    years = np.array([1980, 2025, 2041])
    anc_test_data = np.array([0.1, 0.5, 0.5])*rel_anc_test  # ANC testing probabilities over time
    intv_year = 2027

    # Initialize interventions
    interventions = sc.autolist()
    dxprods = load_syph_dx()
    dxalgos = load_syph_products()

    ####################################################
    # Make GUD screening algorithms (primary-stage ulcers)
    ####################################################
    def all_ulcerative(sim):
        # Only primary-stage visible chancres + background GUD
        testers = (sim.diseases.syph.ulcerative | sim.diseases.gud.symptomatic)
        return testers.uids

    # GUD syndromic management algorithm
    symp_algo = sti.SyphTest(
        rel_test=rel_symp_test,
        product=dxalgos[scenspecs.symp_algo],
        eligibility=all_ulcerative,
        test_prob_data=symp_test_data,
        dt_scale=False,
        name='symp_algo',
        label='symp_algo',
    )

    ####################################################
    # Make secondary rash screening algorithm
    ####################################################
    def secondary_symptomatic(sim):
        # People with visible secondary syphilis rash seeking care
        return (sim.diseases.syph.rash_visible & sim.diseases.syph.secondary).uids

    # Secondary rash: syndromic management (low sensitivity ~10%)
    secondary_algo = sti.SyphTest(
        rel_test=rel_symp_test,
        product=dxprods['syndromic_rash'],  # Syndromic management for rash presentations
        eligibility=secondary_symptomatic,
        test_prob_data=symp_test_data,
        dt_scale=False,
        name='secondary_algo',
        label='secondary_algo',
    )

    interventions += [symp_algo, secondary_algo]

    ####################################################
    # Make newborn testing algorithm (scheduled by ANC when mother is positive)
    ####################################################
    newborn_algo = sti.SyphTest(
        product=dxalgos[scenspecs.newborn_algo],
        eligibility=lambda sim: ss.uids(),  # Eligibility handled via scheduling
        dt_scale=False,
        name='newborn_algo',
        label='newborn_algo',
    )

    # Individual newborn tests routed by the algorithm
    newborn_exam = sti.SyphTest(
        product=dxprods[scenspecs.newborn_test],
        eligibility=lambda sim: sim.interventions['newborn_algo'].outcomes.get('exam', ss.uids()) == sim.interventions['newborn_algo'].ti,
        dt_scale=False,
        name='newborn_exam',
        label='newborn_exam',
    )

    newborn_poc = sti.SyphTest(
        product=dxprods.get('poc_cs', dxprods[scenspecs.newborn_test]),
        eligibility=lambda sim: sim.interventions['newborn_algo'].outcomes.get('poc_cs', ss.uids()) == sim.interventions['newborn_algo'].ti,
        dt_scale=False,
        name='newborn_poc',
        label='newborn_poc',
    )

    ####################################################
    # Make ANC screening algorithms
    ####################################################
    anc_testing = sti.ANCSyphTest(
        product=dxprods['dual'],
        test_prob_data=anc_test_data,
        years=years,
        newborn_test=newborn_algo,  # Schedule newborn test when mother is ANC-positive
        name='anc_screen',
        label='anc_screen',
    )
    interventions += [anc_testing, newborn_algo, newborn_exam, newborn_poc]

    ####################################################
    # Make individual screening interventions
    ####################################################
    syndromic = sti.SyphTest(
        product=dxprods['syndromic_gud'],
        eligibility=lambda sim: sim.interventions['symp_algo'].outcomes['syndromic_gud'] == sim.interventions['symp_algo'].ti,
        dt_scale=False,
        name='syndromic',
        label='syndromic',
    )

    if orig_scenario == 'gud+secondary':
        gud_prod = dxprods['gud2']
    else:
        gud_prod = dxprods['gud']

    gud = sti.SyphTest(
        product=gud_prod,
        eligibility=lambda sim: sim.interventions['symp_algo'].outcomes['gud'] == sim.interventions['symp_algo'].ti,
        dt_scale=False,
        name='gud_test',
        label='gud_test',
    )

    # Positive results on dual test may be given a confirmatory test
    use_confirmation = bool(conf_channels)

    if orig_scenario == 'soc+scaleup':
        def to_confirm(sim):
            p1 = sim.interventions['anc_screen'].outcomes['positive'] == sim.interventions['anc_screen'].ti
            p2 = sim.interventions['dual_hiv'].outcomes['positive'] == sim.interventions['dual_hiv'].ti
            return (p1 | p2).uids
    elif use_confirmation:
        def to_confirm(sim):
            p_anc = sim.interventions['anc_screen'].outcomes['positive'] == sim.interventions['anc_screen'].ti

            if sim.t.now('year') >= intv_year:
                p_kp = sim.interventions['dual_hiv'].outcomes.get('positive', ss.uids()) == sim.interventions['dual_hiv'].ti
                p_plhiv = sim.interventions['plhiv_screen'].outcomes.get('positive', ss.uids()) == sim.interventions['plhiv_screen'].ti

                # Track which positives entered confirmation (for pathway attribution in to_treat)
                sim._conf_origin_kp = p_kp.uids if 'kp' in conf_channels else ss.uids()
                sim._conf_origin_plhiv = p_plhiv.uids if 'plhiv' in conf_channels else ss.uids()

                # KP/PLHIV channels in conf_channels bypass conf_algo — scheduled directly for confirm test
                if 'kp' in conf_channels and len(p_kp.uids):
                    sim.interventions['confirm'].ti_scheduled[p_kp.uids] = sim.ti
                if 'plhiv' in conf_channels and len(p_plhiv.uids):
                    sim.interventions['confirm'].ti_scheduled[p_plhiv.uids] = sim.ti

                # ANC always flows through conf_algo
                # (conf_algo uses conf product if 'anc' in conf_channels, else SOC product)
                return p_anc.uids
            else:
                sim._conf_origin_kp = ss.uids()
                sim._conf_origin_plhiv = ss.uids()
                return p_anc.uids
    else:
        def to_confirm(sim):
            p1 = sim.interventions['anc_screen'].outcomes['positive'] == sim.interventions['anc_screen'].ti
            return p1.uids

    conf_algo = sti.SyphTest(
        product=dxalgos[scenspecs.conf_algo],
        eligibility=to_confirm,
        dt_scale=False,
        name='conf_algo',
        label='conf_algo',
    )

    confirm = sti.SyphTest(
        product=dxprods['confirm'],
        eligibility=lambda sim: sim.interventions['conf_algo'].outcomes['confirm'] == sim.interventions['conf_algo'].ti,
        dt_scale=False,
        name='confirm',
        label='confirm',
    )

    ####################################################
    # Make PLHIV-on-ART syphilis screening
    # NB: must run BEFORE conf_algo so outcomes are available for to_confirm
    ####################################################
    plhiv_yrs = np.array(plhiv_years) if plhiv_years is not None else np.array([2019, 2020, 2030, 2041])
    plhiv_prb = np.array(plhiv_prob) if plhiv_prob is not None else np.array([0.0, 0.1, 0.5, 0.5])

    plhiv_screen = sti.SyphTest(
        product=dxprods['dual'],
        eligibility=lambda sim: (sim.diseases.hiv.on_art & sim.people.alive).uids,
        test_prob_data=plhiv_prb,
        years=plhiv_yrs,
        dt_scale=True,
        name='plhiv_screen',
        label='plhiv_screen',
    )

    testing_intvs = [
        syndromic, gud, plhiv_screen, conf_algo, confirm
    ]
    interventions += testing_intvs

    ####################################################
    # Make treatment intervention
    ####################################################

    # Get all positive adults and offer treatment
    def to_treat(sim):
        p1 = sim.interventions['conf_algo'].outcomes['treat'] == sim.interventions['conf_algo'].ti
        p2 = sim.interventions['syndromic'].outcomes['positive'] == sim.interventions['syndromic'].ti
        p3 = sim.interventions['gud_test'].outcomes['positive'] == sim.interventions['gud_test'].ti
        p4 = sim.interventions['confirm'].outcomes['positive'] == sim.interventions['confirm'].ti
        p5 = sim.diseases.syph.tertiary
        p6 = sim.interventions['secondary_algo'].outcomes['positive'] == sim.interventions['secondary_algo'].ti

        p7 = sim.interventions['dual_hiv'].outcomes.get('positive', ss.uids()) == sim.interventions['dual_hiv'].ti
        p_plhiv = sim.interventions['plhiv_screen'].outcomes.get('positive', ss.uids()) == sim.interventions['plhiv_screen'].ti

        if use_confirmation and sim.t.now('year') >= intv_year:
            # Channels NOT in conf_channels go directly to treatment
            kp_direct = p7 if 'kp' not in conf_channels else ss.uids()
            plhiv_direct = p_plhiv if 'plhiv' not in conf_channels else ss.uids()
            to_treat = (p1 | p2 | p3 | p4 | p5 | p6 | kp_direct | plhiv_direct).uids
        else:
            # Pre-intervention or non-confirmation scenarios: direct to treatment
            to_treat = (p1 | p2 | p3 | p4 | p5 | p6 | p7 | p_plhiv).uids

        # Store pathway flags for the treatment_outcomes analyzer
        # Called during treatment eligibility check, BEFORE states are cleared
        p8_direct = sim.interventions['newborn_algo'].outcomes.get('treat', ss.uids()) == sim.interventions['newborn_algo'].ti
        p8_exam = sim.interventions['newborn_exam'].outcomes.get('positive', ss.uids()) == sim.interventions['newborn_exam'].ti
        p8_poc = sim.interventions['newborn_poc'].outcomes.get('positive', ss.uids()) == sim.interventions['newborn_poc'].ti
        p8 = p8_direct | p8_exam | p8_poc

        if use_confirmation and sim.t.now('year') >= intv_year:
            # Post-intervention: attribute conf_algo/confirm outcomes back to originating pathway
            conf_treated_uids = np.asarray((p1 | p4).uids)
            kp_origin = np.asarray(getattr(sim, '_conf_origin_kp', ss.uids()))
            plhiv_origin = np.asarray(getattr(sim, '_conf_origin_plhiv', ss.uids()))
            kp_from_conf = ss.uids(np.intersect1d(conf_treated_uids, kp_origin))
            plhiv_from_conf = ss.uids(np.intersect1d(conf_treated_uids, plhiv_origin))
            anc_from_conf = ss.uids(np.setdiff1d(np.setdiff1d(conf_treated_uids, kp_origin), plhiv_origin))

            sim._tx_pathways = sc.objdict(
                gud_syndromic=(p2 | p3).uids,
                anc_screen=anc_from_conf,
                secondary_rash=p6.uids,
                kp_screen=kp_from_conf if 'kp' in conf_channels else p7.uids,
                plhiv_screen=plhiv_from_conf if 'plhiv' in conf_channels else p_plhiv.uids,
                newborn=p8.uids,
            )
        else:
            sim._tx_pathways = sc.objdict(
                gud_syndromic=(p2 | p3).uids,
                anc_screen=(p1 | p4).uids,
                secondary_rash=p6.uids,
                kp_screen=p7.uids,
                plhiv_screen=p_plhiv.uids,
                newborn=p8.uids,
            )

        # Store stage at treatment for the treatment_outcomes analyzer
        syph = sim.diseases.syph
        all_eligible = to_treat | p8.uids
        sim._tx_stages = sc.objdict()
        for stage in ['primary', 'secondary', 'early', 'late', 'tertiary']:
            sim._tx_stages[stage] = np.asarray((getattr(syph, stage) & all_eligible).uids)

        return to_treat

    treat = sti.SyphTx(
        years=years,
        eligibility=to_treat,
        name='treat',
        label='treat',
    )

    # Newborn treatment — uses NewbornTreatment which checks congenital state
    def newborn_to_treat(sim):
        p8_direct = sim.interventions['newborn_algo'].outcomes.get('treat', ss.uids()) == sim.interventions['newborn_algo'].ti
        p8_exam = sim.interventions['newborn_exam'].outcomes.get('positive', ss.uids()) == sim.interventions['newborn_exam'].ti
        p8_poc = sim.interventions['newborn_poc'].outcomes.get('positive', ss.uids()) == sim.interventions['newborn_poc'].ti
        return (p8_direct | p8_exam | p8_poc).uids

    newborn_treat = sti.NewbornTreatment(
        years=years,
        eligibility=newborn_to_treat,
        name='newborn_treat',
    )

    interventions += [treat, newborn_treat]

    # Add risk reduction,
    pregnancy_risk = pregnancy_risk_reduction()
    interventions += [pregnancy_risk]

    # Add no intervention scenario
    if orig_scenario == 'no':
        no_interventions = sc.dcp(interventions)
        for intv in no_interventions:
            intv.end = intv_year
        return no_interventions

    return interventions


def make_interventions(scenario='soc', rel_symp_test=1.0, rel_anc_test=1.0, syph_years=None, syph_prob=None, plhiv_years=None, plhiv_prob=None):
    """
    Make interventions for syphilis / HIV coinfection model
    """
    hiv_intvs = make_hiv_intvs(add_dual=True, syph_years=syph_years, syph_prob=syph_prob)
    syph_intvs = make_syph_testing(scenario=scenario, rel_anc_test=rel_anc_test, rel_symp_test=rel_symp_test, plhiv_years=plhiv_years, plhiv_prob=plhiv_prob)
    interventions = sc.autolist(hiv_intvs + syph_intvs)
    return interventions

