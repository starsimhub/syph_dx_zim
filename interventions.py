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
    def __init__(self, product=None, syph_test=None, pars=None, test_prob_data=None, years=None, start=None, eligibility=None, name=None, label=None, **kwargs):
        super().__init__(product=product, pars=pars, test_prob_data=test_prob_data, years=years, start=start, eligibility=eligibility, name=name, label=label, **kwargs)
        if self.eligibility is None:
            self.eligibility = lambda sim: ~sim.diseases.hiv.diagnosed
        self.syph_test = syph_test
        self.syph_start = start if start is not None else 2028

    def step(self, uids=None):
        sim = self.sim
        outcomes = super().step(uids=uids)
        pos_uids = outcomes['positive']
        sim.diseases.hiv.diagnosed[pos_uids] = True
        sim.diseases.hiv.ti_diagnosed[pos_uids] = self.ti

        if self.t.now('year') >= self.syph_start:
            testers = (self.ti_tested == self.ti).uids
            if len(testers) > 0:
                self.syph_test.ti_scheduled[testers] = self.ti

        return outcomes


def get_testing_products():
    """
    Define HIV products and testing interventions
    """

    scaleup_years = np.arange(1990, 2021)  # Years for testing
    years = np.arange(1990, 2041)  # Years for simulation
    n_years = len(scaleup_years)
    fsw_prob = np.concatenate([np.linspace(0, 0.75, n_years), np.linspace(0.75, 0.85, len(years) - n_years)])
    low_cd4_prob = np.concatenate([np.linspace(0, 0.85, n_years), np.linspace(0.85, 0.95, len(years) - n_years)])
    gp_prob = np.concatenate([np.linspace(0, 0.5, n_years), np.linspace(0.5, 0.6, len(years) - n_years)])

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

    return fsw_testing, other_testing, low_cd4_testing, dual_test


def make_hiv_intvs():

    n_art = pd.read_csv(f'data/n_art.csv').set_index('year')
    n_vmmc = pd.read_csv(f'data/n_vmmc.csv').set_index('year')
    fsw_testing, other_testing, low_cd4_testing, dual_test = get_testing_products()
    art = sti.ART(coverage_data=n_art)
    vmmc = sti.VMMC(coverage_data=n_vmmc)
    prep = sti.Prep()

    interventions = [
        fsw_testing,
        other_testing,
        low_cd4_testing,
        dual_test,
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

        # Reset sexual preferences for those who are postpartum
        is_postpartum = self.sim.demographics.pregnancy.postpartum & ~self.sim.demographics.pregnancy.pregnant
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


def make_scen_specs(scenario):
    # Define product set for each scenario
    class ScenSpec:
        def __init__(self, name=None, symp_algo=None, conf_algo=None, newborn_algo=None, newborn_test=None):
            self.name = name
            self.newborn_test = newborn_test

            # Define algorithms
            self.symp_algo = 'symp_testing_'+symp_algo
            self.conf_algo = 'anc_confirm_'+conf_algo
            self.newborn_algo = 'newborn_testing_'+newborn_algo

    scenlist = [
        ScenSpec(name='soc',    symp_algo='soc',  conf_algo='soc',  newborn_algo='soc', newborn_test='exam'),
        ScenSpec(name='gud',    symp_algo='gud',  conf_algo='soc',  newborn_algo='soc', newborn_test='exam'),
        ScenSpec(name='conf',   symp_algo='soc',  conf_algo='conf', newborn_algo='soc', newborn_test='exam'),
        ScenSpec(name='cs',     symp_algo='soc',  conf_algo='soc',  newborn_algo='cs',  newborn_test='dx_cs'),
        ScenSpec(name='no',     symp_algo='soc',  conf_algo='soc',  newborn_algo='soc', newborn_test='exam'),
    ]

    scendict = ss.ndict(scenlist)
    scenspecs = scendict[scenario]

    return scenspecs


def make_syph_testing(scenario='soc'):
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

    orig_scenario = scenario[:]  # Save the original scenario for later use
    if scenario == 'gud+secondary':  # Special case for gud+secondary
        scenario = 'gud'
    if scenario == 'soc+scaleup':  # Special case for soc+scaleup
        scenario = 'soc'

    symp_test_data = pd.read_csv('data/symp_test_prob_soc.csv')  # Risk group/sex/year-specific symptomatic testing
    years = np.array([1990, 2025, 2041])
    anc_test_data = np.array([0.1, 0.5, 0.9])

    # Make scenario specifications - these determine which algorithm gets created and which tests get added
    scenspecs = make_scen_specs(scenario)
    intv_year = 2027

    # Initialize interventions
    interventions = sc.autolist()
    dxprods = load_syph_dx()
    dxalgos = load_syph_products()

    ####################################################
    # Make symptomatic screening algorithms
    ####################################################
    def all_symptomatic(sim):
        testers = (sim.diseases.syph.primary | sim.diseases.gud.symptomatic)
        return testers.uids

    # Determine how symptomatic people are managed
    symp_algo = sti.SyphTest(
        product=dxalgos[scenspecs.symp_algo],  # SOC prior to 2027.
        eligibility=all_symptomatic,
        test_prob_data=symp_test_data,
        dt_scale=False,
        name='symp_algo',
        label='symp_algo',
    )

    interventions += symp_algo

    ####################################################
    # Make ANC screening algorithms
    ####################################################
    anc_testing = sti.ANCSyphTest(
        product=dxprods['dual'],
        test_prob_data=anc_test_data,
        years=years,
        name='anc_screen',
        label='anc_screen',
    )
    interventions += anc_testing

    ####################################################
    # Make individual screening interventions
    ####################################################
    syndromic = sti.SyphTest(
        product=dxprods['syndromic'],
        eligibility=lambda sim: sim.interventions['symp_algo'].outcomes['syndromic'] == sim.interventions['symp_algo'].ti,
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
    if orig_scenario == 'soc+scaleup':
        def to_confirm(sim):
            p1 = sim.interventions['anc_screen'].outcomes['positive'] == sim.interventions['anc_screen'].ti
            p2 = sim.interventions['dual_hiv'].outcomes['positive'] == sim.interventions['dual_hiv'].ti
            return (p1 | p2).uids
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

    testing_intvs = [
        syndromic, gud, conf_algo, confirm
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
        # to_treat = ( p2 | p3 ).uids
        to_treat = (p1 | p2 | p3 | p4 | p5).uids
        return to_treat

    treat = sti.SyphTx(
        years=years,
        eligibility=to_treat,
        # fetus_treat_eff=0.1,
        # treat_eff_reduced=0.05,
        name='treat',
        label='treat',
    )
    interventions += [treat]

    # Add risk reduction,
    pregnancy_risk = pregnancy_risk_reduction()
    interventions += [pregnancy_risk]

    # Add no intervention scenario
    if scenario == 'no':
        no_interventions = sc.dcp(interventions)
        for intv in no_interventions:
            intv.end = intv_year
        return no_interventions

    return interventions


def make_interventions(which='all', scenario='soc'):
    """
    Make interventions for syphilis / HIV coinfection model
    """
    interventions = sc.autolist()

    # HIV interventions
    if which in ['all', 'hiv']:
        hiv_intvs = make_hiv_intvs()
        interventions += hiv_intvs

    # Syphilis testing interventions
    elif which in ['all', 'stis']:
        syph_intvs = make_syph_testing(scenario=scenario)
        interventions += syph_intvs

    else:
        raise NotImplementedError(f'Intervention set "{which}" not recognized')

    return interventions

