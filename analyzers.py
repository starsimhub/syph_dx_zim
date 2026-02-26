"""
Analyzers for the syphilis-HIV coinfection model
"""

# %% Imports and settings
import numpy as np
import sciris as sc
import starsim as ss
import pandas as pd
import networkx as nx

import stisim as sti
import pylab as pl

def count(arr): return np.count_nonzero(arr)


class syph_dalys(ss.Analyzer):
    def __init__(self, start=None, life_expectancy=59, r=0.03, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.start = start
        self.disability_weights = sc.objdict(
            primary=0.006,  # GBD 2019
            tertiary=0.547,  # GBD 2019
            congenital=0.315,  # https://hqlo.biomedcentral.com/articles/10.1186/s12955-024-02234-1
        )
        self.life_expectancy = life_expectancy  # Should typically use country-specific values
        self.dur_tertiary = 5  # Average duration of tertiary syphilis (?)
        self.r = r  # Discount rate, used for YLL
        self.dr = 1/(1+self.r)**np.arange(100)  # Precalculate discount rates

        return

    def init_results(self):
        super().init_results()
        results = [
            ss.Result('yll', dtype=int),
            ss.Result('yld', dtype=int),
            ss.Result('dalys', dtype=int),
            ss.Result('cum_dalys', dtype=int),
        ]
        self.define_results(*results)
        return

    def init_pre(self, sim):
        super().init_pre(sim)
        if self.start is None: self.start = self.sim.t.yearvec[0]

    def get_yld(self, sim):
        """ Method of calculating YLD depends on whether we're using incidence or prevalence / hybrid """
        pass

    def discounted_sum(self, weight, duration_arr):
        """ Calculate NPV for a series of durations """
        duration_arr = duration_arr.astype(int)
        return sum([weight*sum(self.dr[:i]) for i in duration_arr])

    def step(self):
        sim = self.sim
        syph = sim.diseases.syph
        ppl = sim.people
        ti = self.ti

        after_start = (self.t.now('year') >= self.start)

        if after_start:

            # Count deaths. We don't include adult deaths, which need life tables, and there are hardly any of them
            # If you do want them, add | (syph.ti_dead == sim.ti) to the line below
            new_deaths = (syph.ti_nnd == ti)
            if new_deaths.any():
                lex = np.maximum(self.life_expectancy - ppl.age[new_deaths], 0)
                self.results.yll[ti] = self.discounted_sum(1, lex)
            self.results.yld[ti] = self.get_yld(sim)
            self.results.dalys[ti] = self.results.yld[ti] + self.results.yll[ti]

            if np.isnan(self.results.dalys[ti]):
                raise ValueError

        return


class syph_idalys(syph_dalys):
    """
    Syphilis DALYs, calculated as incidence DALYs individuated by sequelae and death.
    Reference: https://pophealthmetrics.biomedcentral.com/articles/10.1186/1478-7954-10-19
    """
    def get_yld(self, sim):
        syph = sim.diseases.syph
        dw = self.disability_weights
        ppl = sim.people

        # Primary
        new_primary = (syph.ti_primary == self.ti)
        if new_primary.any():
            dur_primary = (syph.ti_secondary[new_primary] - self.ti)
            yld_primary = sum(dur_primary*dw.primary)  # Not discounting because duration so short
        else:
            yld_primary = 0

        # YLD tertiary
        new_tertiary = (syph.ti_tertiary == self.ti)
        if new_tertiary.any():
            dur_tertiary = np.maximum(self.life_expectancy - ppl.age[new_tertiary], self.dur_tertiary)
            yld_tertiary = self.discounted_sum(dw.tertiary, dur_tertiary)  # sum(dur_tertiary*dw.tertiary)
        else:
            yld_tertiary = 0

        # YLD congenital
        new_congenital = (syph.ti_congenital == sim.ti)
        if new_congenital.any():
            lex = np.maximum(self.life_expectancy - ppl.age[new_congenital], 0)
            yld_congenital = self.discounted_sum(dw.tertiary, lex)

        else:
            yld_congenital = 0

        return yld_primary + yld_tertiary + yld_congenital


class syph_hdalys(syph_dalys):
    """
    Syphilis DALYs, calculated as hybrid DALYs
    Reference: https://pophealthmetrics.biomedcentral.com/articles/10.1186/1478-7954-10-19
        - YLL are accrued at the point of death
        - YLD are accrued annually
    """
    def get_yld(self, sim):
        syph = sim.diseases.syph
        dw = self.disability_weights
        yld = 0
        for state, weight in dw.items():
            yld += sum(weight * syph[state])
        return yld


class epi_ts(ss.Analyzer):
    """ Save epi results for plotting """
    def __init__(self):
        super().__init__()
        return

    def init_results(self):
        super().init_results()
        results = [
            ss.Result('syph.new_infections_f', dtype=float),
            ss.Result('syph.new_infections_m', dtype=float),
            ss.Result('hiv.new_infections_f', dtype=float),
            ss.Result('hiv.new_infections_m', dtype=float),
        ]
        self.define_results(*results)
        return

    def step(self):
        pass

    def finalize(self):
        super().finalize()
        sim = self.sim
        syph_res = sim.results.syph
        hiv_res = sim.results.hiv

        # Syphilis stats
        self.results['syph.new_infections_f'] = syph_res.new_infections_f
        self.results['syph.new_infections_m'] = syph_res.new_infections_m
        self.results['hiv.new_infections_f'] = hiv_res.new_infections_f
        self.results['hiv.new_infections_m'] = hiv_res.new_infections_m

        return


class transmission_by_stage(ss.Analyzer):
    """ Track which disease stage (primary/secondary/early/late/tertiary) each transmission comes from """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'transmission_by_stage'
        self.stages = ['primary', 'secondary', 'early', 'late', 'tertiary']
        self.transmission_modes = ['sex', 'mtc']

    def init_results(self):
        super().init_results()
        results = sc.autolist()
        for tm in self.transmission_modes:
            for stage in self.stages:
                results += ss.Result(f'new_{tm}_{stage}', dtype=int, scale=True)
        self.define_results(*results)

    def step(self):
        sim = self.sim
        ti = self.ti
        syph = sim.diseases.syph
        new_trans = dict(
            sex=syph.ti_transmitted_sex == ti,
            mtc=syph.ti_transmitted_mtc == ti,
        )
        for tm in self.transmission_modes:
            for stage in self.stages:
                self.results[f'new_{tm}_{stage}'][ti] += count(new_trans[tm] & getattr(syph, stage))


class treatment_outcomes(ss.Analyzer):
    """
    Track treatment outcomes by diagnostic pathway, sex, and HIV status.
    Requires that to_treat() in interventions.py stores sim._tx_pathways.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'treatment_outcomes'
        self.pathways = ['gud_syndromic', 'anc_screen', 'secondary_rash', 'kp_screen']
        self.test_intvs = {
            'gud_syndromic': ['syndromic', 'gud_test'],
            'anc_screen': ['anc_screen'],
            'secondary_rash': ['secondary_algo'],
            'kp_screen': ['dual_hiv'],
        }

    def init_results(self):
        super().init_results()
        results = sc.autolist()

        # Per pathway: treated, success, unnecessary (overtreated), failure
        for pw in self.pathways:
            for oc in ['treated', 'success', 'unnecessary', 'failure']:
                results += ss.Result(f'{pw}_{oc}', dtype=int, scale=True)
                results += ss.Result(f'{pw}_{oc}_f', dtype=int, scale=True)
                results += ss.Result(f'{pw}_{oc}_m', dtype=int, scale=True)
                results += ss.Result(f'{pw}_{oc}_hivpos', dtype=int, scale=True)
                results += ss.Result(f'{pw}_{oc}_hivneg', dtype=int, scale=True)

        # Per pathway: false negatives (tested but missed active infection)
        for pw in self.pathways:
            results += ss.Result(f'{pw}_missed', dtype=int, scale=True)
            results += ss.Result(f'{pw}_missed_f', dtype=int, scale=True)
            results += ss.Result(f'{pw}_missed_m', dtype=int, scale=True)
            results += ss.Result(f'{pw}_missed_hivpos', dtype=int, scale=True)
            results += ss.Result(f'{pw}_missed_hivneg', dtype=int, scale=True)

        # Overall active prevalence (for context)
        results += ss.Result('n_active', dtype=int, scale=True)

        self.define_results(*results)

    def step(self):
        sim = self.sim
        ti = self.ti
        syph = sim.diseases.syph
        hiv = sim.diseases.hiv
        ppl = sim.people

        # Get treatment outcomes (set during treat.step(), before analyzer runs)
        treat = sim.interventions['treat']
        tx_oc = treat.outcomes.get('syph', sc.objdict())
        successful = np.asarray(tx_oc.get('successful', ss.uids()))
        unsuccessful = np.asarray(tx_oc.get('unsuccessful', ss.uids()))
        unnecessary = np.asarray(tx_oc.get('unnecessary', ss.uids()))
        all_treated = np.concatenate([successful, unsuccessful, unnecessary])

        # Get pathway flags (stored by to_treat before states cleared)
        pw_flags = getattr(sim, '_tx_pathways', None)
        if pw_flags is None:
            return

        # Attribute treatments to pathways
        for pw in self.pathways:
            pw_uids = np.asarray(pw_flags.get(pw, ss.uids()))
            if len(pw_uids) == 0 and len(all_treated) == 0:
                continue

            pw_treated = np.intersect1d(all_treated, pw_uids)
            pw_success = np.intersect1d(successful, pw_uids)
            pw_unnecessary = np.intersect1d(unnecessary, pw_uids)
            pw_failure = np.intersect1d(unsuccessful, pw_uids)

            self._record(f'{pw}_treated', pw_treated, ti, ppl, hiv)
            self._record(f'{pw}_success', pw_success, ti, ppl, hiv)
            self._record(f'{pw}_unnecessary', pw_unnecessary, ti, ppl, hiv)
            self._record(f'{pw}_failure', pw_failure, ti, ppl, hiv)

        # False negatives: tested negative but had active syphilis
        # Since these people weren't treated, their active state persists
        for pw, test_names in self.test_intvs.items():
            missed = ss.uids()
            for tn in test_names:
                intv = sim.interventions[tn]
                tested_neg = (intv.ti_negative == ti)
                missed = missed | (tested_neg & syph.active).uids
            self._record(f'{pw}_missed', np.asarray(missed), ti, ppl, hiv)

        # Overall active prevalence
        self.results['n_active'][ti] = count(syph.active)

    def _record(self, prefix, uids, ti, ppl, hiv):
        """Record a count disaggregated by sex and HIV status"""
        self.results[prefix][ti] = len(uids)
        if len(uids):
            self.results[f'{prefix}_f'][ti] = count(ppl.female[uids])
            self.results[f'{prefix}_m'][ti] = count(ppl.male[uids])
            self.results[f'{prefix}_hivpos'][ti] = count(hiv.infected[uids])
            self.results[f'{prefix}_hivneg'][ti] = count(~hiv.infected[uids])


class NetworkSnapshot(ss.Analyzer):
    """
    Capture a snapshot of network properties at a specified year.
    Used for the supplementary network structure figure.
    """
    def __init__(self, year=2020, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.year = year
        self.name = 'network_snapshot'
        self.risk_group_data = None
        self.debut_data = None
        self.lifetime_partners_data = None
        self.partnership_by_age = None
        self.rel_dur_data = None

    def step(self):
        if self.sim.t.yearvec[self.ti] == self.year:
            self._capture_snapshot()

    def _capture_snapshot(self):
        sim = self.sim
        nw = sim.networks.structuredsexual
        ppl = sim.people
        active = nw.participant & ppl.alive

        # Panel C: Risk group composition
        rg_data = {}
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            for rg in [0, 1, 2]:
                rg_data[(sex_label, rg)] = int(((nw.risk_group == rg) & sex_bool & active).count())
            rg_data[(sex_label, 'total')] = int((sex_bool & active).count())
        rg_data[('Female', 'fsw')] = int((nw.fsw & ppl.female & active).count())
        rg_data[('Male', 'client')] = int((nw.client & ppl.male & active).count())
        self.risk_group_data = rg_data

        # Panel A: Lifetime partners for debuted agents only
        debuted = nw.participant & ppl.alive & (ppl.age >= nw.debut)
        lp_data = {}
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            mask = sex_bool & debuted
            lp_data[sex_label] = np.array(nw.lifetime_partners[mask])
        self.lifetime_partners_data = lp_data

        # Panel D: Sexual debut age by sex
        debut_data = {}
        for sex_label, sex_bool in [('Female', ppl.female), ('Male', ppl.male)]:
            mask = sex_bool & debuted
            debut_data[sex_label] = np.array(nw.debut[mask])
        self.debut_data = debut_data

        # Panel E: Female partnership status by age
        age_bins = np.arange(15, 51)
        pba = dict(age_bins=age_bins, prop_stable=[], prop_casual=[])
        for age in age_bins:
            in_age = ppl.female & ppl.alive & (ppl.age >= age) & (ppl.age < age + 1)
            n_total = int(in_age.count())
            if n_total > 0:
                n_stable = int((in_age & (nw.stable_partners >= 1)).count())
                n_casual = int((in_age & (nw.casual_partners >= 1)).count())
                pba['prop_stable'].append(n_stable / n_total)
                pba['prop_casual'].append(n_casual / n_total)
            else:
                pba['prop_stable'].append(np.nan)
                pba['prop_casual'].append(np.nan)
        pba['prop_stable'] = np.array(pba['prop_stable'])
        pba['prop_casual'] = np.array(pba['prop_casual'])
        self.partnership_by_age = pba

    def finalize(self):
        super().finalize()
        # Panel E: Relationship durations by edge type
        nw = self.sim.networks.structuredsexual
        dur_by_type = {0: [], 1: []}  # stable=0, casual=1
        dt_year = self.sim.t.dt_year
        for _, rels in nw.relationship_durs.items():
            for rel in rels:
                etype = rel.get('edge_type', -1)
                if etype in dur_by_type:
                    dur_by_type[etype].append(rel['dur'] * dt_year)
        self.rel_dur_data = dur_by_type


def make_analyzers(which='all', extra_analyzers=None):
    analyzers = sc.autolist()
    if which in ['all', 'stis']:
        analyzers += [
            sti.coinfection_stats('syph', 'hiv', name='coinfection_stats'),
            sti.coinfection_stats('syph', 'hiv', disease1_infected_state_name='active', name='active_coinfection_stats'),
            epi_ts(),
            sti.sw_stats(diseases=['syph', 'hiv']),
            syph_idalys(),
            transmission_by_stage(),
            treatment_outcomes(),
        ]
    analyzers += sc.autolist(extra_analyzers)
    return analyzers
