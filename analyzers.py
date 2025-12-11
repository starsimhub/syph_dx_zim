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


def make_analyzers(which='all', extra_analyzers=None):
    analyzers = sc.autolist()
    if which in ['all', 'stis']:
        analyzers += [
            sti.coinfection_stats('syph', 'hiv'),
            sti.sw_stats(diseases=['syph', 'hiv']),
            syph_idalys(),
        ]
    analyzers += sc.autolist(extra_analyzers)
    return analyzers
