"""
Utils and defaults
"""
import sciris as sc
import numpy as np


def set_font(size=None, font='Libertinus Sans'):
    sc.fonts(add=sc.thisdir(aspath=True) / 'assets' / 'LibertinusSans-Regular.otf')
    sc.options(font=font, fontsize=size)
    return


percentile_pairs = [[.01, .99], [.1, .9], [.25, .75]]  # Order by wide to narrow (for alpha shading in plots)
percentiles = [.5] + [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]


def check_sim_alive(sim):
    """Check that syphilis and HIV are both sustained in the simulation.

    Used by both calibration (as check_fn) and run_msim (post-run filter).
    Syphilis is checked over the full time series to avoid false failures
    from borderline parsets that briefly dip to zero near the end.
    HIV is checked over the last 60 timesteps (~5 years) so we require it
    to be alive in the contemporary period, not just historically.
    """
    if sim is None:
        return False
    if np.sum(sim.results.syph.new_infections) == 0:
        return False
    if np.median(sim.results.hiv.prevalence_15_49[-60:]) < 0.05:
        return False
    return True


def count(arr): return np.count_nonzero(arr)


def get_y(df, which, rname):
    if which == 'single': y = df[rname]
    elif which == 'multi': y = df[(rname, '50%')]
    return y


def plot_single(ax, data, model, data_var_name, var_name, annualize=False, alpha=1, smooth=False):
    ax.scatter(data.time, data[data_var_name], label='Data', color='k')
    if annualize:
        x = np.unique(model['year'])
        y = model.groupby(by='year')[var_name].sum()
    else:
        x = model['timevec']
        y = model[var_name]
    if smooth:
        y = y.rolling(10, min_periods=1).mean()
    ax.plot(x, y, label='Model', alpha=alpha)
    sc.SIticks(ax)
    ax.set_ylim(bottom=0)
    return ax

