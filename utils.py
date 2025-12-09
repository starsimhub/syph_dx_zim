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
percentiles = [percentile for percentile_pair in percentile_pairs for percentile in percentile_pair]


def count(arr): return np.count_nonzero(arr)


def get_y(df, which, rname):
    if which == 'single': y = df[rname]
    elif which == 'multi': y = df[(rname, '50%')]
    return y



def plot_single(ax, data, model, data_var_name, var_name, annualize=False, alpha=1, smooth=False):
    ax.scatter(data.year, data[data_var_name], label='Data', color='k')
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

