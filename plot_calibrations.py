"""
Plot calibrations for the STI model
"""

if __name__ == '__main__':

    which = 'all'  # 'hiv' or 'all'
    from plot_sims import plot_calibrations
    plot_calibrations(dislist=which, start_year=1990)

    print('Done!')


