import matplotlib.pyplot as plt

import importlib

import auxplt

importlib.reload(auxplt)

def gamphi(order, phi, gam, wvl, col):

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    for i in range(len(wvl)):

        ax.plot(phi, gam[i, :], color = col[i])

    ax.set_xlim(0, 90)
    ax.set_ylim(0, 45)

    ax.set_xlabel(r'$\phi$, [deg]', fontsize = 30)
    ax.set_ylabel(r'$\gamma$, [deg]', fontsize = 30)

    ax.tick_params(axis = 'both', labelsize = 20)

    auxplt.savepdf('rainbow/' + order + '/phigam')

def maxgam(order, phi, gam, gam_grid, wvl, col):

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    for i in range(len(wvl)):

        ax.scatter(phi, gam[i, :], color = col[i], s = 7)

    for j in range(len(gam_grid)):

        ax.axhline(y = gam_grid[j], linewidth = 0.5)

    ax.set_xlim(50, 67.75)
    ax.set_ylim(40.5, 42.25)

    ax.set_xlabel(r'$\phi$, [deg]', fontsize = 30)
    ax.set_ylabel(r'$\gamma$, [deg]', fontsize = 30)

    ax.tick_params(axis = 'both', labelsize = 20)

    auxplt.savepdf('rainbow/' + order + '/maxgam')

def hist(order, gam, gam_grid, wvl, col):

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    for i in range(len(wvl)):

        ax.hist(gam[i, :], gam_grid, histtype = 'bar', facecolor = col[i], alpha = 0.5)

    ax.set_xlabel(r'$\gamma$, [deg]', fontsize = 30)
    ax.set_ylabel('Frequency', fontsize = 30)

    ax.tick_params(axis = 'both', labelsize = 20)

    ax.set_xlim(40, 42.5)

    auxplt.savepdf('rainbow/' + order + '/hist')
