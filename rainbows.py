import numpy as np
import matplotlib.pyplot as plt

import importlib
import sys

import first
import auxplt
import plot
import rainbow

importlib.reload(rainbow)
importlib.reload(first)
importlib.reload(auxplt)
importlib.reload(plot)

wvl = np.array([425.0, 475.0, 550.0, 650.0])

col = ['m', 'b', 'g', 'r']

phif = np.linspace(0, 90, 100000)
phic = np.linspace(0, 90, 1000)

gam1f = rainbow.gamma(1, np.deg2rad(phif), wvl)
gam1c = rainbow.gamma(1, np.deg2rad(phic), wvl)

gam_grid_1, g1, I1uw, I1w = first.order(phif, gam1f, wvl, col)

plot.gamphi('1', phif, gam1f, wvl, col)
plot.maxgam('1', phic, gam1c, gam_grid_1, wvl, col)
plot.hist('1', gam1f, gam_grid_1, wvl, col)

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

for i in range(len(wvl)):

    ax.plot(g1, I1w[i, :], color = col[i])

    ax.set_xlim(0, 45)

    ax.set_xlabel(r'$\gamma$, [deg]', fontsize = 30)
    ax.set_ylabel('Intensity, [%]', fontsize = 30)

    ax.tick_params(axis = 'both', labelsize = 20)

    auxplt.savepdf('rainbow/1/inten_w')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

for i in range(len(wvl)):

    ax.plot(g1, I1uw[i, :], color = col[i])

    ax.set_xlim(0, 45)

    ax.set_xlabel(r'$\gamma$, [deg]', fontsize = 30)
    ax.set_ylabel('Intensity, [%]', fontsize = 30)

    ax.tick_params(axis = 'both', labelsize = 20)

    auxplt.savepdf('rainbow/1/inten_uw')
