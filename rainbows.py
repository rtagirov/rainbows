import numpy as np
import matplotlib.pyplot as plt

import importlib
import sys

import first
import auxplt
import rainbow

importlib.reload(rainbow)
importlib.reload(first)
importlib.reload(auxplt)

# rainbow wavelengths accoring to main_article.pdf
wvl = np.array([425.0, 475.0, 550.0, 650.0])

# corresponding colors for the plots
col = ['m', 'b', 'g', 'r']

# angle of incidence phi
# 'f' stands for fine
# 'c' stands for coarse
# fine is for sufficient precision in the calculations
# coarse is for plotting the maxgam figure
phif = np.linspace(0, 90, 100000) # as prescribed in main_article.pdf, see caption for Fig. 5
phic = np.linspace(0, 90, 1000)

# angle of deflection gamma for the first order rainbow
# corresponding to fine and coarse phi grids
# function 'gamma' is described in the rainbow.py module
gam1f = rainbow.gamma(1, np.deg2rad(phif), wvl)
gam1c = rainbow.gamma(1, np.deg2rad(phic), wvl)

# gam_grid_1 --- gamma binning for the first order rainbow
# g1 --- values of gamma corresponding to the gamma binning (calculated as the middle value of the bin)
# I1 --- relative (to incindent) intensity in percents for the first order rainbow
# function 'order' is described in the first.py module
gam_grid_1, g1, I1 = first.order(phif, gam1f, wvl)

# invoking the plotting functions for the first order rainbow from the first.py module
# gamphi is the dependence of deflection angle on the angle incidence
# maxgam is the same plot zoomed in around the maximal value of gamma
# hist is the color histogram corresponding to the gamma(phi) dependence
first.plot_gamphi(phif, gam1f, wvl, col)
first.plot_maxgam(phic, gam1c, gam_grid_1, wvl, col)
first.plot_hist(gam1f, gam_grid_1, wvl, col)

# initiate figure (see auxplt.py)
fig, ax = auxplt.initfig(1, 1, (12.0, 10.0))

# loop over wavelengths/colors
for i in range(len(wvl)):

    ax.plot(g1, I1[i, :], color = col[i])

    ax.set_xlim(0, 45)

    ax.set_xlabel(r'$\gamma$, [deg]', fontsize = 30)
    ax.set_ylabel('Intensity, [%]', fontsize = 30)

#   set the size of tick labels along both axes
    ax.tick_params(axis = 'both', labelsize = 20)

#   invoke savepdf function from auxplt.py module (auxiliary plotting functions)
    auxplt.savepdf('inten')
