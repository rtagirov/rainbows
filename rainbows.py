import numpy as np
import matplotlib.pyplot as plt

import importlib
import sys

import zero
import first
import second
import auxplt
import rainbow

importlib.reload(rainbow)
importlib.reload(zero)
importlib.reload(first)
importlib.reload(second)
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

# coarse is different for zero, first and second orders
phi0c = np.linspace(0, 90, 10000)
phi1c = np.linspace(0, 90, 1000)
phi2c = np.linspace(0, 90, 1000)

# angle of deflection gamma for the zero/first/second order rainbow
# corresponding to fine and coarse phi grids
# function 'gamma' is described in the rainbow.py module
gam0f = rainbow.gamma(0, np.deg2rad(phif), wvl)
gam0c = rainbow.gamma(0, np.deg2rad(phi0c), wvl)

gam1f = rainbow.gamma(1, np.deg2rad(phif), wvl)
gam1c = rainbow.gamma(1, np.deg2rad(phi1c), wvl)

gam2f = rainbow.gamma(2, np.deg2rad(phif), wvl)
gam2c = rainbow.gamma(2, np.deg2rad(phi2c), wvl)

# plotting gamma(phi) dependencies for rainbows of all orders
# -----------------------------------------------

# initiate figure (see auxplt.py)
fig, ax = auxplt.initfig(1, 1, (12.0, 10.0))

# loop over wavelengths/colors
for i in range(len(wvl)):

    ax.plot(phif, gam0f[i, :], color = col[i])
    ax.plot(phif, gam1f[i, :], color = col[i])
    ax.plot(phif, gam2f[i, :], color = col[i])

    ax.set_xlim(0, 90)
    ax.set_ylim(0, 90)

    ax.set_xlabel(r'$\phi$, [deg]', fontsize = 30)
    ax.set_ylabel(r'$\gamma$, [deg]', fontsize = 30)

#   set the size of tick labels along both axes
    ax.tick_params(axis = 'both', labelsize = 20)

    ax.text(42, 20, 'zero order', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 20)
    ax.text(18, 30, 'first order', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 20)
    ax.text(47, 80, 'second order', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 20)

#   invoke savepdf function from auxplt.py module (auxiliary plotting functions)
    auxplt.savepdf('gamphi')
# -----------------------------------------------

# gam_grid_0 --- gamma binning for the zero order rainbow
# gam_grid_1 --- gamma binning for the first order rainbow
# gam_grid_2 --- gamma binning for the second order rainbow
# g0 --- values of gamma corresponding to the gam_grid_0 binning (calculated as the middle value of the bin)
# g1 --- values of gamma corresponding to the gam_grid_1 binning (calculated as the middle value of the bin)
# g2 --- values of gamma corresponding to the gam_grid_2 binning (calculated as the middle value of the bin)
# I0 --- relative (to incindent) intensity in percents for the zero order rainbow
# I1 --- relative (to incindent) intensity in percents for the first order rainbow
# I2 --- relative (to incindent) intensity in percents for the second order rainbow
# function 'order' is described in zero.py/first.py/second.py modules
gam_grid_0, g0, I0 = zero.order(phif, gam0f, wvl)
gam_grid_1, g1, I1 = first.order(phif, gam1f, wvl)
gam_grid_2, g2, I2 = second.order(phif, gam2f, wvl)

# plotting intensities for rainbows of all orders
# -----------------------------------------------

# initiate figure (see auxplt.py)
fig, ax = auxplt.initfig(1, 1, (12.0, 10.0))

# loop over wavelengths/colors
for i in range(len(wvl)):

    ax.plot(g0, I0[i, :], color = col[i])
    ax.plot(g1, I1[i, :], color = col[i])
    ax.plot(g2, I2[i, :], color = col[i])

    ax.set_xlim(0, 90)

    ax.set_xlabel(r'$\gamma$, [deg]', fontsize = 30)
    ax.set_ylabel('Intensity, [%]', fontsize = 30)

#   set the size of tick labels along both axes
    ax.tick_params(axis = 'both', labelsize = 20)

    ax.text(14, 0.23, 'zero order', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 20)
    ax.text(43, 0.13, 'first order', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 20)
    ax.text(55, 0.05, 'second order', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 20)

    ax.text(45.5, 0.08, "Alexander's band", bbox = dict(facecolor = 'red', alpha = 0.25), fontsize = 20, rotation = 90)

#   invoke savepdf function from auxplt.py module (auxiliary plotting functions)
    auxplt.savepdf('inten')
# -----------------------------------------------

# invoking the plotting functions for zero/first/second order
# rainbows from zero.py/first.py/second.py modules
# gamphi is the dependence of deflection angle on the angle incidence
# maxgam is the same plot zoomed in around the maximal value of gamma
# hist is the color histogram corresponding to the gamma(phi) dependence
zero.plot_gamphi(phif, gam0f, wvl, col)
zero.plot_maxgam(phi0c, gam0c, gam_grid_0, wvl, col)

first.plot_gamphi(phif, gam1f, wvl, col)
first.plot_maxgam(phi1c, gam1c, gam_grid_1, wvl, col)

second.plot_gamphi(phif, gam2f, wvl, col)
second.plot_maxgam(phi2c, gam2c, gam_grid_2, wvl, col)

#zero.plot_hist(gam0f, gam_grid_0, wvl, col)
#first.plot_hist(gam1f, gam_grid_1, wvl, col)
second.plot_hist(gam2f, gam_grid_2, wvl, col)
