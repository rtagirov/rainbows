import numpy as np
import matplotlib.pyplot as plt
import sys
import importlib
import rainbow
import auxplt

importlib.reload(auxplt)
importlib.reload(rainbow)

# function 'order' in this module calculates an
# optimal gamma binning for the zero order rainbow;
# this calculation is the first part of the intensity calculation;
# at the end of the function the binning is passed onto
# the function 'intensity' declared in the rainbow.py module;
# the optimal binning is necessary because if the spacing between
# gamma points is equal over all colors then for some colors the calculation
# will have less precision because the maximum gamma for some colors will be close to the
# center of the bin, thus resulting in the intensity underestimation;
# this function constructs the gamma grid in such a way that the maximum
# gamma for each color is always at the end of the bin
# thus ensuring that there is no loss of precision;
# it means, however, that the gamma binning is non-uniform which also results in loss of precision;
# but with the help of array ngb2 (see below), the non-uniformity is minimized
def order(phi, gam, wvl):

#   number of gamma bins for the gamma(phi) dependency
#   corresponding to the red color
    ngb1 = 1000

#   number of gamma bins between:
#   red and green
#   green and blue
#   blue and purple
#   trial and error showed that if we have these numbers of
#   bins in between the colors
#   the spacing between gamma points in the eventual grid will be
#   approximately the same (i.e. differing by no more that 5% for the
#   red bins and those in between the colors)
    ngb2 = np.array([4, 4, 4])

#   initiating gamma grid as a list
    gam_grid = []

#   length of gamma bins up to the maximum of gamma(phi) for red color
    dgam = max(gam[3, :]) / ngb1

#   building gamma grid up to the gamma(phi) maximum for red color
    for j in range(0, ngb1 + 1):

        gam_grid.append(j * dgam)

#   length of gamma bins in between the remaining colors
#   this is an array of length 3:
#   red to green
#   green to blue
#   blue to purple
    dg = np.zeros(len(wvl) - 1)

#   loop over all colors except red
    for i in range(len(wvl) - 1, 0, -1):

        k = len(wvl) - 1 - i

#       dividing the distance between the gamma(phi) maxima
#       by the corresponding number of bins given by the ngb2 array
        dg[k] = (max(gam[i - 1, :]) - max(gam[i, :])) / ngb2[k]

#       now we have to incorporate the new bins into the gamma grid that we already have;
#       to achieve that we need to locate the last element in the gamma grid
#       whenever we transition to the next color pair
        last_gam_idx = len(gam_grid) - 1

#       augmenting the gamma grid with the bins corresponding to the current color pair
        for j in range(1, ngb2[k] + 1):

            gam_grid.append(gam_grid[last_gam_idx] + j * dg[k])

#   eventual number of gamma bins
    ngb = ngb1 + ngb2[0] + ngb2[1] + ngb2[2]

#   initiating the array of gamma values corresponding to each gamma bin
    gamma = np.zeros(ngb)

#   calculating the gamma values corresponding to each bin as the mean of bin edges
    for j in range(ngb):

        gamma[j] = (gam_grid[j + 1] + gam_grid[j]) / 2.0

#   gamma grid is ready; we can calculate the relative (to incident) intensity now (in %)
    inten = rainbow.intensity(0, phi, gam, gam_grid, wvl, ngb)

    return gam_grid, gamma, inten

# plotting the gamma(phi) dependence
# which the rainbow calculation is based on
def plot_gamphi(phi, gam, wvl, col):

#   initiate figure (see auxplt.py)
    fig, ax = auxplt.initfig(1, 1, (12.0, 10.0))

#   loop over wavelengths/colors
    for i in range(len(wvl)):

        ax.plot(phi, gam[i, :], color = col[i])

    ax.set_xlim(0, 90)
    ax.set_ylim(0, 90)

    ax.set_xlabel(r'$\phi$, [deg]', fontsize = 30)
    ax.set_ylabel(r'$\gamma$, [deg]', fontsize = 30)

#   set the size of tick labels along both axes
    ax.tick_params(axis = 'both', labelsize = 20)

#   invoke savepdf function from auxplt.py module (auxiliary plotting functions)
    auxplt.savepdf('phigam_0')

# zoomed in version of the gamma(phi) dependence
# the zoom-in is around the angle of maximum deflection
def plot_maxgam(phi, gam, gam_grid, wvl, col):

#   initiate figure (see auxplt.py)
    fig, ax = auxplt.initfig(1, 1, (12.0, 10.0))

#   loop over wavelengths/colors
    for i in range(len(wvl)):

#       plotting the function with symbols (as opposed to lines)
#       to better illustrate the idea of rays going into and out of the drop
#       within different gamma and phi intervals
        ax.scatter(phi, gam[i, :], color = col[i], s = 7)

#   loop over all mean gamma-bin values 
    for j in range(len(gam_grid)):

#       drawing a horizontal line for each value
        ax.axhline(y = gam_grid[j], linewidth = 0.5)

    ax.set_xlim(88, 90)
    ax.set_ylim(80, 84)

    ax.set_xlabel(r'$\phi$, [deg]', fontsize = 30)
    ax.set_ylabel(r'$\gamma$, [deg]', fontsize = 30)

#   set the size of tick labels along both axes
    ax.tick_params(axis = 'both', labelsize = 20)

#   invoke savepdf function from auxplt.py module (auxiliary plotting functions)
    auxplt.savepdf('maxgam_0')

# plot the rainbow color histogram to qualitatively show the separation of colors
def plot_hist(gam, gam_grid, wvl, col):

#   initiate figure (see auxplt.py)
    fig, ax = auxplt.initfig(1, 1, (12.0, 10.0))

#   loop over all wavelengths/colors
    for i in range(len(wvl)):

#       distribute gam values over bins defined by gam_grid values;
#       alpha is the degree of histogram transparency
        ax.hist(gam[i, :], gam_grid, histtype = 'bar', facecolor = col[i], alpha = 0.5)

    ax.set_xlim(0, 90)

    ax.set_xlabel(r'$\gamma$-bins, [deg]', fontsize = 30)
    ax.set_ylabel(r'Number of $\gamma$-values', fontsize = 30)

#   set the size of tick labels along both axes
    ax.tick_params(axis = 'both', labelsize = 20)

#   invoke savepdf function from auxplt.py module (auxiliary plotting functions)
    auxplt.savepdf('hist_0')
