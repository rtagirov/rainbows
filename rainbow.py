import numpy as np
import more_itertools as mit

import importlib
import auxfunc

importlib.reload(auxfunc)

#package for showing progress bars
from tqdm import tqdm

# calculation of relative (to incindent) intensity
# at each wavelength for a rainbow of order k
# from the gamma(phi) dependence
# gam_grid --- gamma binning
# ngb --- number of gamma bins (calculated in zero.py/first.py)
def intensity(k, phi, gam, gam_grid, wvl, ngb):

# phi to radians
    p = np.deg2rad(phi)

# gamma to radians
    G = np.deg2rad(gam)

# gamma bins grid to radians
    grid = np.deg2rad(gam_grid)

# contribution to intensity from the phi intervals
# corresponding to a given gamma bin
# for the rainbows of the zero, first and second order the
# maximum number of phi intervals is two
    contr1 = np.zeros((len(wvl), ngb))
    contr2 = np.zeros((len(wvl), ngb))

# associated geometrical weights for the contributions
    weigh1 = np.zeros((len(wvl), ngb))
    weigh2 = np.zeros((len(wvl), ngb))

# associated Fresnel weights for the contributions
    fren1 = np.zeros((len(wvl), ngb))
    fren2 = np.zeros((len(wvl), ngb))

# loop over wavelengths/colors
    for i in range(len(wvl)):

# loop over gamma bins
# tqdm creates the progress bar of width corresponding to the width of the terminal
# with the given description (desc) for the progress bar
        for j in tqdm(range(ngb), ncols = auxfunc.term_width(), \
                                  desc = 'Calculating intensity (order ' \
                                          + str(k) + ') ' +
                                         'at ' + str(wvl[i]) + ' nm'):

            p1 = 0.0  # middle phi for the first phi interval
            dp1 = 0.0 # length of the first phi interval

            p2 = 0.0  # middle phi for the second phi inteval
            dp2 = 0.0 # length of the second phi interval

#           indices of all gammas falling into a given gamma bin
#           idx is a tuple, numpy array of indices is given by the first element of this tuple idx[0]
            idx = np.where((G[i, :] > grid[j]) & (G[i, :] <= grid[j + 1]))
 
#           here we take the indices we found and divide them into two groups
#           by means of applying the method from more_itertools module ('mit', imported above) called consecutive_groups;
#           this method takes a list (as opposed to a numpy array) as an input;
#           therefore idx[0], which is a numpy array, must be transformed into list;
#           this can be done with the python method tolist();
#           so eventually idx[0].tolist() is the list of all gamma indices within a given gamma bin
#           then we take this list and apply the consecutive_group method which creates a nested list of groups of indices,
#           i.e. if in the original list idx[0].tolist() we find that two adjacent indices differ by more than 1
#           it means that we have found where the first group of indices ends and the next one begins
#           the groups correspond to different phi intervals rays in which end up in the same gamma bin
#           to sum up: g is a nested list of index groups of the form 
#           g = [[idx_11, idx_12, ... ,idx_1n], [idx_21, idx_22, ... ,idx_2n]];
#           the second list in g can be empty depending on which gamma bin one looks at 
#           as, for example, evidenced by the gamma(phi) plot for the first order rainbow;
#           both lists in g can be empty as well, e.g. when we consider purple color in
#           the first order rainbow and any gamma bin that is located higher
#           than the maximum gamma for the purple color
            g = [list(gr) for gr in mit.consecutive_groups(idx[0].tolist())]

#           this condition means: if list g is NOT empty, i.e. g is NOT equal to [[]]
            if g:

#               g[0] is the group of indices corresponding to the first phi interval
#               contribution is proportional to the number of rays within the phi interval,
#               which is the length of the first list in g, i.e. len(g[0])
                contr1[i, j] = len(g[0])

#               phi corersponding to the first interval is the mean phi of the interval
#               g[0][0] -- first index in the first group
#               g[0][len(g(0)) - 1] --- last index in the first group
#               accordingly:
#               p[g[0][0]] --- first phi angle in the first group of phi angles
#               p[g[0][len(g(0)) - 1]] --- last phi angle in the first group of phi angles
                p1 = (p[g[0][0]] + p[g[0][len(g[0]) - 1]]) / 2.0

#               length of the phi interval is the difference between the last and first phi angle
                dp1 = p[g[0][len(g[0]) - 1]] - p[g[0][0]]

#               geometrical weight corresponding to the first phi interval, see main_article.pdf, Eq. (4)
                weigh1[i, j] = np.cos(p1) * dp1

#               Fresnel weight for a rainbow of a given order at the given wavelength
#               correspoinding to the first phi interval
#               the function fresnel is declared below
                fren1[i, j] = fresnel(k, p1, wvl[i])

#               this condition means: if list g has two elements,
#               i.e. if there are two groups of indices
#               the handling of the second phi interval if it exists
#               is analogous to the first phi interval
                if len(g) == 2:

                    contr2[i, j] = len(g[1])

                    p2 = (p[g[1][0]] + p[g[1][len(g[1]) - 1]]) / 2.0

                    dp2 = p[g[1][len(g[1]) - 1]] - p[g[1][0]]

                    weigh2[i, j] = np.cos(p2) * dp2

                    fren2[i, j] = fresnel(k, p2, wvl[i])

#   the sum of weights has to be equal to one
#   because of the numerical innacuracies the sum of
#   the calculated geometrical weights slightly differs from one
#   here we perform the weight renormalization to correct this
    for i in range(len(wvl)):

        const = np.sum(weigh1[i, :]) + np.sum(weigh2[i, :])

        weigh1[i, :] /= const
        weigh2[i, :] /= const

#   incident intensity is proportional to the sum of the contributions
    inten0 = contr1 + contr2

#   we replace the zero values of incident intensity with NaNs
#   so that we will not get the 0.0 / 0.0 division when calculating the outgoing intensity;
#   if the zero values are replaced with NaNs we will have 0.0 / NaN division which gives NaN
#   and is not seen in the final intensity plot;
#   this is what we want, plus 0.0 / NaN division is not
#   frowned upon by the Python interpreter
    inten0[np.where(inten0 == 0.0)] = np.nan

#   outgoing intensity is the sum of the weighted contributions;
#   we are concerned with the intensity relative to the incident one
#   therefore we divide the sum of the weighted contributions
#   by the incident intensity and multiply by 100 to get the 
#   relative intensity in percents of the incident intensity
    inten = (contr1 * weigh1 * fren1 + contr2 * weigh2 * fren2) * 100.0 / inten0

    inten[np.where(inten == 0.0)] = np.nan

    return inten

# Fresnel weight for a rainbow of a given order as a function of incidence angle and wavelength;
# the weight is given by a product of transmittance and reflectance;
# depending on the rainbow order the product includes different components;
# their number also depends on the rainbow order;
# transmission and reflection coefficients are first calculated for
# p-polarized and s-polarized light (see fresnel_eqs.pdf)
# for the corresponding number of transmissions and reflections
# then the (p + s) / 2 mean is calculated for each coefficient
# to describe the unpolarized light
def fresnel(k, phi, wvl):

#   refraction coefficient for a given wavelength
    n = riH2O(wvl)

#   angle of refraction
    theta = np.arcsin(np.sin(phi) / n)

#   s-polarization (see fresnel_eqs.pdf, p. 13, left column)
#   ----------------------------------------------------------------------
#   first transmission coefficient (into the drop)
    t1s = 2.0 * np.cos(phi) / (np.cos(phi) + n * np.cos(theta))

#   reflection coefficient
    rs = (np.cos(phi) - n * np.cos(theta)) / ((np.cos(phi) + n * np.cos(theta)))

#   second transmission coefficient (out of the drop)
    t2s = 2.0 * n * np.cos(theta) / (n * np.cos(theta) + np.cos(phi))

#   reflectance (see fresnel_eqs.pdf, p. 17)
    Rs = rs**2.0

#   transmittance (see fresnel_eqs.pdf, p. 18)
    T1s = (n * np.cos(theta) / np.cos(phi)) * t1s**2.0
    T2s = (np.cos(phi) / np.cos(theta) / n) * t2s**2.0

#   p-polarization (see fresnel_eqs.pdf, p. 13, right column)
#   ----------------------------------------------------------------------
#   first transmission coefficient (into the drop)
    t1p = 2.0 * np.cos(phi) / (n * np.cos(phi) + np.cos(theta))

#   reflection coefficient
    rp = (n * np.cos(phi) - np.cos(theta)) / (n * np.cos(phi) + np.cos(theta))

#   second transmission coefficient (out of the drop)
    t2p = 2.0 * n * np.cos(theta) / (n * np.cos(phi) + np.cos(theta))

#   reflectance (see fresnel_eqs.pdf, p. 17)
    Rp = rp**2.0

#   transmittance (see fresnel_eqs.pdf, p. 18)
    T1p = (n * np.cos(theta) / np.cos(phi)) * t1p**2.0
    T2p = (np.cos(phi) / np.cos(theta) / n) * t2p**2.0
#   ----------------------------------------------------------------------

#   calculating the (p + s) / 2 means of reflectance and transmittance to describe the unpolarized light
    R = 0.5 * (Rp + Rs)

    T1 = 0.5 * (T1p + T1s)
    T2 = 0.5 * (T2p + T2s)

#   fraction of light that gets out of the drop after two transmissions and k reflections
    return T1 * R**k * T2

# angle of deflection as a function of wavelength/color and angle of incidence
# for rainbow of a given order
# phi is given in radians
def gamma(order, phi, wvl):

# gamma is different for each wavelength/color, hence 2D array
    g = np.zeros((len(wvl), len(phi)))

# loop over wavelengths/colors
    for i in range(len(wvl)):

# refraction index for given wavelength/color
        n = riH2O(wvl[i])

        if order == 0:

#           follows from geometrical considerations
            g[i, :] = 2.0 * (phi - np.arcsin(np.sin(phi) / n))

        if order == 1:

#           Eq. (3) in main_article.pdf
            g[i, :] = 4.0 * np.arcsin(np.sin(phi) / n) - 2.0 * phi

        if order == 2:

#           follows from geometrical considerations
            g[i, :] = np.pi + 2.0 * phi - 6.0 * np.arcsin(np.sin(phi) / n)

# returning the value of gamma in degrees 
    return np.rad2deg(g)

# refractive index (ri) of water (or steam) as a function of 
# wavelength (nm), water density (kg / m^-3) and temperature (Celcius)
# see ref_ind_water.pdf in the 'papers' folder, Eq. (7)
# validity:
# 0 < T < 225 Celcius
# 0 < d < 1060 kg / m^3
# 0.2 < w < 2.5 micrometer
# here, default values of water density (1000 kg / m^3) and temperature (15 deg Celsius) are taken
def riH2O(w, d = 1000, T = 15):

# nm -> micrometers
    w *= 1e-3

# Celsius -> Kelvins
    T += 273.15

    d0 = 1e+3

    T0 = 273.15

    w0 = 0.589

    d /= d0

    T /= T0

    w /= w0

    a = np.array([0.243905091,
                  9.53518094e-3,
                 -3.64358110e-3,
                  2.65666426e-4,
                  1.59189325e-3,
                  2.45733798e-3,
                  0.897478251,
                 -1.63066183e-2])

    wuv = 0.2292020

    wir = 5.432937

    c = d * (a[0] +
             a[1] * d +
             a[2] * T +
             a[3] * w**2 * T +
             a[4] / w**2 +
             a[5] / (w**2 - wuv**2) +
             a[6] / (w**2 - wir**2) +
             a[7] * d**2)

    return np.sqrt((2.0 * c + 1) / (1 - c))
