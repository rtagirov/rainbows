import numpy as np
import more_itertools as mit

import importlib
import auxfunc

importlib.reload(auxfunc)

from tqdm import tqdm

# refractive index of water (or steam) as a function of 
# density (kg / m^-3),
# temperature (Celcius),
# wavelength (nm)
# validity:
# 0 < T < 225 Celcius
# 0 < d < 1060 kg / m^-3
# 0.2 < w < 2.5 micrometer
def riH2O(w, d = 1000, T = 20):

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

# angle of deflection for rainbow of a given order
def gamma(order, phi, wvl):

    g = np.zeros((len(wvl), len(phi)))

    for i in range(len(wvl)):

        n = riH2O(wvl[i])

        if order == 1: g[i, :] = 4.0 * np.arcsin(np.sin(phi) / n) - 2.0 * phi
    
    return np.rad2deg(g)

def intensity(order, phi, gam, gam_grid, wvl, ngb):

    p = np.deg2rad(phi)

    G = np.deg2rad(gam)

    grid = np.deg2rad(gam_grid)

    weigh1 = np.zeros((len(wvl), ngb))
    weigh2 = np.zeros((len(wvl), ngb))
    contr1 = np.zeros((len(wvl), ngb))
    contr2 = np.zeros((len(wvl), ngb))

    fren1 = np.zeros((len(wvl), ngb))
    fren2 = np.zeros((len(wvl), ngb))

    for i in range(len(wvl)):

        for j in tqdm(range(ngb), ncols = auxfunc.term_width(), \
                                  desc = 'Calculating intensity ' +
                                         'at ' + str(wvl[i]) + ' nm'):

#        for j in range(ngb):

            p1 = 0.0
            dp1 = 0.0

            p2 = 0.0
            dp2 = 0.0

            idx = np.where((G[i, :] > grid[j]) & (G[i, :] <= grid[j + 1]))
 
            g = [list(gr) for gr in mit.consecutive_groups(idx[0].tolist())]

            if g:

                contr1[i, j] = len(g[0])

                p1 = (p[g[0][0]] + p[g[0][len(g[0]) - 1]]) / 2.0

                dp1 = p[g[0][len(g[0]) - 1]] - p[g[0][0]]

                weigh1[i, j] = np.cos(p1) * dp1

                fren1[i, j] = fresnel(order, p1, wvl[i])

                if len(g) == 2:

                    contr2[i, j] = len(g[1])

                    p2 = (p[g[1][0]] + p[g[1][len(g[1]) - 1]]) / 2.0

                    dp2 = p[g[1][len(g[1]) - 1]] - p[g[1][0]]

                    weigh2[i, j] = np.cos(p2) * dp2

                    fren2[i, j] = fresnel(order, p2, wvl[i])

    for i in range(len(wvl)):

        const = np.sum(weigh1[i, :]) + np.sum(weigh2[i, :])

        weigh1[i, :] /= const
        weigh2[i, :] /= const

    inten0 = contr1 + contr2

    inten = (contr1 * weigh1 * fren1 + contr2 * weigh2 * fren2) * 100.0 / inten0

    return inten0, inten

def fresnel(order, phi, wvl):

    n = riH2O(wvl)

    theta = np.arcsin(np.sin(phi) / n)

    if order == 1:

        t1s = 2.0 * np.cos(phi) / (np.cos(phi) + n * np.cos(theta))

        rs = (np.cos(phi) - n * np.cos(theta)) / ((np.cos(phi) + n * np.cos(theta)))

        t2s = 2.0 * n * np.cos(theta) / (n * np.cos(theta) + np.cos(phi))

        t1p = 2.0 * np.cos(phi) / (n * np.cos(phi) + np.cos(theta))

        rp = (n * np.cos(phi) - np.cos(theta)) / (n * np.cos(phi) + np.cos(theta))

        t2p = 2.0 * n * np.cos(theta) / (n * np.cos(phi) + np.cos(theta))

        Rp = rp**2.0

        T1p = (n * np.cos(theta) / np.cos(phi)) * t1p**2.0
        T2p = (np.cos(phi) / np.cos(theta) / n) * t2p**2.0

        Rs = rs**2.0

        T1s = (n * np.cos(theta) / np.cos(phi)) * t1s**2.0
        T2s = (np.cos(phi) / np.cos(theta) / n) * t2s**2.0

        R = 0.5 * (Rp + Rs)

        T1 = 0.5 * (T1p + T1s)

        T2 = 0.5 * (T2p + T2s)

        return T1 * R * T2
