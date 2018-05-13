import numpy as np
import sys
import importlib
import rainbow

importlib.reload(rainbow)

def order(phif, gamf, wvl, col):

    ngb1 = 1000
    ngb2 = np.array([12, 12, 11])

    gam_grid = []

    dgam = max(gamf[0, :]) / ngb1

    for j in range(0, ngb1 + 1):

        gam_grid.append(j * dgam)

    dg = np.zeros(len(wvl) - 1)

    for i in range(0, len(wvl) - 1):

        dg[i] = (max(gamf[i + 1, :]) - max(gamf[i, :])) / ngb2[i]

        last_gam_idx = len(gam_grid) - 1

        for j in range(1, ngb2[i] + 1):

            gam_grid.append(gam_grid[last_gam_idx] + j * dg[i])

    ngb = ngb1 + ngb2[0] + ngb2[1] + ngb2[2]

    gamma = np.zeros(ngb)

    for j in range(ngb):

        gamma[j] = (gam_grid[j + 1] + gam_grid[j]) / 2.0

    inten0, inten = rainbow.intensity(1, phif, gamf, gam_grid, wvl, ngb)

    inten0[np.where(inten0 == 0.0)] = np.nan

    inten[np.where(inten == 0.0)] = np.nan

    return gam_grid, gamma, inten0, inten
