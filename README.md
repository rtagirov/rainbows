This version calculates first order rainbow.
The calculation is done for purple, blue, green and red colors.
The results of the original paper are reproduced only to some extent.
In particular, for reasons unknown, there is more green in the rainbow than blue.
The values of relative intensity are a bit higher than in the original paper.

rainbows.py --- main program

rainbow.py  --- various rainbow functions

                --> riH2O - refractive index of water (see 'papers' directory)

                --> gamma - angle of deflection for a rainbow of a given order

                --> intensity - calculation of intensity from the gamma(phi) function

                --> fresnel - transmission and reflection coefficients for unpolarized light given by Frensnel equations

first.py    --- builds optimal gamma grid for the first order rainbow and calls 'intensity' from rainbow.py

auxfunc.py  --- auxiliary functions

auxplt.py   --- auxiliary functions for plotting
