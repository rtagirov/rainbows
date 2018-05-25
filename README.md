This version calculates zero and first order rainbows.
The calculation is done for purple, blue, green and red colors.
The results of the original paper, main_article.pdf, i.e. for the first order rainbow, are reproduced only to some extent.
In particular, for reasons unknown, there is more green in the first order rainbow than blue.
Moreover, the values of relative intensity are a bit higher than in main_article.pdf.

rainbows.py --- main program; calls functions from the module zero.py/first.py

                plots:

                --> gamma(phi) dependence

                --> color histogram

                --> relative (to incident) intensity

rainbow.py  --- various rainbow functions

                --> intensity - derivation of relative (to incident) intensity from the gamma(phi) function

                --> fresnel - reflectance and transmittance for unpolarized light given by Frensnel equations

                --> gamma - angle of deflection for a rainbow of a given order

                --> riH2O - refraction index of water (see 'papers' directory)

first.py    --- builds optimal gamma grid for the first order rainbow, calls 'intensity' from rainbow.py;
                plotting functions for gamma(phi) dependence and the color histogram are declared here

auxplt.py   --- plotting auxiliary functions

auxsys.py   --- system auxiliary functions

auxfunc.py  --- miscellaneous auxiliary functions

The Plan:

1. Understanding the main idea of the paper (main_article.pdf), i.e. transformation of gamma(phi) dependence into relative intensity

2. Understanding the main workings of the existing code

3. Updating the code so that it describes the zero order rainbow formation

   This will involve:

        a. geometrical considerations of the zero order rainbow formation

        b. calculation of gamma for zero order rainbow using function gamma in rainbow.py,
           then invoking the zero order gamma in rainbows.py

        c. writing the module zero.py based on the module first.py,
           i.e. gamma grid optimization for the zero order rainbow in zero.py
           (analogous to how it is done in first.py) plus the plots:

           c1. plotting gamma(phi) for the zero order rainbow

           c2. plotting the color histogram for the zero order rainbow

        d. writing the zero order branch for the 'fresnel' function in rainbow.py

        e. invoking zero.py from rainbows.py

        e. adding the intensity of the zero order rainbow to the existing plot for the first order intensity

4. Repeating steps in 3 for the second order rainbow
   (creation of second.py module + second order branch in function fresnel from rainbow.py module)

5. Plot the frensnel weights for all rainbow orders as functions of phi and/or gamma
