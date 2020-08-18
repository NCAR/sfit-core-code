    FCIA and SCIA PSEUDO-LINELISTS FOR COLLISION-INDUCED ABSORPTION

                  G. C. Toon and A. Kleinboehl
                   Jet Propulsion Laboratory


INTRODUCTION

This file gives an introduction to the latest version (April 20, 2006) of
the pseudo-linelists for collision induced absorption. Separate pseudo-linelists
are given for foreign-collision induced absorption (FCIA) and self-collision
induced absorption (SCIA). They contain pseudolines for O2 and N2. In the
mid-infrared they cover the fundamental O2 band around 1550 cm-1, the
fundamental N2 band around 2330 cm-1, and the first overtone band of N2
around 4630 cm-1. In near-infrared pseudo-lines are given for the O2
bands around 7900, 9400 cm-1, and 13250 cm-1. Each linelist consists
of 5909 pseudolines with a spacing of 1 cm-1.


ABSORPTION CALCULATIONS

The pseudo-linelists are in the HITRAN format. However, as they describe
collision induced absorption, the line strength is given in the unit
cm-1/molec2/cm-5. To calculate the correct absorption, the SCIA absorption
coefficients derived from the pseudolines must be multiplied by the number
density times the volume mixing ratio of the considered gas, and the FCIA
absorption coefficients must be multiplied by the number density times
(1 - volume mixing ratio).
   Absorption = S_scia * d^2 * v^2  +   S_fcia * d^2 * v*(1-v)
The sum of both contributions can then be used to
calculate optical thicknesses or transmittances. Note that the units of the
line strengths and hence the calculation of the absorption is different
from previous versions of the CIA linelists, in which units were given in
cm-1/molec/cm-2/atm and absorption was calculated by multiplicating
the coefficients with the pressure.


DERIVATION OF THE LINELISTS

For the O2 and N2 fundamental bands tabulated values from empirical
models were converted to HITRAN units (line strength and ground state
energy). For the N2 overtone band absorptions were read from a figure
and converted to HITRAN units. In all three cases it was necessary to
fit an exponential function to the linear temperature dependence of
(T/T_hitran), which originates from the rotational partition function.
This fit has an accuracy of 2% in a temperature range between 190-295 K.

The O2 fundamental band is based on the measurements and the empirical model
given in Tab. 1 by Thibault et al. (1997), which gives parameters for O2-air
collision induced absorption. Separation in FCIA and SCIA was done by assuming
a temperature independent effectiveness of
O2-O2 collisions = 0.94 * effectiveness of O2-N2 collisions
estimated from data in Tab. 1 of Orlando et al. (1991). The linelist has been
extrapolated to zero line strengths at 2*FWHM from the maximum absorption
towards lower wavenumbers and 2.5*FWHM towards higher wavenumbers. This covers
the range from 1275-1905 cm-1, giving a total of 631 pseudo-lines.

The N2 fundamental band is based on the measruements and the empirical model
given in Tab. 1 by Lafferty et al. (1996). They give their empirical model for
N2-N2 collisions, which could be directly converted to HITRAN units for SCIA.
N2-O2 collisions were calculated from the N2-N2 collisions using their
temperature dependent conversion factor, which parameterizes a linear
temperature dependence based on Menoux et al. (1993). For this, an exponential
function has been fitted to the product of (T/T_hitran) and the N2-O2
conversion factor, the accuracy of which is 0.3% between 190-295 K. The new
linelist has been extrapolated to zero line strengths at 2*FWHM from the
maximum absorption towards lower wavenumbers and 2.5*FWHM towards higher
wavenumbers. This covers the range from 2030-2705 cm-1, giving a total of
676 pseudo-lines.

The N2 first overtone band is based on absorption data given by Shapiro and
Gush (1966). The absorption was digitized from Fig. 8 of the paper, and
converted to pseudoline strengths assuming similar ground state energies as
for the fundamental band. Shapiro and Gush (1966) give the absorption for
N2-N2 collisions. Separation to SCIA and FCIA was done assuming an
effectiveness of N2-O2 collisions = 0.92 * effectiveness of N2-N2 collisions
as suggested by Menoux et al. (1993) for the fundamental N2 band (Tab. 4),
however, no possible temperature dependence has been taken into account.
The new linelist covers a range of 4330-4930 cm-1 with 601 pseudo-lines.

For the O2 CIA in the 7900 and 9400 cm-1 bands pseudolines were derived
by fitting laboratory spectra from Smith and Newnham (2000). These lab
spectra covered 198 to 295K and O2 vmrs from 21% to 75%. The pseudo-lines
are able to reproduce the O2 slant column amounts in the lab spectra to
within 1% rms (2.4% max).

The O2 A-band CIA pseudo-linelist was derived by fitting the absorptions
given in Fig. 7 of the paper "Line-Mixing and collision-induced absorption
by O2 in the A-band" by H.Tran, C. Boulet, and J.-M.Hartmann.  Since no
lab spectra were available, we simply fitted the two curves in the paper.


ACKNOWLEDGMENTS

We greatly acknowledge the work by F. Mills who created earlier versions
of CIA-linelists for the O2 and N2 fundamental bands and whose excellently
documented routines made the production of the new linelists an almost easy
task. We also would like to thank P. Wennberg for digitizing the figure of the
N2 first overtone absorption and providing a first pseudo-linelist as well as
transmission calculations for a balloon flight.


REFERENCES

Shapiro, M. M., and H. P. Gush, The collision-induced fundamental and first
overtone bands of oxygen and nitrogen, Can. J. Phys., 44, 949-963, 1966.

Orlando, J. J., G. S. Tyndall, K. E. Nickerson, and J. G. Calvert,
The temperature dependence of collision-induced absorption by oxygen
near 6 um, J. Geophys. Res., 96, 20755-20760, 1991.

Menoux, V., R. Le Doucen, C. Boulet, A. Roblin, and A. M. Bouchardy,
Collision-induced absorption in the fundamental band of N2: temperature
dependence for N2-N2 and N2-O2 pairs, Appl. Opt., 32, 263-268, 1993.

Lafferty, W. J., A. M. Solodov, A. Weber, W. B. Olson, and J.-M. Hartmann,
Infrared collision-induced absorption by N2 near 4.3 um for atmospheric
applications: measurements and empirical modeling, Appl. Opt., 35, 5911-5917,
1996.

Thibault, F., V. Menoux, R. Le Doucen, L. Rosenmann, J.-M. Hartmann, and
C. Boulet, Infrared collision-induced absorption by O2 near 6.4 um for
atmospheric applications: measurements and empirical modeling, Appl. Opt., 36,
563-567, 1997.

Smith, K.M. and Newnham D.A., Near infrared absorption cross-sections and
integrated absorptions of molecular oxygen (O2, O2-O2, and O2-N2),
JGR, 105, 7383-7396, 2000 

Tran, H., C. Boulet, and J.-M. Hartmann, Line-mixing and collision-indiced
absorption by O2 in the A-band, submitted to JQSRT?
