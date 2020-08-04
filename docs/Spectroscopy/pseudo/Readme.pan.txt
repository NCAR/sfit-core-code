
             PAN PSEUDO-LINELIST

        A. Kleinboehl and G. C. Toon


INTRODUCTION.
This document gives information on the pseudo-linelist for
PAN (CH3C(O)OONO2) derived at JPL in October 2005.
The linelist was created based on 5 data sets of absorption
cross-sections from laboratory measurements taken at
the Rutherford Appleton Laboratory (RAL). The measurements and
the absorption cross sections are described in Allen et al. (2005a,b).

The cross-section data sets were derived from spectra measured in
regions between 550 and 2200 cm-1 at three different temperatures
and pressures between 0.19 and 2.20 mbar. Each measurement used
the same cell of 26.1 cm length. The resolution of the spectra
was 0.25 cm-1, and they were given with a spectral point spacing
of about 0.1 cm-1. The temperatures and frequency ranges for the
different measurement sets are given in the following table.

 #  Temp [K] Freq.-range [cm-1]
-------------------------------
 1    295         550-1650
 2    295        1650-2200
 3    273         550-1400
 4    250         550-1400
 5    250        1600-2200
-------------------------------


DESCRIPTION.
The cross-sections were converted back into transmittance
spectra using the cell length and an average gas concentration.
The resulting laboratory transmittance spectra were then
simultaneously fitted (using the GFIT algorithm) by iteratively
adjusting the strengths and ground-state energies of the pseudo-lines.

Due to the resolution of the laboratory spectra of 0.25 cm-1
a pseudo-line spacing of 0.1 cm-1 was chosen. Fitting was performed
in the bands with centers at 794 cm-1, 1163 cm-1, 1302 cm-1, and 1741 cm-1,
which correspond to the four strongest bands of PAN in the mid-infrared.
Because there were still small residuals of water vapor contamination in
the RAL spectra the water vapor amount was fitted in a narrow frequency window
and taken into account during the fitting of the PAN bands.
At each line frequency, an effective strength and ground-state energy
was derived by simultaneous non-linear least squares fitting to the spectra
that covered the frequency band of interest. In analogy to other organics
(e. g. CH3OH) an airbroadened halfwidth of 0.1 cm-1/atm and a selfbroadened
halfwidth of 0.4 cm-1/atm was assumed. The result of the fitting process is
a pseudo-linelist that covers three frequency regions, 760-870 cm-1,
1110-1340 cm-1, and 1700-1780 cm-1. It contains 4203 pseudo-lines.


PARTITION FUNCTION.
The rotational partition function for PAN was assumed to be (296/T)^2.
The vibrational partition function was calculated in the way it had been
done for the ATMOS experiment, as described e. g. by Norton and Rinsland (1991).
The following 27 vibrational frequencies (in cm-1) were used (J. Francisco,
private communication), all degeneracies were set to 1:

3164, 3121, 3058, 1880, 1806, 1475, 1471, 1400, 1352, 1172, 1065,
999, 984, 828, 806, 736, 727, 616, 585, 495, 373, 327, 316, 100,
96, 82, 24

It should be noted that the derived ground-state energies are higher than
for most other molecules. This is likely to be related to the large number
of low-lying vibrational modes, such that very few PAN molecules are in
the ground-state at the temperatures at which the spectra were measured.


ACCURACY.
To estimate how well the pseudo-linelist represents the RAL spectra,
test retrievals were performed in which the laboratory spectra were fitted
using the pseudo-linelist. The retrieved scale factors for the PAN
abundances in the different spectra are tabulated below.

         Scale factors retrieved in freq. region
 #	760-870 cm-1  1110-1340 cm-1  1700-1780 cm-1
-----------------------------------------------------
 1         1.0089         1.0048 
 2                                        0.9989
 3         0.9903         0.9947
 4         1.0036         1.0022
 5                                        1.0005
-----------------------------------------------------

The pseudolines correctly represent the RAL spectra to better than 1%
of the given PAN amount in all bands. However, it has to be noted that
additional uncertainties may arise from the extrapolation of the temperature
dependence of the pseudo-lines from laboratory temperatures to atmospheric
temperatures. This may in particular affect the pseudo-linelist in the
region from 1700-1780 cm-1 because it was created with spectra available
at two different temperatures temperatures only.


ACKNOWLEDGMENTS.
We would like to thank Grant Allen for providing the files with the PAN
cross sections, and we are grateful to Joe Francisco for providing the
fundamental frequencies for PAN prior to publication.


REFERENCES.

Allen, G., J. J. Remedios, D. A. Newnham, K. M. Smith, and P. S. Monks,
Improved mid-infrared cross-sections for peroxyacetyl nitrate (PAN) vapour,
Atmos. Chem. Phys., 5, 47-56, 2005a.

Allen, G., J. J. Remedios, and K. M. Smith, Low temperature mid-infrared
cross-sections for peroxyacetyl nitrate (PAN) vapour,
Atmos. Chem. Phys. Discuss., 5, 5669-5685, 2005b.

Norton, R. H. and C. P. Rinsland, ATMOS data processing and science
analysis methods, Appl. Opt., 30, 389-400, 1991. 
