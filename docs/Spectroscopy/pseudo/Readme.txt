      DERIVATION OF PSEUDO-LINES FROM LABORATORY CROSS-SECTIONS

                        G.C.Toon, B.Sen, and A.Kleinboehl
                     Jet Propulsion Laboratory

INTRODUCTION
The problem of how best to interpolate/extrapolate laboratory cross-
section spectra in temperature and pressure has no simple solution.
We have chosen the approach of deriving "pseudo-linelists" to represent
the absorption of heavy molecules (e.g. CCl4, CFC-11, CFC-12, CFC-113,
HCFC-22, ClNO3, HCFC-142b, N2O5, CF4, SF6) for which only cross-section
spectra exist. These pseudo-linelist have been sucessfully used in the
analysis of MkIV balloon spectra (e.g., Toon et al., JGR, 104, 26,779,
1999) and ATMOS Version 3 data (Irion et al., 2002).

This directory contains various pseudo-linelists which were derived
at JPL by performing spectral fits to laboratory transmittance spectra.
These lab spectra were NOT measured here at JPL, they were from various
sources that are tabulated below. The lab transmittance spectra were
re-created using the published temperature- and pressure-dependent
cross-sections, and knowledge of the ILS of the laboratory FTIR
spectrometer, and cell length.

These pseudo-linelists are in the HITRAN format and can be used to compute
spectra in the same manner as any other linelist. The directory also
contains the file, "isotopomer.dat", of molecular parameters (vibration
frequencies, temperature-dependences of the rotational partition function)
assumed in the derivation of the pseudo-linelists, along with a FORTRAN
program "isotopomer.f" defining the paramaters and illustrating how to
read them.

Each pseudo-linelist was derived by fitting all of the relevant
laboratory transmittance spectra simultaneously while solving for
the 296K strength and the Ground State Energy (E") of each pseudo-
line. The pressure-broadened half-width (PBHW) and its temperature
dependence were determined "manually", by trying various values
and selecting the ones that gave the best overall fit.  Generally,
for gases without sharp absorption features, the goodness of fit
was insensitive to the choice of PBHW, whereas for gases like
CFC-12 and HCFC-22 which have sharp Q-branches, the right choice
of PBHW is important.  Note that for some gases (e.g. CFC-12) the
resulting value for the temperature-dependence of the PBHW (0.0)
is well outside the normal range (0.5 to 0.8).  All lines in a
given absorption band were assumed to have the same PBHW and
temperature dependence.

The idea of using pseudo-lines to represent broad featureless
absorption bands is not new. However, whereas previously workers
minimized the number of lines needed by ascribing them an
exaggerated PBHW, we achieve the same goal by giving each
pseudo-line an exaggerated Doppler width.  The advantage of this
latter approach is that it allows the correct PBHW to be employed,
so that a realistic pressure-dependence can still be simulated,
even in cases when all of the laboratory transmittance spectra
were measured at low pressure (e.g. CF4).

These lists are not intended to supplant proper quantum-mechanically-
based linelists. They were derived primarily as a convenient means of
interpolating (and extrapolating) the laboratory cross-sections to
temperatures and pressures where actual measurements are unavailable
(I could not think of a realistic way of doing this directly from the
cross-sections). However, in deriving and using these pseudo-linelists,
several additional advantages became apparent:
1) Since the pseudo-linelists are in the HITRAN format, they can be
   accessed in exactly the same manner as all the regular gases, avoiding
   special code to read the raw cross-section spectra and interpolate
   them to the desired temperatures and pressures.
2) Fitting a physically-based function to the laboratory cross-sections
   also serves as a quality control measure: Since we are typically trying
   to determine just two unknowns (S & E") from 4-30 spectra, the problem
   is over-determined and so performing the fit provides an assessment of
   the consistency of the various laboratory transmittance spectra. This
   makes it possible to identify and reject laboratory spectra which are
   inconsistent with the others, or even to quantify biases between
   different sets of laboratory spectra, perhaps measured under very
   different conditions.  Furthermore, the retrieval of unphysical
   (i.e. -ve) values of S and E" provides a warning that serious problems exist.
3) Fitting the laboratory spectra provides an opportunity to remove
   instrumental artifacts. For example, the laboratory cross-sections
   will always be convolved with the Instrument Line Shape (ILS) of
   the laboratory spectrometer.  In making a pseudo-linelist, the
   effects of this ILS is removed, since it is included in the forward
   model which calculates the cross-sections from the pseudo-lines.
   This is particularly important if the laboratory transmittance spectra
   are measured at a worse spectral resolution than the atmospheric spectra.
   Another example of removal of an instrumental artifact, is that of
   channel fringes in the laboratory transmittance spectra. If these
   are properly fitted they cannot propagate into the pseudo-linelist.
4) Fitting the laboratory transmittance spectra provides an opportunity
   to remove absorption lines not belonging to the gas of interest. For
   example, laboratory spectra acquired over the 1300 to 1900 cm-1 region
   are often contaminated by H2O absorption lines. However, by fitting
   H2O along with the gas of interest, during the analysis of the
   laboratory spectra, propagation of H2O artifacts into the pseudo-
   linelist of the gas of interest will be greatly reduced.
5) Several different laboratory data-sets, even with widely different
   measurement conditions and spectral resolutions, can easily be
   assimilated into a single pseudo-linelist.
6) At the end of the fitting process, the pseudo list can be checked by
   comparing the forward model calculation (which uses the pseudo-lines)
   with the laboratory transmittance spectra. Of course, the agreement
   will not be perfect since the fit was overconstrained, but the
   differences are usually <1%.
7) Absorption spectrum derived from the pseudolines are guaranteed to
   be continuous and differentiable functions of pressure and temperature
   (unlike some bivariate interpolation schemes), which helps minimize
   artifacts in the retrieved vmr profiles.
8) Since all the pseudo-lines in a given band are assumed to have the
   same PBHW and Doppler widths, only one evaluation of the Voigt
   lineshape per atmospheric level is necessary to compute the
   absorption spectrum resulting from all the psuedo-lines (provided
   that this lineshape is stored). Thus, the speed of using the pseudo-
   linelists is competitive with 2-D interpolation in the raw cross-
   sections (assuming one knew a good way of doing this).

SPACING OF PSEUDO-LINES
The choice of line spacing for the pseudo-lines was somewhat arbitrary.
We tried to make it as wide as possible to minimize the total number of
lines, yet still preserve any structure observed in the laboratory
transmittance spectra that would also be apparent in an atmospheric
spectrum.

Typically, the line spacing was chosen to be similar to the resolution
of the laboratory transmittance spectra.  Note that the positions and
spacing of the pseudo-lines are completely independent of the spectral
frequencies in the laboratory spectra. This fact makes it possible to
simultaneously fit different sets of laboratory spectra.

Most of the pseudo-linelists are spaced at 0.01 cm-1, which is ten times
larger than an actual Doppler widths of most heavy gases. Although this
would not be a problem in the troposphere where the pressure broadening
would cause the pseudo-lines to overlap, in the upper stratosphere a high
resolution computed spectrum would show narrow lines with large gaps
between. To avoid this problem one must artificially increase the Doppler
width until it approximately matches the line spacing. A convenient way
of doing this is to set the molecular weight to an artificially small
value (e.g. 1), however, this has drawbacks if one wants to use pseudo-
lines, together with real quantum-mechanical lines, of the same gas in
the same interval. Therefore, in all of the pseudo linelists, we have
defined the isotope number to be zero.  This allows pseudo lines to be
easily distinguished from real lines (which have isotope numbers 1-9),
and could allow the line-by-line code to explitly set the Doppler width
equal to the pseudo-line spacing whenever it encounters pseudo-lines,
avoiding the need to fudge the molecular weight. This is especially
helpful for gases like ClNO3 for which a proper quantum-mechanically
derived linelist (requiring the actual molecular weight) exists for the
region around the 780 cm-1 Q-branch, but pseudo lines must be used for
other regions.

LINE STRENGTH
The following expression for line strength was assumed in the derivation
of the pseudo-linelists:

S(T)=S(296).(296/T)^tdrpf.Qvib(T)/Qvib(296).SE(T)/SE(296).exp(hcE"(1/296-1/T)

where S(296) is the line strength (cm-1/molec/cm-2) at 296K.
tdrpf is the Temperature Dependence of the Rotational Partitiion Function,
and is commonly denoted by the symbol "Beta".
Qvib(T) = II [1-exp(-h.c.Vj/kT)]^(-Gj) is the vibrational partition function
and the product is performed over all the vibrational frequencies, Vj, which
are read from the file "isotopomer.dat", along with their degeneracies, Gj.
SE(T) = [1-exp(-h.c.Vi/kT)]  is the correction for the Stimulated Emission,
Vi being the center frequency (cm-1) of the line in question.
The term (296/T)^tdrpf is commonly known as the rotational partition function
and tdrpf is usually (1.0, 1.5, or 2.0)
The term exp(hcE"(1/296-1/T) is simply the Boltzmann factor, E" being the
ground-state energy (cm-1).

Most forward models should already have the code to compute the above
expression because it is needed for the lighter gases. So it should be
a simple matter to extend this capability to the heavy gases. Note that
the same expression for S(T) was used for the fitting of the laboratory
transmittance spectra, so the derived values of S(296) and E" values will
only correctly reproduce the laboratory spectra provided that the user
employs the same expressions.

Pseudo-linelists for collision induced absorption are given for foreign-
collision induced absorption (fcia) and self-collision induced absorption
(scia). They contain pseudo-lines for O2 and N2 with the line strength
units (cm-1/molec2/cm-5). Note that in the case of CIA calculations,
the scia absorption coefficients derived from the pseudolines must be
multiplied by the number density times the volume mixing ratio of the
considered gas, and the fcia absorption coefficients must be multiplied
by the number density times (1 - volume mixing ratio). The sum of both
contributions can then be used to calculate optical thicknesses or
transmittances.

Finally, we want to make it clear that one should not expect the forward 
calculation made using these pseudo-lines to agree perfectly with individual
laboratory transmittance spectra, since the pseudo-lines were derived from
an over-determined fit to ALL of the laboratory spectra. Differences will
arise from noise on the laboratory spectra, and uncertainties in the measurement
conditions (T,P,vmr), in addition to inadequacies in the pseudo-line approach.

In the table below we summarize the gases and spectral intervals for which
we have computed pseudo-linelists. We estimated the maximum error in absorber
found in re-fitting the laboratory spectra using the final pseudo-linelist.
Note that at specific frequencies, the error in the computed absorption
coefficient may well exceed these tabulated values.

 File      GAS       Interval   Spacing   Lines  Error  Measurer
============================================================================
cf4.h92    CF4       1250 - 1290   0.0025  16093   4%   Nemtchinov & Varanasi

f12.h92    CFC-12     850 -  950   0.010   10000   2%   Varanasi
           CFC-12    1050 - 1200   0.010   15000   1%   Varanasi

f11.h92    CFC-11     810 -  880   0.010    7000   7%   Varanasi
           CFC-11    1050 - 1120   0.010    7000   6%   Varanasi

ccl4.h92   CCl4       750 -  812   0.010    6201   3%   Varanasi & Nemtchinov

f22.h92    CHF2Cl     776 -  850   0.00742  9977   5%   Varanasi & McDaniel
           CHF2Cl    1080 - 1150   0.010    7001   2%   McDaniel
           CHF2Cl    1290 - 1335   0.010    4501   2%   McDaniel

f113.h92   CFC-113    786 -  990   0.500     408   8%   McDaniel (omitted 203K)

sf6.h92    SF6        925 -  955   0.010    3001   2%   Varanasi

f142b.h92  HCFC-142b  870 - 1270   0.010   40000   4%   Newnham

clno3.h92  ClNO3      690 -  880   0.01    19000   2%   Birk & Wagner
           ClNO3      965 - 1005   0.01     4000   4%   Birk & Wagner
           ClNO3     1090 - 1130   0.01     4000   4%   Birk & Wagner
           ClNO3     1215 - 1330   0.01    11500   3%   Birk & Wagner
           ClNO3     1680 - 1790   0.07142  1540   ?    Ballard (only 2 spectra)

n2o5.h92   N2O5       547 -  610   0.160     373   4%   NCAR
           N2O5       709 -  775   0.210     315   2%   NCAR
           N2O5      1194 - 1281   0.350     266   2%   NCAR
           N2O5      1663 - 1793   0.480     271   2%   NCAR

c2h6.101   C2H6      1350 - 1496   0.0025  58401   1%   Sharpe (PNNL), Brown (JPL) 	

ch3cn.101  CH3CN      870 - 1650   0.05    15601   1%   Rinsland et al.(2005) 	

pan.101    PAN        760 - 870    0.1 	    1101   1%   Allen et al.(2005) 	
                     1110 - 1340   0.1 	    2301   1%   Allen et al.(2005)
                     1700 - 1780   0.1 	     801   1%   Allen et al.(2005) 	

ch3cho.101 CH3CHO    1000 - 1900   0.05    18001   2%   Sharpe (PNNL) 	
                     2600 - 2900   0.005   60001   6%   Sharpe (PNNL), Wennberg (CIT)

============================================================================


COLLISION INDUCED ABSORPTION LINELISTS:

Files: fcia_20060420.101, scia_20060420.101

GAS  Interval     Spacing  Lines  Measurer
============================================================================
O2   1275 -  1905  1.0      631   Thibault et al.(1997)
N2   2030 -  2705  1.0      676   Lafferty et al.(1996)
N2   4330 -  4930  1.0      601   Shapiro and Gush (1966)
O2   7100 - 10100  1.0     3001   Smith and Newnham (2000)
O2  12751 - 13750  1.0     1000   Tran et al.
============================================================================

Please direct any questions or comments to:
       Geoff Toon        818 354 8259    Geoffrey.C.Toon@jpl.nasa.gov
or     Armin Kleinboehl  818 393 6421    Armin.Kleinboehl@jpl.nasa.gov
