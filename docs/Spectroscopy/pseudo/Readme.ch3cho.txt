
              CH3CHO PSEUDO-LINELIST

           A. Kleinboehl and G. C. Toon


INTRODUCTION.
This document gives information on the pseudo-linelist for
Acetaldehyde (CH3CHO) derived at JPL in November 2005.
The linelist was created based on three data sets of absorption
cross-sections from laboratory measurements taken at the
Pacific Northwest National Laboratory (PNNL) by Steven Sharpe
and co-workers, and two absorption spectra taken at the California
Institute of Technology (CIT) by Paul Wennberg and co-workers.
The PNNL data consists of pressure broadened measurements at
three different temperatures and moderate resolution. The CIT
data is available only at room temperature but with high
resolution in both a pressure broadened and a non-pressure
broadened case. The characteristic quantities for the different
measurement sets are given in the following table.

 #  Source Temp. p_ch3cho p_tot Freq.-range  Res.  Spacing l_cell 
-----------------------------------------------------------------
 1   PNNL   323     2.2    760    510-6500  0.1125  0.060   19.96
 2   PNNL   298     2.2    760    510-6500  0.1125  0.060   19.96
 3   PNNL   278     2.2    760    510-6500  0.1125  0.060   19.96
 4   CIT    293    10.7    10.7  2500-3500  0.008   0.004   10.30
 5   CIT    293    10.7    513   2500-3500  0.008   0.004   10.30
-----------------------------------------------------------------
Temp.       - Temperature in K
p_ch3cho    - CH3CHO partial pressure in torr (assumed for PNNL)
P_tot       - Total pressure in torr
Freq.-range - Frequency range in cm-1
Res.        - Resolution in cm-1
Spacing     - Spectral point spacing in cm-1
l_cell      - Cell length in cm


DESCRIPTION.
The PNNL cross-sections were converted back into transmittance
spectra using the cell length and an average gas concentration.
The CIT spectra were ratioed using an empty cell spectrum taken in
the same set-up as the absorption spectra. All resulting
laboratory transmittance spectra were then simultaneously fitted
(using the GFIT algorithm) in appropriate frequency bands by
iteratively adjusting the strengths and ground-state energies
of the pseudo-lines.

Due to the different frequency coverage and resolution of the PNNL and
CIT spectra, different pseudo-line spacings were used for fitting different
frequency bands. In a region between 1000 and 1900 cm-1, where only PNNL
spectra were available, a pseudo-line spacing of 0.05 cm-1 was used.
In a region between 2600 and 2900 cm-1, where both PNNL and CIT spectra
were available, a finer pseudo-line spacing of 0.005 cm-1 was used.
Fitting was performed in the regions 1000-1250 cm-1, 1250-1650 cm-1,
1650-1900 cm-1, and 2600-2900 cm-1. Because there were small residuals of
water vapor and methane in the CIT spectra the water vapor and methane amounts
were fitted in a narrow frequency window and taken into account during the
fitting of the region at 2600-2900 cm-1. At each line frequency, an effective
strength and ground-state energy was derived by simultaneous non-linear least
squares fitting to the spectra that covered the frequency band of interest.
In analogy to other organics (e. g. CH3OH) an airbroadened halfwidth of
0.1 cm-1/atm and a selfbroadened halfwidth of 0.4 cm-1/atm was assumed.
In the regions between 1000-1900 cm-1 both a continuum level and a continuum
tilt were fitted simultaneously to the spectra, in the 2600-2900 region only
a continuum level was adjusted.

The result of the fitting process is a pseudo-linelist that covers the
frequency region 1000-1900 cm-1 and the region 2600-2900 cm-1, containing
78002 pseudo-lines in total. It has to be noted that at 2900 cm-1 the CH3CHO
absorption is still above zero, which may result in a discontinuity if
atmospheric spectra in a frequency window that reaches beyond 2900 cm-1 are
fitted with the present linelist.


PARTITION FUNCTION.
The rotational partition function for acetaldehyde was assumed to be (296/T)^1.5.
The vibrational partition function was calculated in the way it had been
done for the ATMOS experiment, as described e. g. by Norton and Rinsland (1991).
The following 15 vibrational frequencies (in cm-1) were used,
all degeneracies were set to 1:

3005, 2917, 2822, 1743, 1441, 1400, 1352, 1113, 919, 509,
2967, 1420, 867, 763, 150


ACCURACY.
To estimate how well the pseudo-linelist represents the laboratory spectra,
test retrievals were performed in which the laboratory spectra were fitted
using the pseudo-linelist. The retrieved scale factors for the acetaldehyde
abundances in the different spectra are tabulated below.

                    Scale factors retrieved in freq. region
 #     1000-1250 cm-1  1250-1650 cm-1  1650-1900 cm-1  2600-2900 cm-1
----------------------------------------------------------------------
 1        0.9982          1.0026          1.0062          1.0360
 2        1.0014          0.9945          0.9886          1.0388                 
 3        0.9994          1.0022          1.0053          1.0588
 4                                                        0.9924
 5                                                        0.9826
-----------------------------------------------------------------------

In the 1000-1900 cm-1 region the pseudolines correctly represent the PNNL
spectra to better than 1.5% of the given CH3CHO amount. In the region from
2600-2900 cm-1, where both PNNL and CIT spectra are available, the spectra
can be fitted to better than 6%, indicating a slight inconsistency between
the PNNL and the CIT measurements.


ACKNOWLEDGMENTS.
We would like to thank S. W. Sharpe, T. J. Johnson, and R. L. Sams
for the acetaldehyde cross-sections from PNNL, as well as P. O. Wennberg and
R. A. Washenfelder for the acetaldehyde spectra from CIT.


REFERENCE.

Norton, R. H. and C. P. Rinsland, ATMOS data processing and science
analysis methods, Appl. Opt., 30, 389-400, 1991. 
