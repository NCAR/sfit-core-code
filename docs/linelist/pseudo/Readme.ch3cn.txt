
            CH3CN PSEUDO-LINELIST

        A. Kleinboehl and G. C. Toon


INTRODUCTION.
This document gives information on the CH3CN pseudo-linelist
derived at JPL in May 2005. The linelist was created based on
on 29 laboratory spectra taken at the Pacific Northwest National
Laboratory (PNNL) and provided by Steven Sharp. These are not
the same as the CH3CN cross-sections on the HITRAN website, which
have been mutilated by converting all negative values to zero.
The measurements and the absorption cross sections
(incl. assignments of major bands) are described by
Rinsland et al. (2005).

The measurement conditions for each of these spectra are tabulated below.
Each measurement used the same cell of 8.1576 m length. Each spectrum
covers a region between 600 and 6500 cm-1 with a resolution of 0.1125 cm-1
and a spectral point spacing of 0.0603 cm-1. 

 #  File          Temp   P_tot  P_ch3cn
----------------------------------------
 1  "CH3CNA.D01"  298.7  750.5  0.137679
 2  "CH3CNA.D02"  298.7  749.4  0.068739
 3  "CH3CNA.D03"  298.7  749.0  0.206106
 4  "CH3CNA.D04"  298.7  747.8  0.411552
 5  "CH3CNA.D05"  298.7  747.8  0.274368
 6  "CH3CNA.D06"  298.7  747.9  0.686012
 7  "CH3CNA.D07"  298.7  748.1  0.548956
 8  "CH3CNA.D08"  298.7  748.1  1.029293
 9  "CH3CNA.D09"  298.7  748.2  1.235317
10  "CH3CNA.D10"  298.7  748.3  0.864837
11  "CH3CNA.D11"  298.7  748.6  0.343327
12  "CH3CNB.D01"  276.1  752.7  0.069041
13  "CH3CNB.D02"  276.1  752.7  0.096658
14  "CH3CNB.D03"  276.1  752.4  0.138028
15  "CH3CNB.D04"  276.1  752.3  0.207014
16  "CH3CNB.D05"  276.1  752.3  0.276019
17  "CH3CNB.D06"  276.1  752.2  0.413974
18  "CH3CNB.D07"  276.1  752.2  0.551965
19  "CH3CNB.D08"  276.1  752.0  0.689773
20  "CH3CNB.D09"  276.1  751.9  1.034521
21  "CH3CNC.D01"  324.1  750.0  0.068794
22  "CH3CNC.D02"  324.1  750.0  0.137588
23  "CH3CNC.D03"  324.1  749.9  0.275139
24  "CH3CNC.D04"  324.1  749.9  0.687846
25  "CH3CNC.D05"  324.1  750.1  0.206409
26  "CH3CNC.D06"  324.1  750.1  0.412818
27  "CH3CNC.D07"  324.1  750.2  0.550497
28  "CH3CNC.D08"  324.1  750.4  0.344152
29  "CH3CNC.D09"  324.1  750.6  1.032733

Temp - Temperature in K
P_tot - total pressure in torr
P_ch3cn - CH3CN partial pressure in torr


DESCRIPTION.
First, the cross-sections were converted back into transmittance
spectra from knowledge of the cell length and gas concentrations.
The resulting laboratory transmittance spectra were then
simultaneously fitted (using the GFIT algorithm) by iteratively
adjusting the strengths and ground-state energies of the pseudo-lines.

Due to the resolution of the laboratory spectra of 0.1125 cm-1
a pseudo-line spacing of 0.05 cm-1 was considered to be appropriate.
Fitting was performed in the frequency regions around 900 cm-1
(where the nu_4 band is located), around 1050 cm-1 (where the
nu_7 band is located), and around 1450 cm-1 (for the nu_3, nu_6, and
nu_7+nu_8 bands). These regions include the two bands with the strongest
absorption features. A zero level offset of 0.2% has been assumed throughout,
based on fits to spectra in which the absorption feature at 1463 cm-1
was saturated. The result of the fitting process is a continuous
pseudo-linelist containing 15601 lines between 870 and 1650 cm-1.


CALCULATION OF S, E", and ABHW.
At each line frequency, an effective strength and ground-state energy
was derived by simultaneous non-linear least squares fitting to the 29
spectra. Furthermore, the ABHW was calculated from the ground-state
energy using the formula ABHW = 0.04 * (E" + 2000)/(E" + 1000),
giving a ABHW of 0.08 cm-1/atm for E"->0 and a ABHW of 0.04 cm-1/atm
for E"->oo. This formulation seemed to be the most approriate to fit
the feature at 1042 cm-1, which is the narrowest feature in the
considered frequency region. These widths are smaller than those measured
by Drouin (2003) in the microwave region, but we found that using
larger widths produced significantly poorer fits to the sharp spectral
features. As part of the fitting, the strengths and ground-state energies
were both constrained to be positive.


PARTITION FUNCTION.
The rotational partition function for CH3CN was assumed to be (296/T)^1.5.
The vibrational partition function was calculated in the way it had been
done for the ATMOS experiment, as described e. g. by Norton and Rinsland (1991).
The following vibrational frequencies and degeneracies were assumed:

freq. | 2954  2267  1385   920  3009  1448  1041   362
deg.  |    1     1     1     1     2     2     2     2


ACCURACY.
To estimate how well the pseudo-linelist represents the PNNL spectra,
test retrievals were performed in which the laboratory spectra were fitted
using the pseudo-linelist. The retrieved scale factors for the CH3CN
abundances in the different spectra are tabulated below.

         Scale factors retrieved in freq. region
#	900 cm-1	1050 cm-1	1450 cm-1
-------------------------------------------------
1	1.0073		0.9984		1.0065
2	1.0180		0.9993		1.0071
3	1.0043		0.9938		1.0013
4	1.0089		0.9964		1.0032
5	1.0073		0.9973		1.0037
6	1.0058		0.9953		0.9998
7	1.0018		0.9925		0.9975
8	1.0013		0.9928		0.9949
9	1.0001		0.9930		0.9946
10	0.9968		0.9880		0.9906
11	1.0033		0.9905		0.9969
12	0.9847		0.9733		0.9820
13	0.9993		0.9935		1.0019
14	1.0013		1.0003		1.0079
15	0.9885		0.9830		0.9908
16	0.9985		0.9963		1.0032
17	1.0040		1.0015		1.0083
18	1.0034		1.0004		1.0068
19	0.9930		0.9908		0.9972
20	1.0016		1.0006		1.0064
21	0.9968		0.9974		1.0086
22	1.0020		0.9944		1.0052
23	0.9987		0.9928		1.0045
24	0.9999		0.9952		1.0048
25	0.9987		0.9934		1.0054
26	0.9918		0.9891		1.000
27	0.9993		0.9928		1.0035
28	1.0087		1.0009		1.0119
29	0.9958		0.9940		1.0013
-----------------------------------------------
mean	1.00072		0.99403		1.00158
stddev	0.00655		0.00588		0.00655

The pseudolines correctly represent the PNNL spectra
to within 0.7% of the given CH3CN amount in all bands
analyzed. The main exception to this is spectrum 12
which appears to contain ~2% less CH3CN than advertised.


ACKNOWLEDGMENTS.
We would like to thank Steven Sharpe and Curtis Rinsland for
providing CH3CN cross-sections prior to publication, and
Aaron Goldman and Isabelle Kleiner for valuable discussions
about the CH3CN partition function.


REFERENCES.

Drouin, B. J., Temperature dependent rotational transition lineshape
parameters for O3, O2, SO2, CH3CN and CO, International Symposium on
Molecular Spectroscopy, Columbus, OH, USA, June 16-20, 2003.

Norton, R. H. and C. P. Rinsland, ATMOS data processing and science
analysis methods, Appl. Opt., 30, 389-400, 1991. 

Rinsland, C. P., S. W. Sharpe, and R. L. Sams, Temperature-dependent
infrared absorption cross-sections of methyl cyanide (acetonitrile),
JQSRT, 96, 271-280, 2005.
