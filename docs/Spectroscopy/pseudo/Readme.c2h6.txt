            Ethane PSEUDO-LINELIST

There are two C2H6 pseudo linelists:
c2h6_2953_3018.101
c2h6_1350_1496.101

The former is a hodge-podge described in:
http://mark4sun.jpl.nasa.gov/report/c2h6_spectroscopy.pdf
Basically it consists of:
 - a quantum-mechanically-derived linelist for the PQ3-branch
   at 2976 cm-1 which was included in HITRAN2004
 - linelists kludged by Linda Brown 20 years ago for ATMOS
   covering the 9 strongest PQ-branches from 2973-3001 cm-1
   These were in HITRAN2000 but not in HITRAN2004
 - linelists kludged by Geoff Toon for the 2953-2973 and
   3001-3018 cm-1 regions containing weaker PQ-branches

The latter linelist is described below.

INTRODUCTION.
This document gives information on a C2H6 pseudo-linelist
derived at JPL in May 2005. The linelist was created based on
on 3 laboratory spectra taken at the Pacific Northwest National
Laboratory (PNNL), provided by Steven Sharp, and one high-resolution
spectrum recorded at Kitt Peak National Solar observatory by
Linda Brown.

The pseudo-linelist covers the 1350-1496 cm-1 region containing
the nu_6 band, which is missing from HITRAN_2004.

[HITRAN_2004 contains C2H6 lines for the nu_12 (822 cm-1), nu_1
(2954 cm-1), and nu_10 (2985) bands, but not the nu_6 (1379 cm-1)].

The measurement conditions for these spectra are tabulated below.

File            Temp   P_tot  P_c2h6  l_cell spacing  res.
-----------------------------------------------------------
"c2h6pnnl_50C"  323.2  760.5   2.14    0.20  0.0603  0.1125
"c2h6pnnl_25C"  298.2  760.5   2.14    0.20  0.0603  0.1125
"c2h6pnnl_05C"  278.2  760.5   2.14    0.20  0.0603  0.1125
"r850530R0.017" 299.0    5.0   5.00    0.25  0.0029  0.0060
-----------------------------------------------------------
Temp - Temperature in K
P_tot - total pressure in torr
P_ch3cn - CH3CN partial pressure in torr
l_cell - cell length in m
spacing - spectral point spacing in cm-1
res. - spectral resolution in cm-1


DESCRIPTION.
The laboratory transmittance spectra were simultaneously
fitted (using the GFIT algorithm) by iteratively adjusting
the strengths and ground-state energies of the pseudo-lines.

Due to the high resolution of the Kitt Peak spectrum
a pseudo-line spacing of 0.0025 cm-1 was chosen.
An air broadened halfwidth of 0.068 cm-1/atm and a
self broadened halfwidth of 0.300 cm-1/atm
was assumed for all pseudolines.

The current version of the pseudo-linelist covers a frequency
range between 1350 and 1496 cm-1. This is only part of the
ethane band around 1500 cm-1. This region was chosen as
the Kitt Peak spectrum was contaminated with a significant
amount of water, blacking out several spectral regions
beyond 1496 cm-1.


ACKNOWLEDGMENTS.
We would like to thank Steven Sharpe and Linda Brown for the
C2H6 cross-sections and transmittance spectra, respectively.
