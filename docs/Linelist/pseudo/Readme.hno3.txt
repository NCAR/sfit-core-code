Supplemental HNO3 Pseudo-linelist

This HNO3 linelist is intended to supplement HITRAN,
which extends only to 1770 cm-1 and therefore misses
the nu1 fundamental, plus dozens of overtone and
combinations bands.

This supplemental linelist covers:
2600-2670 cm-1  2nu3 
3506-3595 cm-1  nu1
4004-4007 cm-1  nu1+nu9 (Q-branch only)
These are the three most noticeable missing HNO3 bands in MkIV balloon spectra.

The linelist was derived by fitting MkIV balloon spectra and
PNNL lab spectra. This is described with figures in:
http://mark4sun/report/HNO3_Spectroscopy_Evaluation.pdf

Description.
 Generated fake HNO3 linelists to cover the nu1 band centered
 at 3551 cm-1 and the 2nu3 band centered at 2645 cm-1, both of
 which are missing from HITRAN.  I did this by
 exploiting the similarity of the nu1 and nu2 bands. I first
 generated the following list of the positions of the manifolds
 in the nu1, nu2, and nu3 bands

 j   f_nu1   f_nu2   f_2nu3
-51 3529.17  1686.26    0.00
-50 3529.61  1686.75    0.00
-49 3530.03  1687.22    0.00
-48 3530.45  1687.71    0.00
-47 3530.89  1688.20    0.00
-46 3531.31  1688.68    0.00
-45 3531.74  1689.15    0.00
-44 3532.16  1689.63    0.00
-43 3532.58  1690.12    0.00
-42 3533.01  1690.60    0.00
-41 3533.43  1691.06    0.00
-40 3533.86  1691.55    0.00
-39 3534.30  1692.02    0.00
-38 3534.72  1692.49    0.00
-37 3535.14  1692.97    0.00
-36 3535.57  1693.43    0.00
-35 3535.99  1693.90    0.00
-34 3536.42  1694.36    0.00
-33 3536.84  1694.83    0.00
-32 3537.26  1695.29    0.00
-31 3537.69  1695.76    0.00
-30 3538.11  1696.21    0.00
-29 3538.54  1696.66    0.00
-28 3538.97  1697.12    0.00
-27 3539.39  1697.59    0.00
-26 3539.81  1698.03    0.00
-25 3540.225 1698.48    0.00
-24 3540.65  1698.95 2634.52
-23 3541.07  1699.39 2634.99
-22 3541.50  1699.84 2635.46
-21 3541.92  1700.29 2635.89
-20 3542.34  1700.73 2636.35
-19 3542.76  1701.18 2636.80
-18 3543.17  1701.62 2637.255
-17 3543.58  1702.06 2637.70
-16 3544.00  1702.51 2638.155
-15 3544.40  1702.94 2638.60
-14 3544.82  1703.38 2639.04
-13 3545.26  1703.81 2639.495
-12 3545.68  1704.24 2639.93
-11 3546.11  1704.69 2640.375
-10 3546.52  1705.11 2640.81
 -9 3546.94  1705.54 2641.24
 -8 3547.38  1705.97 2641.675
 -7 3547.80  1706.39 2642.11
 -6 3548.22  1706.82 2642.53
 -5 3548.65  1707.25 2642.96
 -4 3549.07  1707.68 2643.38
 -3 3549.49  1708.10 2643.80
 -2 3549.91  1708.52 2644.22
 -1 3550.33  1708.96 2644.65
  0 3550.74  1709.40 2645.07
  1 3551.16  1709.82 2645.49
  2 3551.58  1710.23 2645.90
  3 3551.99  1710.64 2646.31
  4 3552.41  1711.04 2646.71
  5 3552.83  1711.45 2647.12
  6 3553.24  1711.85 2647.53
  7 3553.66  1712.25 2647.93
  8 3554.08  1712.65 2648.33
  9 3554.50  1713.06 2648.725
 10 3554.91  1713.46 2649.125
 11 3555.32  1713.87 2649.52
 12 3555.73  1714.28 2649.91
 13 3556.14  1714.67 2650.30
 14 3556.55  1715.07 2650.69
 15 3556.96  1715.47 2651.075
 16 3557.39  1715.87 2651.46
 17 3557.79  1716.26 2651.85
 18 3558.24  1716.65 2652.23
 19 3558.64  1717.04 2652.61
 20 3559.04  1717.43 2653.01
 21 3559.44  1717.82 2653.38
 22 3559.89  1718.21 2653.77
 23 3560.29  1718.60 2654.175
 24 3560.72  1718.98    0.00
 25 3561.13  1719.37    0.00
 26 3561.55  1719.75    0.00
 27 3561.98  1720.12    0.00
 28 3562.38  1720.51    0.00
 29 3562.78  1720.89    0.00
 30 3563.20  1721.27    0.00
 31 3563.61  1721.63    0.00
 32 3564.02  1722.01    0.00
 33 3564.43  1722.39    0.00
 34 3564.84  1722.76    0.00
 35 3565.25  1723.12    0.00
 36 3565.66  1723.49    0.00
 37 3566.07  1723.86    0.00
 38 3566.48  1724.23    0.00
 39 3566.89  1724.58    0.00
 40 3567.31  1724.94    0.00
 41 3567.72  1725.30    0.00
 42 3568.12  1725.65    0.00
 43 3568.53  1726.01    0.00
 44 3568.94  1726.37    0.00
 45 3569.35  1726.72    0.00
 46 3569.76  1727.07    0.00
 47 3570.17  1727.41    0.00
 48 3570.57  1727.77    0.00
 49 3570.97  1728.11    0.00
 50 3571.37  1728.46    0.00
 51 3571.78  1728.81    0.00
 52 3572.18  1729.16    0.00
 53 3572.58  1729.50    0.00
 54 3572.99  1729.84    0.00

I then wrote a fortran program (linelist/f77/hno3.f) that fitted a
polynomial to the nu1 and 2nu3 frequencies as functions of the nu2
manifold frequencies, .e.g.

 f_nu1-3550.7 = a0+a1*(f_nu2-1709.5)+a2*(f_nu2-1709.5)^2+a3*(f_nu2-1709.5)^3

It turned out that a 4'th order polynomial produced a fit with a
rms deviation of 0.014 cm-1. Higher order polynomials failed to
improve significantly upon this.

The coefficients of the chosen polynomial were
a0 = 0.162455291
a1 = 1.00235796
a2 = 0.00362555636
a3 = 2.84212583E-05

Note that the a1 coefficient is very close to 1.0, indicating
that the average spacing of the manifolds is very similar in
the nu1 and nu2 bands.

I then used the polynomial to transform the frequencies of
the nu2 HNO3 linelist into a nu1 linelist. At the same time,
I scaled the strengths by a factor 0.20.  I then had to do
some massaging because the nu1 manifolds were wider than their
nu2 counterparts and because the nu1 Q-branch is relatively
stronger and has a different shape. This involved moving
individual HNO3 lines which were producing bumps in the
residuals to nearby frequencies where there were dips in
the residuals.

Also kludged a HNO3 pseudo-linelist for the nu1+nu9 Q-branch
at 4006.5 cm-1. 
