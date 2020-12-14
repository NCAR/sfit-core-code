# v2 = Following Launch of the C2H6 global study (see doc file Launch of
#      the C2H6 global study_v2.docx and LAUNCH folder)

# Filenames Section

 file.in.stalayers                   = ./station.layers
 file.in.sbdflt                       = /idon.t.kno.junk/sbdefault.ctl
 file.in.refprofile                  = ./reference.prf
 file.in.spectrum                    = ./t15asc.4
 file.in.modulation_fcn              = ./ils_poly2.mf
 file.in.phase_fcn                   = ./ils_poly2.phs
 file.in.sa_matrix                   =
 file.in.isotope                     =
 file.in.solarlines                  = ../solar.dat
 file.in.linelist                    = ./02972.601722-02990.908278.hbin
# file.in.linelist                    = ./02972.601722-02990.908278-2008.hbin
# file.in.linelist                    = ./02972.601722-02990.908278.hbin


 file.out.pbpfile                    = pbpfile
 file.out.statevec                   = statevec
 file.out.k_matrix                   = k.output
 file.out.kb_matrix                  = kb.output
 file.out.sa_matrix                  = sa.complete
 file.out.retprofiles                = rprfs.table
 file.out.aprprofiles                = aprfs.table
 file.out.ab_matrix                  = ab.output
 file.out.summary                    = summary
 file.out.parm_vectors               = parm.output
 file.out.seinv_vector               = seinv.output
 file.out.sainv_matrix               = sainv.complete
 file.out.smeas_matrix               = smeas.target
 file.out.shat_matrix                = shat.complete
 file.out.ssmooth_matrix             = ssmooth.target
 file.out.ak_matrix                  = ak.target
 file.out.solarspectrum              = solar.output
 file.out.g_matrix		             = d.complete

# Definition for retrieval gases

 gas.layers                              = 44

# Profiles
 gas.profile.list                        = C2H6 O3
 gas.profile.C2H6.correlation            = T
 gas.profile.C2H6.correlation.type       = 2
 gas.profile.C2H6.correlation.width      = 4.0
 gas.profile.C2H6.correlation.minalt     = 0.0
 gas.profile.C2H6.correlation.maxalt     = 120.0
 gas.profile.C2H6.logstate               = F

 gas.profile.C2H6.scale                  = 2.7
 gas.profile.C2H6.sigma                  =
0.2673 0.2949 0.3279 0.3492 0.3719 0.3962 0.4218
0.4504 0.4822 0.5171 0.5556 0.5976 0.6415 0.6884
0.7332 0.7715 0.8032 0.8468 0.8904 0.9588 1.0754
1.2669 1.0427 0.6787 0.4584 0.2974 0.2512 0.2235
0.2027 0.2028 0.2272 0.2246 0.2403 0.2592 0.2847
0.3103 0.3287 0.3387 0.3418 0.3392 0.3376 0.3469
0.3413 0.3242


gas.profile.O3.correlation = T
gas.profile.O3.correlation.type   = 2
gas.profile.O3.correlation.width  = 10.0
gas.profile.O3.correlation.minalt = 0.0
gas.profile.O3.correlation.maxalt = 120.0
gas.profile.O3.logstate    = F
gas.profile.O3.scale       = 1.0
gas.profile.O3.sigma       =
0.0792 0.0844  0.0901  0.0964  0.1034  0.1111  0.1195
0.1283  0.1377  0.1466  0.1543  0.1606  0.1644 0.1667
0.1684  0.1709  0.1728  0.1754  0.1782 0.1803  0.1833
0.1857  0.1890  0.1925  0.1952 0.1990  0.2020  0.2064
0.2111  0.2150  0.2203 0.2247  0.2309  0.2375  0.2433
0.2510  0.2576 0.2668  0.2776  0.2863  0.2995  0.3147
0.3464  0.3563

# actually use column scale for H2O (below)
gas.profile.H2O.correlation         = T
gas.profile.H2O.correlation.type    = 2
gas.profile.H2O.correlation.width   = 4.0
gas.profile.H2O.correlation.minalt  = 0.0
gas.profile.H2O.correlation.maxalt  = 120.0
gas.profile.H2O.logstate            = F
gas.profile.H2O.scale               = 1
gas.profile.H2O.sigma               =
0.080178    0.088465    0.098374    0.104765    0.111571    0.118864
0.126547    0.135113    0.144673    0.155126    0.166667    0.179284
0.192450    0.206529    0.219971    0.231455    0.240966    0.246598
0.250000    0.252646    0.256307    0.259161    0.263117    0.267261
0.270501    0.275010    0.278543    0.283473    0.288675    0.292770
0.298511    0.303046    0.309295    0.315789    0.320879    0.328148
0.333870    0.342104    0.350955    0.357981    0.368161    0.376296
0.388173    0.401286

# Columns
 gas.column.list                     = CH4 CH3CL H2O

 gas.column.H2O.scale                = 1.0
 gas.column.H2O.sigma                = 1.0
 gas.column.O3.scale                 = 1.0
 gas.column.O3.sigma                 = 0.1
 gas.column.CH4.scale                = 1.0
 gas.column.CH4.sigma                = 0.4
 gas.column.CH3CL.scale              = 1.0
 gas.column.CH3CL.sigma              = 0.1


# Forward model parameters

 fw.raytonly                         = F
 fw.isotope_separation               = F
 fw.delnu                            = 0.1000
 fw.lshapemodel                      = 1
 fw.linemixing                       = F
 fw.linemixing.gas                   = CO2
 fw.tips                             = F
 fw.solar_spectrum                   = T
 fw.pressure_shift                   = T
 fw.apod_fcn                         = F
 fw.apod_fcn.type                    = 2
 fw.apod_fcn.order                   = 2
 fw.phase_fcn                        = F
 fw.phase_fcn.type                   = 2
 fw.phase_fcn.order                  = 2
 fw.emission                         = F


# Retrieval parameter

 rt                                  = T
 rt.lm                               = F
 rt.lm.gamma_start                   = 1.0e5
 rt.lm.gamma_inc                     = 10.0
 rt.lm.gamma_dec                     = 10.0
 rt.convergence                      = 0.01
 rt.tolerance                        = 0.05
 rt.max_iteration                    = 25
 rt.ifcalc_se                        = T
 rt.dwshift                          = F
 rt.wshift                           = T
 rt.wshift.type                      = 2
 rt.wshift.apriori                   = 0.000
 rt.wshift.sigma                     = 0.100
 rt.slope                            = T
 rt.slope.apriori                    = 0.000
 rt.slope.sigma                      = 0.100
 rt.curvature                        = F
 rt.curvature.apriori                = 0.000
 rt.curvature.sigma                  = 0.100
 rt.phase                            = F
 rt.phase.apriori                    = 0.000
 rt.phase.sigma                      = 0.200
 rt.apod_fcn                         = F
 rt.apod_fcn.apriori                 = 1.000
 rt.apod_fcn.sigma                   = 0.050
 rt.phase_fcn                        = F
 rt.phase_fcn.apriori                = 1.000
 rt.phase_fcn.sigma                  = 0.200
 rt.solshift                         = T
 rt.solshift.apriori                 = 0.000
 rt.solshift.sigma                   = 0.010
 rt.solstrnth                        = F
 rt.solstrnth.apriori                = 0.000
 rt.solstrnth.sigma                  = 0.010
 rt.temperature                      = F


# Kb derivative calculations

 kb                                  = F
 kb.temperature                      = T
 kb.slope                            = F
 kb.curvature                        = T
 kb.solshft                          = F
 kb.solstrnth                        = F
 kb.phase                            = F
 kb.wshift                           = F
 kb.apod_fcn                         = T
 kb.phase_fcn                        = T
 kb.zshift                           = F
 kb.sza                              = T
 kb.omega                            = T
 kb.line                             = T
 kb.line.type                        = 1
 kb.line.gas                         = retrieval


# Microwindows and their parameters

 band                                = 1 2 3

 band.1.nu_start                     = 2976.66
 band.1.nu_stop                      = 2977.059
 band.1.zshift                       = F
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.200
 band.1.calc_point_space             = 0.0005
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 180.00
 band.1.omega                        = 2.3923
 band.1.apodization_code             = 0
 band.1.gasb                         = C2H6 H2O O3 CH4 CH3CL
 band.1.tempretb                     = F

 band.2.nu_start                     = 2983.2000
 band.2.nu_stop                      = 2983.50
 band.2.zshift                       = F
 band.2.zshift.type                  = 1
 band.2.zshift.apriori               = 0.000
 band.2.zshift.sigma                 = 0.200
 band.2.calc_point_space             = 0.0005
 band.2.wave_factor                  = 1.000
 band.2.max_opd                      = 180.00
 band.2.omega                        = 2.3923
 band.2.apodization_code             = 0
 band.2.gasb                         = C2H6 H2O O3 CH4 CH3CL
 band.2.tempretb                     = F

 band.3.nu_start                     = 2986.4500
 band.3.nu_stop                      = 2986.850
 band.3.zshift                       = F
 band.3.zshift.type                  = 1
 band.3.zshift.apriori               = 0.000
 band.3.zshift.sigma                 = 0.200
 band.3.calc_point_space             = 0.0005
 band.3.wave_factor                  = 1.000
 band.3.max_opd                      = 180.00
 band.3.omega                        = 2.3923
 band.3.apodization_code             = 0
 band.3.gasb                         = C2H6 H2O O3 CH4 CH3CL
 band.3.tempretb                     = F


# Spectrum Section

 sp.snr =
 sp.snr.1.nu_start = 2976.7
 sp.snr.1.nu_stop  = 2976.9
 sp.snr.1.snr  = 800.0


# Output Files Section

 out.level                           = 1
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.g_matrix			                = T
 out.k_matrix                        = T
 out.shat_matrix                     = T
 out.retprofiles                     = T
 out.aprprofiles                     = T
 out.ab_matrix                       = F
 out.ak_matrix                       = T
 out.summary                         = T
 out.pbpfile                         = T
 out.channel                         = F
 out.parm_vectors                    = T
 out.seinv_vector                    = T
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = T
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F
