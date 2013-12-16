# Filenames Section

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = t15asc.4
# file.in.spectrum                    =  t15asc.4.1.2529
 file.in.modulation_fcn              = ils.dat
 file.in.phase_fcn                   = ils.dat
 file.in.sa_matrix                   = sa.input
 file.in.isotope                     = isotope.input
 file.in.solarlines                  = ../../linelist/solar/120621/solar.dat
 file.in.linelist                    = ./02723.588856-02930.241144.hbin

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

# Definition for retrieval gases

 gas.layers                          = 47

 gas.profile.list                    = HCL O3 CH4
 # O3

 gas.profile.HCL.scale               = 1.0D0
 gas.profile.HCL.correlation         = T
 gas.profile.HCL.correlation.type    = 2
 gas.profile.HCL.correlation.width   = 6.0D0
 gas.profile.HCL.correlation.minalt  = 0.0D0
 gas.profile.HCL.correlation.maxalt  = 120.0D0
 gas.profile.HCL.logstate            = F
 gas.profile.HCL.sigma               =
    0.026726    0.029488    0.032791    0.034922    0.037190    0.039621
    0.042182    0.045038    0.048224    0.051709    0.055556    0.059761
    0.064150    0.068843    0.073324    0.077152    0.080322    0.082199
    0.083333    0.084215    0.085436    0.086387    0.087706    0.089087
    0.090167    0.091670    0.092848    0.094491    0.096225    0.097590
    0.099504    0.101015    0.103098    0.105263    0.106960    0.109383
    0.111290    0.114035    0.116985    0.119327    0.122720    0.125432
    0.129391    0.133762    0.137283    0.142538    0.148446

 gas.profile.O3.scale                = 1.0D0
 gas.profile.O3.correlation          = T
 gas.profile.O3.correlation.type     = 1
 gas.profile.O3.correlation.width    = 3.0D0
 gas.profile.O3.correlation.minalt   = 0.0D0
 gas.profile.O3.correlation.maxalt   = 120.0D0
 gas.profile.O3.logstate             = F
 gas.profile.O3.sigma                =
    0.026726    0.029488    0.032791    0.034922    0.037190    0.039621
    0.042182    0.045038    0.048224    0.051709    0.055556    0.059761
    0.064150    0.068843    0.073324    0.077152    0.080322    0.082199
    0.083333    0.084215    0.085436    0.086387    0.087706    0.089087
    0.090167    0.091670    0.092848    0.094491    0.096225    0.097590
    0.099504    0.101015    0.103098    0.105263    0.106960    0.109383
    0.111290    0.114035    0.116985    0.119327    0.122720    0.125432
    0.129391    0.133762    0.137283    0.142538    0.148446

 gas.profile.CH4.scale                = 1.0D0
 gas.profile.CH4.correlation          = T
 gas.profile.CH4.correlation.type     = 1
 gas.profile.CH4.correlation.width    = 3.0D0
 gas.profile.CH4.correlation.minalt   = 0.0D0
 gas.profile.CH4.correlation.maxalt   = 120.0D0
 gas.profile.CH4.logstate             = F
 gas.profile.CH4.sigma                =
    0.026726    0.029488    0.032791    0.034922    0.037190    0.039621
    0.042182    0.045038    0.048224    0.051709    0.055556    0.059761
    0.064150    0.068843    0.073324    0.077152    0.080322    0.082199
    0.083333    0.084215    0.085436    0.086387    0.087706    0.089087
    0.090167    0.091670    0.092848    0.094491    0.096225    0.097590
    0.099504    0.101015    0.103098    0.105263    0.106960    0.109383
    0.111290    0.114035    0.116985    0.119327    0.122720    0.125432
    0.129391    0.133762    0.137283    0.142538    0.148446

 gas.column.list                     = NO2 N2O HDO H2O

 gas.column.H2O.scale                = 1.0d0
 gas.column.H2O.sigma                = 0.1d0

 gas.column.O3.scale                 = 1.0d0
 gas.column.O3.sigma                 = 0.1d0

 gas.column.CH4.scale                = 1.0d0
 gas.column.CH4.sigma                = 0.1d0

 gas.column.HDO.scale                = 0.5d0
 gas.column.HDO.sigma                = 0.1d0

 gas.column.N2O.scale                = 1.0d0
 gas.column.N2O.sigma                = 0.01d0

 gas.column.NO2.scale                = 1.0d0
 gas.column.NO2.sigma                = 0.1d0


# Forward model parameters

 fw.raytonly                         = F
 fw.isotope_separation               = T
 fw.delnu                            = 1.0D0
 fw.lshapemodel                      = 0
 fw.linemixing                       = T
 fw.linemixing.gas                   = CO2
 fw.solar_spectrum                   = T
 fw.pressure_shift                   = T
 fw.apod_fcn                         = T
 fw.apod_fcn.type                    = 4
 fw.apod_fcn.order                   = 0
 fw.phase_fcn                        = T
 fw.phase_fcn.type                   = 4
 fw.phase_fcn.order                  = 0
 fw.emission                         = F
# fw.emission.T_infinity              = 6000.0
# fw.emission.object                  = E
# fw.emission.normalized              = T

# Retrieval parameter

 rt                                  = T
 rt.lm                               = T
 rt.lm.gamma_start                   = 10000.
 rt.lm.gamma_dec                     = 10.0
 rt.lm.gamma_inc                     = 10.0
 rt.convergence                      = 0.05
 rt.tolerance                        = 0.05
 rt.max_iteration                    = 21
 rt.wshift                           = T
 rt.wshift.type                      = 3
 rt.wshift.apriori                   = 0.000
 rt.wshift.sigma                     = 0.100
 rt.dwshift                          = T
 rt.slope                            = T
 rt.slope.apriori                    = 0.000
 rt.slope.sigma                      = 0.100
 rt.curvature                        = F
 rt.curvature.apriori                = 0.000
 rt.curvature.sigma                  = 0.100
 rt.phase                            = F
 rt.phase.apriori                    = 0.000
 rt.phase.sigma                      = 0.000
 rt.apod_fcn                         = F
 rt.apod_fcn.apriori                 = 0.000
 rt.apod_fcn.sigma                   = 0.000
 rt.phase_fcn                        = F
 rt.phase_fcn.apriori                = 0.000
 rt.phase_fcn.sigma                  = 0.000
 rt.solshift                         = T
 rt.solshift.apriori                 = 0.000
 rt.solshift.sigma                   = 0.100
 rt.ifcalc_se                        = F
 rt.temperature                      = F
 rt.temperature.sigma                =
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0


# Microwindows and their parameters

 band                                = 1 2 3

 band.1.nu_start                     = 2727.630
 band.1.nu_stop                      = 2727.900
 band.1.zshift                       = F
 band.1.zshift.type                  = 0
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.200
 band.1.beam                         = 0
 band.1.calc_point_space             = 0.800E-03
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 257.143
 band.1.omega                        = 2.273
 band.1.apodization_code             = 0
 band.1.snr                          = 200.000
 band.1.gasb                         = HCL O3 H2O HDO

 band.2.nu_start                     = 2775.500
 band.2.nu_stop                      = 2775.900
 band.2.zshift                       = F
 band.2.zshift.type                  = 0
 band.2.zshift.apriori               = 0.000
 band.2.zshift.sigma                 = 0.200
 band.2.beam                         = 0
 band.2.calc_point_space             = 0.800E-03
 band.2.wave_factor                  = 1.000
 band.2.max_opd                      = 257.143
 band.2.omega                        = 2.273
 band.2.apodization_code             = 0
 band.2.snr                          = 200.000
 band.2.gasb                         = HCL O3 N2O H2O

 band.3.nu_start                     = 2925.600
 band.3.nu_stop                      = 2926.200
 band.3.zshift                       = F
 band.3.zshift.type                  = 0
 band.3.zshift.apriori               = 0.000
 band.3.zshift.sigma                 = 0.200
 band.3.beam                         = 0
 band.3.calc_point_space             = 0.800E-03
 band.3.wave_factor                  = 1.000
 band.3.max_opd                      = 257.143
 band.3.omega                        = 2.273
 band.3.apodization_code             = 0
 band.3.snr                          = 200.000
 band.3.gasb                         = HCL CH4 NO2
 band.3.tempretb                     = F


# spectrum parameters

 sp.snr                              =

 sp.snr.1.nu_start                   = 2727.0
 sp.snr.1.nu_stop                    = 2975.0
 sp.snr.1.snr                        = 200.0

 sp.snr.2.nu_start                   = 2727.4
 sp.snr.2.nu_stop                    = 2727.7
 sp.snr.2.snr                        = 100.0


# Output Files Section

 out.level                           = 1
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.k_matrix                        = T
 out.shat_matrix                     = T
 out.retprofiles                     = T
 out.aprprofiles                     = T
 out.ab_matrix                       = F
 out.ak_matrix                       = F
 out.summary                         = T
 out.pbpfile                         = T
 out.channel                         = F
 out.parm_vectors                    = T
 out.seinv_vector                    = T
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = F
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F

# Kb derivative calculations

 kb                                  = F
 kb.temperature                      = F
 kb.slope                            = T
 kb.curvature                        = T
 kb.solshft                          = T
 kb.solstrnth                        = T
 kb.phase                            = T
 kb.dwshift                          = T
 kb.wshift                           = T
 kb.apod_fcn                         = T
 kb.phase_fcn                        = T
 kb.zshift                           = T
 kb.sza                              = T
 kb.omega                            = T
 kb.max_opd                          = T
 kb.line                             = T
 kb.line.type                        = 1
 kb.line.gas                         = HCL

