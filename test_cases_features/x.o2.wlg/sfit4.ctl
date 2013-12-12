# edited for version d0.9.2 / 20130502 input parser

# Filenames Section

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = new_spectrum
 file.in.modulation_fcn              = ils_poly.dat
 file.in.phase_fcn                   = ilp_poly.dat
 file.in.sa_matrix                   = sa_38.input
 file.in.isotope                     = isotope.input
 file.in.solarlines                  = /Users/jamesw/FDP/sfit/400/sfit-core/linelist/solar/120621/solar.dat
 file.in.linelist                    = ./07760.768891-08009.231109.hbin

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

 gas.layers                          = 48

 gas.profile.list                    = O2

 gas.profile.O2.scale               = 1.0D0
 gas.profile.O2.correlation         = T
 gas.profile.O2.correlation.type    = 1
 gas.profile.O2.correlation.width   = 1000.0D0
 gas.profile.O2.correlation.minalt  = 0.0D0
 gas.profile.O2.correlation.maxalt  = 120.0D0
 gas.profile.O2.logstate            = F
 gas.profile.O2.sigma               =
    0.006736 0.007016 0.007284 0.007476 0.007729 0.007972
    0.008149 0.008380 0.008548 0.008769 0.008986 0.009142
    0.009349 0.009500 0.009699 0.009899 0.010050 0.010247
    0.010392 0.010583 0.010770 0.010909 0.011091 0.011225
    0.011402 0.011576 0.011705 0.011874 0.012000 0.012166
    0.012450 0.012961 0.013638 0.014526 0.015588 0.016733
    0.018000 0.019339 0.020736 0.022204 0.023707 0.025239
    0.026889 0.028636 0.030496 0.033912 0.037417 0.037417




 gas.column.list                     =
 #H2O O2CIA

 gas.column.O2.scale                        = 1.00000
 gas.column.O2.sigma                        = 0.04560

 gas.column.O2CIA.sigma                     = 0.00001
 gas.column.O2CIA.scale                     = 1.00000

 gas.column.H2O.sigma                       = 1.00000
 gas.column.H2O.scale                       = 1.00000

# Forward model parameters

 fw.raytonly                         = F
 fw.isotope_separation               = F
 fw.delnu                            = 1.000
 fw.lshapemodel                      = 0
 fw.linemixing                       = T
 fw.linemixing.gas                   = CO2
 fw.solar_spectrum                   = T
 fw.pressure_shift                   = T
 fw.apod_fcn                         = T
 fw.apod_fcn.type                    = 2
 fw.apod_fcn.order                   = 4
 fw.phase_fcn                        = T
 fw.phase_fcn.type                   = 2
 fw.phase_fcn.order                  = 4
 fw.emission                         = F


# Retrieval parameter

 rt                                  = F
 rt.lm                               = F
 rt.convergence                      = 0.01
 rt.tolerance                        = 0.050
 rt.max_iteration                    = 13
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
 rt.phase.sigma                      = 0.100
 rt.apod_fcn                         = F
 rt.apod_fcn.apriori                 = 0.000
 rt.apod_fcn.sigma                   = 0.100
 rt.phase_fcn                        = F
 rt.phase_fcn.apriori                = 0.000
 rt.phase_fcn.sigma                  = 0.100
 rt.solshift                         = F
 rt.solshift.apriori                 = 0.000
 rt.solshift.sigma                   = 0.010
 rt.solstrnth                        = F
 rt.solstrnth.apriori                = 0.000
 rt.solstrnth.sigma                  = 0.500
 rt.temperature                      = F


# Microwindows and their parameters

 band                                = 1
 band.1.gasb                         = O2
 #H2O O2CIA
 band.1.nu_start                     = 7765.000
 band.1.nu_stop                      = 8005.000
 band.1.zshift                       = F
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sa                    = 0.200
 band.1.beam                         = 0
 band.1.calc_point_space             = 5.000E-03
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 45.01
 band.1.omega                        = 1.9
 band.1.apodization_code             = 0
 band.1.snr                          = 200.000

# Spectrum parameters

 sp.snr                              = 1
 sp.snr.1.nu_start                   = 7765.0
 sp.snr.1.nu_stop                    = 8005.30
 sp.snr.1.snr                        = 200.0


# Output Files Section

 out.level                           = 3
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.k_matrix                        = T
 out.shat_matrix                     = F
 out.retprofiles                     = T
 out.aprprofiles                     = T
 out.ab_matrix                       = F
 out.ak_matrix                       = T
 out.summary                         = T
 out.pbpfile                         = T
 out.channel                         = T
 out.parm_vectors                    = T
 out.seinv_vector                    = T
 out.sainv_matrix                    = T
 out.smeas_matrix                    = T
 out.ssmooth_matrix                  = T
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = T
 out.levmardet                       = T
 out.xscdetail                       = T

# Kb derivative calculations

 kb                                  = F
 kb.temperature                      = F
 kb.slope                            = T
 kb.curvature                        = T
 kb.solshft                          = T
 kb.solstrnth                        = F
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
 kb.line.gas                         = O2

