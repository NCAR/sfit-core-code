 file.in.spectrum                             = spectrum
 file.in.stalayers                            = station.layers
 file.in.refprofile                           = reference.prf
 file.in.solarlines                           = solar.dat
 file.in.linelist= 02759.361722-02786.058278.hbin
 file.in.isotope                             = isotope.input.h2o
 file.out.pbpfile= pbpfile
 file.out.statevec= statevec
 file.out.k_matrix= k.out
 file.out.kb_matrix= kb.out
 file.out.sa_matrix= sa.out
 file.out.retprofiles= rprfs.table
 file.out.aprprofiles= aprfs.table
 file.out.ab_matrix                           =ab.output
 file.out.summary= summary
 file.out.parm_vectors                        =parm.output
 file.out.seinv_vector= seinv.out
 file.out.sainv_matrix                        =sainv.complete
 file.out.smeas_matrix                        =smeas.target
 file.out.shat_matrix= shat.out
 file.out.ssmooth_matrix                      =ssmooth.target
 file.out.ak_matrix= ak.out
 file.out.solarspectrum                       =solar.output
 
 gas.profile.list                             = H2CO CH4 HDO O3
 gas.layers                                   = 48
 gas.profile.H2CO.correlation =  T
 gas.profile.H2CO.correlation.type =  6
 gas.profile.H2CO.correlation.lambda =  100.0
 gas.profile.H2CO.logstate = F
 gas.profile.H2CO.scale       = 1.0
 gas.profile.H2CO.sigma               =         # if TIK
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.profile.HDO.correlation =  T
 gas.profile.HDO.correlation.type =  6
 gas.profile.HDO.correlation.lambda =  100.0	
 gas.profile.HDO.logstate = F
 gas.profile.HDO.scale       = 1.0
 gas.profile.HDO.sigma               =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.profile.CH4.correlation =  T
 gas.profile.CH4.correlation.type =  6
 gas.profile.CH4.correlation.lambda =  100.0	
 gas.profile.CH4.logstate = F
 gas.profile.CH4.scale       = 1.0
 gas.profile.CH4.sigma               =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.profile.O3.correlation =  T
 gas.profile.O3.correlation.type =  6
 gas.profile.O3.correlation.lambda =  10.0	
 gas.profile.O3.logstate = F
 gas.profile.O3.scale       = 1.0
 gas.profile.O3.sigma               =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 
 gas.column.list             =  N2O
 gas.column.N2O.scale            = 1.0d0
 gas.column.N2O.sigma            = 1.0d1
 
 fw.isotope_separation                        =T
 fw.delnu                                     =0.1
 fw.lshapemodel                               =0
 fw.linemixing                                =T
 fw.linemixing.gas                            =CO2
 fw.solar_spectrum                            =T
 fw.pressure_shift                            =T
 fw.apod_fcn                                  =F
 fw.apod_fcn.type                             =2
 fw.apod_fcn.order                            =1 
 fw.phase_fcn                                 =F
 fw.phase_fcn.type                            =2
 fw.phase_fcn.order                           =1 
 fw.emission                                  =F
 fw.tips = F

 rt                                           =T
 rt.lm                                        = T
 rt.lm.gamma_start			      = 1e5
 rt.lm.gamma_inc			      = 10
 rt.lm.gamma_dec			      = 10
 rt.convergence                               =0.01
 rt.max_iteration                             =30
 rt.dwshift                                   =F
 rt.wshift                                    =T
 rt.wshift.type                               =2
 rt.wshift.apriori                            =0.0
 rt.wshift.sigma                              =0.1
 rt.slope                                     =T
 rt.slope.apriori                             =0.0
 rt.slope.sigma                               =0.1
 rt.curvature                                 =F
 rt.curvature.apriori                         =0.0
 rt.curvature.sigma                           =0.1
 rt.phase                                     =F
 rt.phase.apriori                             =0.0
 rt.phase.sigma                               =0.1
 rt.apod_fcn                                  =F
 rt.apod_fcn.apriori                          =1.0
 rt.apod_fcn.sigma                            =0.05
 rt.phase_fcn                                 =F
 rt.phase_fcn.apriori                         =0.0
 rt.phase_fcn.sigma                           =0.1
 rt.solshift                                  =T
 rt.solshift.apriori                          =0.0
 rt.solshift.sigma                            =0.1
 rt.solstrnth                                 =T
 rt.solstrnth.apriori                         =1.0
 rt.solstrnth.sigma                           =0.1
 rt.temperature                               =F
 rt.temperature.sigma                         =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 kb= F
 kb.temperature= T
 kb.slope                                     =T
 kb.curvature= T
 kb.solshft= T
 kb.solstrnth                                 =T
 kb.phase= T
 kb.wshift                                    =T
 kb.apod_fcn                                  =T
 kb.phase_fcn                                 =T
 kb.zshift= T
 kb.sza= T
 kb.omega= T
 kb.line                                      =T
 kb.line.type                                 =1
 kb.line.gas                                  = H2CO
 
 band                        =  1 2 3 4
 
 band.1.nu_start             =  2763.42
 band.1.nu_stop              =  2764.17
 band.1.zshift               =  F
 band.1.beam                 =  #5 6
 band.1.beam.model           = PS
 band.1.beam.5.apriori       = 2.0e-3 0.9 2763.81 0.0
 band.1.beam.5.sigma         = 0.0 0.0 0.0 0.0
 band.1.beam.6.apriori       = 0.08e-3 0.11 2763.81 0.0
 band.1.beam.6.sigma         = 0.0 0.0 0.0 0.0
 band.1.calc_point_space     =  0.500E-03
 band.1.wave_factor= 0.999998990000
 band.1.max_opd= 180.0000
 band.1.omega= 1.1962
 band.1.apodization_code     =  0
 band.1.gasb                 = H2CO CH4 O3 N2O HDO
 
 band.2.nu_start             =  2765.65
 band.2.nu_stop              =  2766.01
 band.2.zshift               =  F
 band.2.beam                 =  #5 6
 band.2.beam.model           = PS
 band.2.beam.5.apriori       = 2e-3 0.9 2763.81 0.0
 band.2.beam.5.sigma         = 0.0 0.0 0.0 0.0
 band.2.beam.6.apriori       = 0.08e-3 0.11 2763.81 0.0
 band.2.beam.6.sigma         = 0.0 0.0 0.0 0.0
 band.2.calc_point_space     =  0.500E-03
 band.2.wave_factor= 0.999998990000
 band.2.max_opd= 180.0000
 band.2.omega= 1.1962
 band.2.apodization_code     =  0
 band.2.gasb                 = H2CO CH4 O3 N2O HDO
 
 band.3.nu_start             =  2778.15
 band.3.nu_stop              =  2779.1
 band.3.zshift               =  F
 band.3.beam                 =  #5 6
  band.3.beam.model           = PS
 band.3.beam.5.apriori       = 2e-3 0.9 2763.81 0.0
 band.3.beam.5.sigma         = 0.0 0.0 0.0 0.0
 band.3.beam.6.apriori       = 0.08e-3 0.11 2763.81 0.0
 band.3.beam.6.sigma         = 0.0 0.0 0.0 0.0
 band.3.calc_point_space     =  0.500E-03
 band.3.wave_factor= 0.999998990000
 band.3.max_opd= 180.0000
 band.3.omega= 1.1962
 band.3.apodization_code     =  0
 band.3.gasb                 = H2CO CH4 O3 N2O HDO
 
 band.4.nu_start             =  2780.65
 band.4.nu_stop              =  2782.0
 band.4.zshift               =  F
 band.4.beam                 =  #5 6
 band.4.beam.model           = PS
 band.4.beam.5.apriori       = 2e-3 0.9 2763.81 0.0
 band.4.beam.5.sigma         = 0.0 0.0 0.0 0.0
 band.4.beam.6.apriori       = 0.08e-3 0.11 2763.81 0.0
 band.4.beam.6.sigma         = 0.0 0.0 0.0 0.0
 band.4.calc_point_space     =  0.500E-03
 band.4.wave_factor= 0.999998990000
 band.4.max_opd= 180.0000
 band.4.omega= 1.1962
 band.4.apodization_code     =  0
 band.4.gasb                 = H2CO CH4 O3 N2O HDO
 
 sp.snr                                       = 2 3
 sp.snr.1.nu_start                            =2750.0      # Fix SNR = 500
 sp.snr.1.nu_stop                             =2900.0
 sp.snr.1.snr                                 =500.0
 sp.snr.2.nu_start                            =2780.967    #O3
 sp.snr.2.nu_stop                             =2780.993
 sp.snr.2.snr                                 =1.0
 sp.snr.3.nu_start                            =2781.42     #CH4 small line right of the big one
 sp.snr.3.nu_stop                             =2781.48
 sp.snr.3.snr                                 =1.0
 
 out.level= 1
 out.gas_spectra= T
 out.gas_spectra.type                         =1
 out.sa_matrix= T
 out.statevec= T
 out.k_matrix= T
 out.shat_matrix= T
 out.retprofiles= T
 out.aprprofiles= T
 out.ab_matrix                                =F
 out.ak_matrix= T
 out.summary                                  =T
 out.pbpfile= T
 out.channel                                  =F
 out.parm_vectors                             =T
 out.seinv_vector= T
 out.sainv_matrix                             =F
 out.smeas_matrix                             =F
 out.ssmooth_matrix                           =F
 out.raytrace                                 =F
 out.raytrace.type                            =1
 out.solarspectrum                            =F
 out.levmardet                                =F
 out.xscdetail                                =F
 
 out.g_matrix= T
 file.out.g_matrix= g.out
