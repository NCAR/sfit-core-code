 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.input
 file.in.solarlines= solar.dat
 file.in.linelist= 00863.441722-00874.058278.hbin
 
 # Definition for retrieval gases
 
 gas.layers                  =              48
 gas.profile.list            = HNO3
 gas.column.list             = H2O OCS NH3
 gas.profile.HNO3.correlation          =               F
 gas.profile.HNO3.logstate             =               F
 gas.profile.HNO3.scale                =               1.0
 gas.profile.HNO3.sigma                =
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3   0.3   0.3
 0.3   0.3   0.3
 gas.column.H2O.scale               =             1.0
 gas.column.H2O.sigma               = 1.0
 gas.column.OCS.scale               = 1.0
 gas.column.OCS.sigma               = 1.0
 gas.column.NH3.scale               = 1.0
 gas.column.NH3.sigma               = 1.0
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum           = F
 fw.pressure_shift           = T
 fw.apod_fcn                 = F
 fw.phase_fcn                = F
 fw.isotope_separation       = F
 
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =               F
 rt.convergence              =           0.1
 rt.max_iteration            =              15
 rt.wshift                   =               T
 rt.wshift.type              = 3
 rt.wshift.apriori           =           0.000
 rt.wshift.sigma             =           0.100
 rt.slope                    =               T
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           0.100
 rt.curvature                =               T
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma          =           0.100
 rt.phase                    =               T
 rt.phase.apriori            =           0.000
 rt.phase.sigma              =           0.200
 rt.dwshift                  =              F
 
 # Microwindows and their parameters
 
 band                        =   1
 band.1.nu_start             =        867.5
 band.1.nu_stop              =        870.0
 band.1.zshift		     =	      F
 band.1.beam                 =               1
 band.1.beam.model           = IP
 band.1.beam.1.apriori       = 0.005 0.92 868.0 0.0
 band.1.beam.1.sigma         = 1.0 0.0 1.0 0.0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor          =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 3.5885
 band.1.apodization_code =               0
 band.1.gasb                 = HNO3 H2O OCS NH3
 
 out.level= 1
 out.smeas_matrix = T
 
 kb= T
 kb.slope = F
 kb.curvature = F
 kb.solshft= T
 kb.solstrnth= T
 kb.temperature= T
 kb.phase= F
 kb.omega= T
 kb.zshift= T
 kb.sza= T
 out.gas_spectra= F
 out.k_matrix= T
 out.ak_matrix= T
 out.g_matrix= T
 out.sa_matrix= T
 out.shat_matrix= T
 out.seinv_vector= T
 out.retprofiles= T
 out.aprprofiles= T
 out.pbpfile= T
 out.statevec= T
 file.out.pbpfile= pbpfile
 file.out.statevec= statevec
 file.out.ak_matrix= ak.out
 file.out.k_matrix= k.out
 file.out.g_matrix= g.out
 file.out.kb_matrix= kb.out
 file.out.sa_matrix= sa.out
 file.out.shat_matrix= shat.out
 file.out.seinv_vector= seinv.out
 file.out.retprofiles= rprfs.table
 file.out.aprprofiles= aprfs.table
 file.out.summary= summary
