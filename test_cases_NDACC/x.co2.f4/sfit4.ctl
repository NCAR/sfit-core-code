 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.sa_matrix              = ny_alesund.dat
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.input
 file.in.solarlines= ../../linelist/solar/120621/solar.dat
 file.in.linelist= 02616.491422-02634.008578.hbin
 
 # Definition for retrieval gases
 
 
 gas.layers                  =              41
 gas.isoflag                 =               F
 gas.profile.list                 = CO2
 gas.column.list                  = H2O HDO CH4
 gas.profile.CO2.correlation        =               T
 gas.profile.CO2.correlation.type        =               1
 gas.profile.CO2.correlation.width                =         1000.0
 gas.profile.CO2.correlation.minalt                =             0.0
 gas.profile.CO2.correlation.maxalt                =           120.0
 gas.profile.CO2.logstate             =               F
 gas.profile.CO2.scale = 1.0		
 gas.profile.CO2.sigma                =
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1
 gas.column.H2O.scale               = 1.0
 gas.column.H2O.sigma               = 1.0
 gas.column.HDO.scale               = 1.0
 gas.column.HDO.sigma               = 1.0
 gas.column.CH4.scale               = 1.0
 gas.column.CH4.sigma               = 1.0
 gas.column.CO2.scale               = 1.0
 gas.column.CO2.sigma               = 1.0
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum           =               T
 fw.pressure_shift           = T
 fw.apod_fcn                 = F
 fw.phase_fcn                = F
 fw.isotope_separation       = T
 
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =               F
 rt.lm.gamma_start	     =           1.0e3
 rt.lm.gamma_inc             =            10.0
 rt.lm.gamma_dec             =            10.0
 rt.convergence              =             0.1
 rt.max_iteration            =              17
 rt.wshift = T
 rt.wshift.type                   =               3
 rt.wshift.apriori           =          0.0
 rt.wshift.sigma                =           0.100
 rt.slope                    = T
 rt.slope.apriori            =           0.000
 rt.slope.sigma                 =           0.100
 rt.curvature                = F
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma             =           0.100
 rt.phase                    = T
 rt.phase.apriori            = 0.0
 rt.phase.sigma              = 0.1
 rt.dwshift                   =               F
 rt.solshift = F
 rt.solshift.apriori            =             0.0
 rt.solshift.sigma                 =             1.0
 rt.temperature              =               T
 rt.temperature.sigma        = 
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001
 
 # Microwindows and their parameters
 
 band                        =   1 2 3 4
 band.1.nu_start             =        2620.55
 band.1.nu_stop              =        2621.1
 band.1.zshift               =               F
 band.1.zshift.apriori       =             0.0
 band.1.zshift.sigma            =             0.2
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.900E-03
 band.1.wave_factor               =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.3923
 band.1.apodization_code     =               0
 band.1.gasb                 = H2O HDO CH4 CO2
 band.1.tempretb = T
 
 band.2.nu_start             =        2626.4
 band.2.nu_stop              =        2626.85
 band.2.zshift               =               F
 band.2.zshift.apriori       =             0.0
 band.2.zshift.sa            =             0.2
 band.2.beam                 =               0
 band.2.calc_point_space                   =       0.900E-03
 band.2.wave_factor               =           1.000
 band.2.max_opd= 180.0000
 band.2.omega= 2.3923
 band.2.apodization_code                  =               0
 band.2.gasb                 = H2O HDO CH4 CO2
 

 band.3.nu_start             =        2627.1
 band.3.nu_stop              =        2627.6
 band.3.zshift               =               F
 band.3.zshift.apriori       =             0.0
 band.3.zshift.sa            =             0.2
 band.3.beam                 =               0
 band.3.calc_point_space     =       0.900E-03
 band.3.wave_factor        =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 2.3923
 band.3.apodization_code                  =               0
 band.3.gasb                 = H2O HDO CH4 CO2
 
 band.4.nu_start             =        2629.275
 band.4.nu_stop              =        2629.95
 band.4.zshift               =               F
 band.4.zshift.apriori       =             0.0
 band.4.zshift.sa            =             0.2
 band.4.beam                 =               0
 band.4.calc_point_space     =       0.900E-03
 band.4.wave_factor        =           1.000
 band.4.max_opd= 180.0000
 band.4.omega= 2.3923
 band.4.apodization_code                  =               0
 band.4.gasb                 = H2O HDO CH4 CO2
 
 out.level= 1
 out.gas_spectra= T
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
 kb= T
 kb.slope= T
 kb.curvature= T
 kb.solshft= T
 kb.solstrnth= T
 kb.temperature= T
 kb.phase= T
 kb.omega= T
 kb.max_opd= T
 kb.zshift= T
 kb.sza= T
