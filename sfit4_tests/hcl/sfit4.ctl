 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.input
 file.in.solarlines= solar.dat
 file.in.linelist= 02723.689356-02930.040644.hbin
 
 # Definition for retrieval gases
 
 
 gas.layers                  =              48
 gas.isoflag                 =               F
 gas.profile.list                 = HCL O3
 gas.column.list                  = HDO CH4 H2O N2O NO2 CO2
 gas.profile.HCL.correlation        =               F
 gas.profile.HCL.correlation.type        =               1
 gas.profile.HCL.correlation.width                =             5.0
 gas.profile.HCL.correlation.minalt                =             0.0
 gas.profile.HCL.correlation.maxalt                =           120.0
 gas.profile.HCL.logstate             =               F
 gas.profile.HCL.scale = 1.0		
 gas.profile.HCL.sigma                =
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 gas.profile.O3.correlation        =               F
 gas.profile.O3.correlation.type        =               1
 gas.profile.O3.correlation.width                =             5.0
 gas.profile.O3.correlation.minalt                =             0.0
 gas.profile.O3.correlation.maxalt                =           120.0
 gas.profile.O3.logstate             =               F
 gas.profile.O3.scale = 1.0		
 gas.profile.O3.sigma                =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.profile.H2O.correlation        =               F
 gas.profile.H2O.correlation.type        =               1
 gas.profile.H2O.correlation.width                =             5.0
 gas.profile.H2O.correlation.minalt                =             0.0
 gas.profile.H2O.correlation.maxalt                =           120.0
 gas.profile.H2O.logstate             =               F
 gas.profile.H2O.scale = 1.0		
 gas.profile.H2O.sigma                =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.column.HDO.scale               =            0.00001
 gas.column.HDO.sigma               = 1.0
 gas.column.CH4.scale               =             1.0
 gas.column.CH4.sigma               = 1.0
 gas.column.N2O.scale               =             1.0
 gas.column.N2O.sigma               = 1.0
 gas.column.O3.scale = 1.0
 gas.column.O3.sigma               = 1.0
 gas.column.NO2.scale               =             1.0
 gas.column.NO2.sigma               = 1.0
 gas.column.CO2.scale               =             1.0
 gas.column.CO2.sigma               = 1.0
 gas.column.H2O.scale               =             1.0
 gas.column.H2O.sigma               = 1.0
 
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
 rt.curvature                = T
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma             =           0.100
 rt.dwshift                   =               F
 rt.solshift = T
 rt.solshift.apriori            =             0.0
 rt.solshift.sigma                 =             1.0
 rt.temperature              =               F
 
 # Microwindows and their parameters
 
 band                        =   1 2 3
 band.1.nu_start             =        2727.7300
 band.1.nu_stop              =        2727.8300
 band.1.zshift               =               F
 band.1.zshift.apriori       =             0.0
 band.1.zshift.sigma            =             0.2
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.900E-03
 band.1.wave_factor               =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 1.1962
 band.1.apodization_code     =               0
 band.1.gasb                 = HCL O3 HDO CH4 CO2
 
 band.2.nu_start             =        2775.7000
 band.2.nu_stop              =        2775.8000
 band.2.zshift               =               F
 band.2.zshift.apriori       =             0.0
 band.2.zshift.sa            =             0.2
 band.2.beam                 =               0
 band.2.calc_point_space                   =       0.900E-03
 band.2.wave_factor               =           1.000
 band.2.max_opd= 180.0000
 band.2.omega= 1.1962
 band.2.apodization_code                  =               0
 band.2.gasb                 = HCL O3 N2O CH4 CO2 HDO
 
 band.3.nu_start             =        2925.800
 band.3.nu_stop              =        2926.00
 band.3.zshift               =               F
 band.3.zshift.apriori       =             0.0
 band.3.zshift.sa            =             0.2
 band.3.beam                 =               0
 band.3.calc_point_space     =       0.900E-03
 band.3.wave_factor        =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 1.1962
 band.3.apodization_code                  =               0
 band.3.gasb                 = HCL CH4 NO2 H2O
 
 out.level= 1
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
 kb= T
 kb.slope = T
 kb.curvature = T
 kb.solshft = T
 kb.solstrnth= T
 kb.temperature= T
 kb.phase= T
 kb.omega= T
 kb.zshift= T
 kb.sza= T
