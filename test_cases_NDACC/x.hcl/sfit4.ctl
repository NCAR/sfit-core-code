 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.input
<<<<<<< HEAD
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 02723.671422-02930.058578.hbin
 file.out.ak_matrix= ak.out
 file.out.ab_matrix             = AB.out
 file.out.summary= summary
 file.out.k_matrix= k.out
=======
 file.in.solarlines             = ../../linelist/solar/120621/solar.dat
 file.in.linelist               = 02723.671422-02930.058578.hbin
 file.out.ak_matrix             = ak.target
 file.out.ab_matrix             = AB.out
 file.out.summary               = summary
 file.out.k_matrix              = k.out
>>>>>>> Pre-Release
 file.out.smeas_matrix          = SMEAS.out
 
 # Definition for retrieval gases
 
 
 gas.layers                  =              41
<<<<<<< HEAD
 gas.profile.list                 = HCL
 gas.column.list                  = CH4 NO2
=======
 gas.profile.list                 = HCL 
 gas.column.list                  = CH4 NO2 O3 N2O HDO
>>>>>>> Pre-Release
 gas.profile.HCL.correlation        =               F
 gas.profile.HCL.correlation.type        =               2
 gas.profile.HCL.correlation.width        =               4.0
 gas.profile.HCL.correlation.minalt        =               0.0
 gas.profile.HCL.correlation.maxalt        =               120.0
 gas.profile.HCL.logstate             =               F
 gas.profile.HCL.scale = 1.0		
 gas.profile.HCL.sigma                =
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5
 0.5
 gas.column.HDO.logstate            =               F
 gas.column.HDO.scale               = 1.0
 gas.column.HDO.sigma               = 1.0
 gas.column.CH4.logstate            =               F
 gas.column.CH4.scale               =             1.0
 gas.column.CH4.sigma               = 1.0
<<<<<<< HEAD
 gas.column.H2O.logstate            =               F
 gas.column.H2O.scale               =             1.0
 gas.column.H2O.sigma               = 1.0
=======
 gas.column.N2O.logstate            =               F
 gas.column.N2O.scale               =             1.0
 gas.column.N2O.sigma               = 1.0
>>>>>>> Pre-Release
 gas.column.O3.logstate            =               F
 gas.column.O3.scale = 1.0
 gas.column.O3.sigma               = 1.0
 gas.column.NO2.logstate            =               F
 gas.column.NO2.scale               =             1.0
 gas.column.NO2.sigma               = 1.0
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
<<<<<<< HEAD
 fw.lshapemodel              =              0 
=======
 fw.lshapemodel              =               0
>>>>>>> Pre-Release
 fw.solar_spectrum           =               T
 fw.pressure_shift          =               T
 fw.isotope_separation                    =               T
 # Retrieval parameter
 
<<<<<<< HEAD
 rt                          =              F 
=======
 rt                          =               T
>>>>>>> Pre-Release
 rt.lm                       =               F
 rt.lm.gamma_start	     =           1.0e3
 rt.lm.gamma_inc             =            10.0
 rt.lm.gamma_dec             =            10.0
 rt.convergence              =             0.1
 rt.max_iteration            =              17
 rt.wshift = T
 rt.wshift.type                   =               3
<<<<<<< HEAD
 rt.wshift.apriori           =          4.8848E-06
=======
 rt.wshift.apriori           =          4.8316E-06
>>>>>>> Pre-Release
 rt.wshift.sigma                =           0.100
 rt.slope                    = T
 rt.slope.apriori            =           0.000
 rt.slope.sigma                 =           0.100
<<<<<<< HEAD
 rt.curvature                = T
=======
 rt.curvature                = F
>>>>>>> Pre-Release
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma             =           0.100
 rt.dwshift                   =               F
 rt.temperature              =               F
 
<<<<<<< HEAD
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
 
 # Microwindows and their parameters
 
 band                        =    1 2 3
=======
 # Microwindows and their parameters
 
 band                        = 1 2 3
>>>>>>> Pre-Release
 band.1.nu_start             =        2727.7300
 band.1.nu_stop              =        2727.8300
 band.1.zshift               =               F
 band.1.zshift.apriori       =             0.0
 band.1.zshift.sigma            =             0.2
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.900E-03
 band.1.wave_factor               =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.3932
 band.1.apodization_code     =               0
<<<<<<< HEAD
 band.1.gasb                 = HCL  O3 HDO
=======
 band.1.gasb                 = HCL O3 HDO
>>>>>>> Pre-Release
 band.2.nu_start             =        2775.7000
 band.2.nu_stop              =        2775.8000
 band.2.zshift               =               F
 band.2.zshift.apriori       =             0.0
 band.2.zshift.sa            =             0.2
 band.2.beam                 =               0
 band.2.calc_point_space     =       0.900E-03
 band.2.wave_factor          =           1.000
 band.2.max_opd= 180.0000
<<<<<<< HEAD
 band.2.omega= 1.9139
=======
 band.2.omega= 2.3932
>>>>>>> Pre-Release
 band.2.apodization_code     =               0
 band.2.gasb                 = HCL O3 N2O
 band.3.nu_start             =        2925.800
 band.3.nu_stop              =        2926.00
 band.3.zshift               =               F
 band.3.zshift.apriori       =             0.0
 band.3.zshift.sa            =             0.2
 band.3.beam                 =               0
 band.3.calc_point_space     =       0.900E-03
 band.3.wave_factor          =           1.000
 band.3.max_opd= 180.0000
<<<<<<< HEAD
 band.3.omega= 1.9139
 band.3.apodization_code     =               0
 band.3.gasb                 = HCL CH4 NO2
 
 out.level= 1
 out.gas_spectra= F
 out.ak_matrix= T
 out.k_matrix= T
 out.g_matrix= T
 out.sa_matrix= T
 out.retprofiles= T
 out.aprprofiles= T
 file.out.g_matrix= g.out
 file.out.kb_matrix= kb.out
 file.out.sa_matrix= sa.out
 file.out.retprofiles= rprfs.table
 file.out.aprprofiles= aprfs.table
=======
 band.3.omega= 2.3932
 band.3.apodization_code     =               0
 band.3.gasb                 = HCL CH4 NO2 
 
 out.level = 1
 out.gas_spectra = F

>>>>>>> Pre-Release
