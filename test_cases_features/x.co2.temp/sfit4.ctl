 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = solar.dat
 file.in.linelist               = 00895.941722-01004.058278.hbin
 
 # Definition for retrieval gases
 
 gas.layers                  =              41
 gas.profile.list            = CO2 H2O O3
 gas.column.list             = 
 gas.profile.CO2.correlation          =               T
 gas.profile.CO2.correlation.type     =               1
 gas.profile.CO2.correlation.lambda     =               1e6	
 gas.profile.CO2.correlation.minalt  =               0.0
 gas.profile.CO2.correlation.maxalt  =               150.0
 gas.profile.CO2.correlation.width    =               1000.0
 gas.profile.CO2.logstate             =               F
 gas.profile.CO2.scale                =               1.0
 gas.profile.CO2.sigma                =
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 0.01   0.01   0.01   0.01   0.01
 1.0
 gas.profile.H2O.correlation          =               T
 gas.profile.H2O.correlation.type     =               6
 gas.profile.H2O.correlation.lambda     =               1e6		
 gas.profile.H2O.correlation.minalt  =               0.0
 gas.profile.H2O.correlation.maxalt  =               150.0
 gas.profile.H2O.correlation.width    =               1000.0
 gas.profile.H2O.logstate             =               F
 gas.profile.H2O.scale                =               1.0
 gas.profile.H2O.sigma                =
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.0
 gas.profile.O3.correlation         =               F
 gas.profile.O3.correlation.type     =               6
 gas.profile.O3.correlation.lambda     =               1e6
 gas.profile.O3.correlation.minalt  =               0.0
 gas.profile.O3.correlation.maxalt  =               150.0
 gas.profile.O3.correlation.width    =               1000.0
 gas.profile.O3.logstate            =               F
 gas.profile.O3.scale               =             1.0
 gas.profile.O3.sigma               =
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.00000   1.00000   1.00000   1.00000   1.00000
 1.0
 gas.column.CO2.scale               = 1.0
 gas.column.CO2.sigma               = 1.0
 gas.column.O3.scale               = 1.0
 gas.column.O3.sigma               = 1.0
 gas.column.H2O.scale               = 1.0
 gas.column.H2O.sigma               = 1.0
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum           = T
 fw.pressure_shift           = T
 fw.apod_fcn                 = F
 fw.phase_fcn                = F
 fw.isotope_separation       = F
 fw.tips = F

 # Retrieval parameter
 
 rt                          =              T 
 rt.lm                       =               T
 rt.lm.gamma_start            =               1.0e5
 rt.lm.gamma_inc                       =  10.0
 rt.lm.gamma_dec                       =            10.0
 rt.convergence              =           0.1
 rt.max_iteration            =              15
 rt.wshift                   =               T
 rt.wshift.type              = 3
 rt.wshift.apriori           =           0.000
 rt.wshift.sigma             =           0.100
 rt.slope                    =               T
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           0.100
 rt.curvature                =               F
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma          =           0.100
 rt.solshift = F
 rt.solshift.apriori = 0.0
 rt.solshift.sigma   = 0.1
 rt.phase                    =               T
 rt.phase.apriori            =           0.000
 rt.phase.sigma              =           0.200
 rt.dwshift                  =               F
 rt.temperature              = T
 rt.temperature.lambda       = 1e5
 
 # Microwindows and their parameters
 
 band                        =   1
 band.1.nu_start             =        920.00
 band.1.nu_stop              =        960.00
 band.1.zshift		     =	      F
 band.1.beam                 =               
 band.1.calc_point_space     =       1.00E-3
 band.1.wave_factor          =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.3923
 band.1.apodization_code =               0
 band.1.gasb                 = CO2 H2O O3
 band.1.tempretb             = T  

 out.gas_spectra= T

 kb= F
 kb.slope= T
 kb.curvature= T
 kb.solshft= T
 kb.solstrnth= T
 kb.temperature= T
 kb.phase= T
 kb.omega= T
 kb.zshift= T
 kb.sza= T

 out.level = 1
 out.retprofiles                     = T
 out.aprprofiles                     = T
