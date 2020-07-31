 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.inp
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 01895.821722-01916.908278.hbin
 
 # Definition for retrieval gases
 
 gas.layers                  =              41

 gas.column.list                         = O3 CO2 H2O N2O 
 gas.profile.list                         = NO
 
 gas.profile.NO.correlation                =               T
 gas.profile.NO.correlation.type                =               2
 gas.profile.NO.correlation.width                =             4.0
 gas.profile.NO.correlation.minalt               =             0.0
 gas.profile.NO.correlation.maxalt                =           120.0
 gas.profile.NO.logstate             =               F
 gas.profile.NO.scale = 1.0
 gas.profile.NO.sigma                =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0

 gas.profile.O3.correlation                =              F 
 gas.profile.O3.logstate             =               F
 gas.profile.O3.sigma                =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0
 
 gas.column.O3.scale = 1.0
 gas.column.O3.sigma                = 1.0
 
 gas.column.H2O.scale                = 0.9
 gas.column.H2O.sigma                = 0.9
 
 gas.column.CO2.scale                = 0.8
 gas.column.CO2.sigma                = 0.8
 
 gas.column.N2O.sigma                = 0.7
 gas.column.N2O.scale = 0.7
 
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum	     =               T
 fw.pressure_shift           =               F
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.isotope_separation       =               F
 # Retrieval parameter

 kb = T
 rt                          =               T
 rt.lm                       =               T
 rt.lm.gamma_start	     =           1.0e3
 rt.lm.gamma_inc             =            10.0
 rt.lm.gamma_dec             =            10.0
 rt.convergence              =             0.1
 rt.max_iteration            =              17
 rt.wshift                   =               T
 rt.wshift.type              =               3
 rt.wshift.apriori           =             0.0
 rt.wshift.sigma                =             0.1
 rt.slope                    = T
 rt.slope.apriori            =           0.000
 rt.slope.sigma                 =           0.100
 rt.curvature                = F
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma             =           0.100
 rt.solshift                    = T
 rt.solshift.apriori            =             0.0
 rt.solshift.sigma                 =             1.0
 rt.dwshift                   =               F
 rt.temperature              =               F
 
 # Microwindows and their parameters
 
 band                        =   1  2 3
 band.1.nu_start             =  1899.880
 band.1.nu_stop              =        1900.150
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor           =           1.000
 band.1.max_opd    = 180.0000
 band.1.omega  = 2.3923
 band.1.apodization_code                  =               0
 band.1.gasb                 = NO CO2 H2O N2O O3
 band.2.nu_start             =  1903.050
 band.2.nu_stop              =        1903.260
 band.2.zshift               =               F
 band.2.beam                 =               0
 band.2.calc_point_space                   =       0.500E-03
 band.2.wave_factor               =           1.000
 band.2.max_opd = 180.0000
 band.2.omega= 2.3923
 band.2.apodization_code                   =               0
 band.2.gasb                 = NO CO2 H2O O3
 band.3.nu_start             =  1912.742
 band.3.nu_stop              =        1912.850
 band.3.zshift               =               F
 band.3.beam                 =               0
 band.3.calc_point_space                    =       0.500E-03
 band.3.wave_factor               =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 2.3923
 band.3.apodization_code                  =               0
 band.3.gasb                 = NO CO2 H2O N2O
 
 out.level = 2
 out.gas_spectra = T
 out.gas_spectra.type = 1
