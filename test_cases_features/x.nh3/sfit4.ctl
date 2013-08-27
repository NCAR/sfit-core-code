 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.inp
 file.in.solarlines             = ../../linelist/solar/120621/solar.dat
 file.in.linelist               = 00903.541722-01107.758278.hbin
 
 
 
 # Definition for retrieval gases
 
 gas.layers                  =              50
 gas.isoflag                 =               F
 gas.profile.list            = NH3 O3 
 gas.column.list             = H2O N2O CO21 CO22 CO23 HNO3
 gas.profile.NH3.correlation =               F
 gas.profile.NH3.logstate = F
 gas.profile.NH3.scale       = 1.0               
 gas.profile.NH3.sigma               =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.profile.HNO3.correlation         = F
 gas.profile.HNO3.scale               =          1.0
 gas.profile.HNO3.sigma               =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.column.H2O.scale        =        1.0
 gas.column.H2O.sigma        = 1.0
 gas.column.HNO3.scale        =        1.0
 gas.column.HNO3.sigma        = 1.0
 gas.profile.H2O.correlation         = F
 gas.profile.H2O.scale        =        1.0
 gas.profile.H2O.sigma        =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.profile.O3.correlation =               F
 gas.profile.O3.scale               = 1.0
 gas.profile.O3.sigma               =
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
 gas.column.CO21.scale          =  1.0
 gas.column.CO21.sigma          =  1.0
 gas.column.CO22.scale               =               1.0
 gas.column.CO22.sigma               = 1.0
 gas.column.CO23.scale               = 1.0
 gas.column.CO23.sigma               = 1.0
 gas.column.N2O.scale              =  1.0
 gas.column.N2O.sigma               = 1.0
 
 kb =F
 kb.curvature = T
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               3
 fw.linemixing = F
 fw.pressure_shift           =               T
 fw.emission                 =               F
 fw.isotope_separation       =               T
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =               F
 rt.lm.gamma_start	     =           1.0e5
 rt.lm.gamma_inc             =            10.0
 rt.lm.gamma_dec             =            10.0
 rt.convergence              =             0.1
 rt.max_iteration            =              17
 rt.wshift                   =               T
 rt.wshift.type = 3
 rt.wshift.apriori           =             0.0
 rt.wshift.sigma                =             0.1
 rt.slope                    = T
 rt.slope.apriori            =           0.000
 rt.slope.sigma                 =           0.100
 rt.curvature                = F
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma             =           0.100
 rt.solshift                    = F
 rt.solshift.apriori = 0.0
 rt.solshift.sigma = 0.1
 
 # Microwindows and their parameters
 
 band                        =  1 2 3 4 5 6
 band.1.nu_start             =        907.6
 band.1.nu_stop              =        908.6
 band.1.zshift               =               F
 band.1.zshift.type          =               2
 band.1.zshift.apriori       =             0.0
 band.1.zshift.sigma            =             0.2
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor               =           1.000
 band.1.max_opd = 180.0000
 band.1.omega= 4.0670
 band.1.apodization_code                  =               0
 band.1.gasb                 = NH3   HNO3 H2O CO21 CO22 CO23
 
 band.2.nu_start             =        930.43
 band.2.nu_stop              =        931.12
 band.2.zshift               =               F
 band.2.zshift.type          =               2
 band.2.zshift.apriori       =             0.0
 band.2.zshift.sigma            =             0.2
 band.2.beam                 =               0
 band.2.calc_point_space                   =       0.500E-03
 band.2.wave_factor               =           1.000
 band.2.max_opd= 180.0000
 band.2.omega= 4.0670
 band.2.apodization_code                  =               0
 band.2.gasb                 = NH3   N2O CO21 CO23
 
 band.3.nu_start             =        964.94
 band.3.nu_stop              =        965.72
 band.3.zshift               =               F
 band.3.zshift.type          =               2
 band.3.zshift.apriori       =             0.0
 band.3.zshift.sidma            =             0.2
 band.3.beam                 =               0
 band.3.calc_point_space                   =       0.500E-03
 band.3.wave_factor               =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 4.0670
 band.3.apodization_code                  =               0
 band.3.gasb                 = NH3 O3     H2O    N2O CO21 CO23
 
 band.4.nu_start             =        966.97
 band.4.nu_stop              =        967.67
 band.4.zshift               =              F
 band.4.zshift.type = 2
 band.4.zshift.apriori       =             0.0
 band.4.zshift.sigma            =             0.2
 band.4.beam                 =               0
 band.4.calc_point_space     =       0.500E-03
 band.4.wave_factor               =           1.000
 band.4.max_opd = 180.0000
 band.4.omega= 4.0670
 band.4.apodization_code     =               0
 band.4.gasb                 = NH3 H2O   O3 CO21 CO22 CO23
 
 band.5.nu_start             =        1046.2
 band.5.nu_stop              =        1046.6
 band.5.zshift               =              F
 band.5.zshift.type = 1
 band.5.zshift.apriori       =             0.0
 band.5.zshift.sigma            =             0.2
 band.5.beam                 =               0
 band.5.calc_point_space                   =       0.500E-03
 band.5.wave_factor               =           1.000
 band.5.max_opd= 180.0000
 band.5.omega= 4.0670
 band.5.apodization_code                  =               0
 band.5.gasb                 = NH3  O3 CO21 CO23
 
 band.6.nu_start             =        1103.3
 band.6.nu_stop              =        1103.7
 band.6.zshift               =               F
 band.6.zshift.type               =               2
 band.6.zshift.apriori       =             0.0
 band.6.zshift.sigma            =             0.2
 band.6.beam                 =               0
 band.6.calc_point_space                   =       0.500E-03
 band.6.wave_factor               =           1.000
 band.6.max_opd= 180.0000
 band.6.omega= 4.0670
 band.6.apodization_code                  =               0
 band.6.gasb                 = NH3 H2O O3
 
 out.level = 1
 out.gas_spectra = T
 out.gas_spectra.type = 2	
