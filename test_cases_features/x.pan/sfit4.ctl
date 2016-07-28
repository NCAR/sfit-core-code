 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.inp
 file.in.solarlines             = ../../linelist/solar/120621/solar.dat
 file.in.linelist               = 00755.947222-00874.052778.hbin
 
 
 
 # Definition for retrieval gases
 
 gas.layers                  =              50
 gas.isoflag                 =               F
 gas.profile.list            = PAN 
 gas.column.list             = H2O CO2 O3 HNO3 CH4
 gas.profile.PAN.correlation =               F
 gas.profile.PAN.logstate = F
 gas.profile.PAN.scale       = 1.0               
 gas.profile.PAN.sigma               =
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
 gas.column.H2O.scale          =  1.0
 gas.column.H2O.sigma          =  1.0
 gas.column.CO2.scale               =               1.0
 gas.column.CO2.sigma               = 1.0
 gas.column.O3.scale               = 1.0
 gas.column.O3.sigma               = 1.0
 gas.column.N2O.scale              =  1.0
 gas.column.N2O.sigma               = 1.0
 gas.column.CH4.scale              =  1.0
 gas.column.CH4.sigma               = 1.0
 
 kb =F
 kb.curvature = T
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               3
 fw.linemixing = F
 fw.pressure_shift           =               T
 fw.emission                 =               F
 fw.isotope_separation       =               F
 # Retrieval parameter
 
 rt                          =               F
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
 
 band                        =  1
 band.1.nu_start             =       760.0
 band.1.nu_stop              =       870.0
 band.1.zshift               =               F
 band.1.zshift.type          =               2
 band.1.zshift.apriori       =             0.0
 band.1.zshift.sigma            =             0.2
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-01
 band.1.wave_factor               =           1.000
 band.1.max_opd = 180.0000
 band.1.omega= 4.0670
 band.1.apodization_code                  =               0
 band.1.gasb                 = PAN H2O CO2 O3 HNO3 CH4
 
 
 out.level = 1
 out.gas_spectra = T
 out.gas_spectra.type = 1	
