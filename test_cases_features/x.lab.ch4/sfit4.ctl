 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.isotope                = isotope.input
 file.in.spectrum               = CH4-Mix-HBr-DLR-0300mbar-296K-040m-BL.dpt
 file.in.solarlines             = solar.dat
 file.in.linelist               = 02609.641672-02925.658328.hbin
 
 # Definition for retrieval gases
 
 gas.layers                  =              1
 gas.column.list                         = CH4
 gas.profile.list                         = 
 gas.column.CH4.logstate = F
 gas.column.CH4.scale               = 1.0
 gas.column.CH4.sigma               = 1.0
 gas.column.HBR.logstate = F
 gas.column.HBR.scale               = 1.0
 gas.column.HBR.sigma               = 1.0
 gas.column.CO2.logstate = F
 gas.column.CO2.scale               = 1.0
 gas.column.CO2.sigma               = 1.0
 gas.column.H2O.logstate = F
 gas.column.H2O.scale               = 1.0
 gas.column.H2O.sigma               = 1.0
   
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.isotope_separation       =               T
 fw.lab = T
 fw.lab.length = 40.0
 fw.lab.pressure = 300.0
 fw.lab.temperature = 296.0
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                               = F
 rt.lm.gamma_start                   = 1.0e5
 rt.lm.gamma_inc                     = 10.0
 rt.lm.gamma_dec                     = 10.0
 rt.convergence                      = 0.01
 rt.tolerance                        = 0.050
 rt.max_iteration                    = 30
 rt.dwshift                          = F
 rt.wshift                           = T
 rt.wshift.type                      = 3
 rt.wshift.apriori                   = 0.000
 rt.wshift.sigma                     = 0.100
 rt.slope                            = F
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
 rt.apod_fcn.sigma                   = 0.000
 rt.phase_fcn                        = F
 rt.phase_fcn.apriori                = 0.000
 rt.phase_fcn.sigma                  = 0.000
 rt.solshift                         = F
 rt.solstrnth                        = F
 rt.temperature                      = T
 rt.temperature.sigma                = 0.1

 # Microwindows and their parameters
 
 band                        =   1 3 4
 band.1.nu_start             =        2613.700
 band.1.nu_stop              =        2615.400
 band.1.zshift               =              F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-04
 band.1.wave_factor          =      1.0
 band.1.max_opd		     =        180.0000
 band.1.omega	             =		2.409
 band.1.apodization_code     =               0
 band.1.gasb                 = CH4
 band.1.tempretb = T
 band.3.nu_start             =        2835.500
 band.3.nu_stop              =        2835.800
 band.3.zshift               =              F
 band.3.beam                 =               0
 band.3.calc_point_space     =       0.500E-04
 band.3.wave_factor          =           1.0
 band.3.max_opd		     = 	      180.0000
 band.3.omega		     =          2.409
 band.3.apodization_code     =               0
 band.3.gasb                 = CH4  
 band.3.tempretb = T

 band.4.nu_start             =        2921.000
 band.4.nu_stop              =        2921.600
 band.4.zshift               =              F
 band.4.beam                 =               0
 band.4.calc_point_space     =       0.500E-04
 band.4.wave_factor          =           1.0
 band.4.max_opd		     =        180.0000
 band.4.omega                = 2.409
 band.4.apodization_code     =               0
 band.4.snr= 570.3035
 band.4.gasb                 = CH4 
 band.4.tempretb = T
 out.level = 1