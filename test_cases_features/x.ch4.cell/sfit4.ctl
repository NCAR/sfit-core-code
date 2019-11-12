 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = CH4_AIR.ref
 file.in.isotope                = isotope.input
# file.in.spectrum                =  CH4-Mix-DLR-1000mbar-296K-040m.dpt
 file.in.spectrum                =  CH4-Mix-DLR-0030mbar-296K-040m.dpt
 file.in.solarlines             = solar.dat
# file.in.linelist               = 02390.441722-02709.558278.hbin.kit01.all
# file.in.linelist               = 02390.441722-03008.558278.hbin.009
#  file.in.linelist               = 02390.441722-02709.558278.hbin.kit01.all
 file.in.linelist               = 02916.941672-02925.658328.hbin 
# Definition for retrieval gases
 
 gas.layers                  =              0
 gas.column.list                         = CH4 
 gas.profile.list                         = 
 gas.column.CH4.logstate = F
 gas.column.CH4.scale               = 1.0
 gas.column.CH4.sigma               = 1.0
 gas.column.CH42.logstate = F
 gas.column.CH42.scale               = 1.0
 gas.column.CH42.sigma               = 1.0
 gas.column.CH43.logstate = F
 gas.column.CH43.scale               = 1.0
 gas.column.CH43.sigma               = 1.0
 gas.column.HBR.logstate = F
 gas.column.HBR.scale               = 1.0
# Definition for retrieval gases
 
   
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
# fw.lshapemodel.sdv          =               T
 fw.linemixing = T
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               F
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.isotope_separation       =               F
## CH4-DLR-00p02mbar-296K-040m.dpt
# fw.lab = T
# fw.lab.length = 40.0
# fw.lab.pressure = 0.02
# fw.lab.temperature = 296.0
## CH4-Mix-DLR-1000mbar-296K-040m.dpt
#  fw.lab = T
#  fw.lab.length = 40.0
#  fw.lab.pressure = 1000.0
#  fw.lab.temperature = 296.0
## CH4-Mix-DLR-0030mbar-296K-040m.dpt
# fw.lab = T
# fw.lab.length = 40.0
# fw.lab.pressure = 30.0
# fw.lab.temperature = 296.0
## CH4-Mix-DLR-0100mbar-296K-040m.dpt
 # fw.lab = T
 # fw.lab.length = 40.0
 # fw.lab.pressure = 100.0
 # fw.lab.temperature = 296.0
 # Retrieval parameter

 cell = 1
 cell.1.temperature = 296
 cell.1.pressure = 30.54
 cell.1.gas = CH4
 cell.1.vmr = 0.00498
 cell.1.path = 4017.5

 rt                          = F
 rt.lm                               = F
 rt.lm.gamma_start                   = 1.0e5
 rt.lm.gamma_inc                     = 10.0
 rt.lm.gamma_dec                     = 10.0
 rt.convergence                      = 0.1
 rt.max_iteration                    = 30
 rt.dwshift                          = F
 rt.wshift                           = T
 rt.wshift.type                      = 3
 rt.wshift.apriori                   = 0.0
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
 rt.temperature                      = F
 rt.temperature.sigma                = 0.1

 # Microwindows and their parameters
 
 band                        =   4
 band.1.nu_start             =        2900.000
 band.1.nu_stop              =        2950.000
 band.1.zshift               =              F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor          =      1.0
 band.1.max_opd		     =        180.0000
 band.1.omega	             =		2.3923
 band.1.apodization_code     =               0
 band.1.gasb                 = CH4
 band.1.tempretb =T

# band.2.nu_start             =        2794.9000
# band.2.nu_stop              =        2795.100
 band.2.nu_start             =        2790.000
 band.2.nu_stop              =        2810.000
 band.2.zshift               =              F
 band.2.beam                 =               0
 band.2.calc_point_space     =       0.500E-03
 band.2.wave_factor          =      0.9999999
 band.2.max_opd		     =        180.0000
 band.2.omega	             =		2.4096
 band.2.apodization_code     =               0
 band.2.gasb                 = CH4
 band.2.tempretb =T

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
 band.4.omega                = 	2.3923
 band.4.apodization_code     =               0
 band.4.gasb                 = CH4 
 band.4.tempretb = T
 out.level = 1
 out.gas_spectra = T
