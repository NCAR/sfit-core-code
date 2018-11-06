
 # General

 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = /usr/local/share/linelist_v0.9/solar/120621/solar.dat
 file.in.linelist               = 00780.871900-00811.128100.hbin

 # Definition for retrieval gases

 gas.layers                  =              41
 gas.profile.list            = CCL4 H2O O3  
 gas.column.list             = CO2
 gas.profile.CCL4.correlation              =               T
 gas.profile.CCL4.correlation.type              =               2
 gas.profile.CCL4.correlation.width       =          4.000
 gas.profile.CCL4.correlation.minalt               =           0.000
 gas.profile.CCL4.correlation.maxalt               =         120.000
 gas.profile.CCL4.logstate           =               F
 gas.profile.CCL4.scale           =               1.0
 gas.profile.CCL4.sigma              =
   0.02580   0.03520   0.03730   0.03960   0.04220
   0.04500   0.04820   0.05170   0.05560   0.05980
   0.06420   0.06880   0.07330   0.07720   0.08030
   0.08220   0.08330   0.08420   0.08540   0.08640
   0.08770   0.08910   0.09020   0.09170   0.09280
   0.09450   0.09620   0.09760   0.09950   0.10100
   0.10310   0.10550   0.10730   0.10990   0.11190
   0.11490   0.11810   0.12070   0.12440   0.12440
   0.12440
 gas.profile.H2O.correlation               =               F
 gas.profile.H2O.logstate            =               F
 gas.profile.H2O.scale               = 1.0
 gas.profile.H2O.sigma               =
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000
 gas.profile.O3.correlation                =               F
 gas.profile.O3.scale                = 1.0
 gas.profile.O3.sigma                = 
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000   1.00000   1.00000   1.00000   1.00000
   1.00000
 gas.column.CO2.scale               = 1.0
 gas.column.CO2.sigma               = 1.0


 # Forward model parameters

 fw.delnu                    =           0.10000
 fw.lshapemodel              =               4
 fw.solar_spectrum           =               T
 fw.linemixing               =              T
 fw.linemixing.gas           =            CO2
 fw.pressure_shift           =    T 

 # Retrieval parameter

 rt                          =              T
 rt.lm                       =               F
 rt.convergence              =           0.1
 rt.max_iteration            =              15
 rt.wshift                   =               T
 rt.wshift.type              = 3
 rt.wshift.apriori           =           0.000
 rt.wshift.sigma             =           0.010
 rt.slope                    =               T
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           1.000
 rt.curvature                =              F 
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma          =           1.000
 rt.dwshift                  =   F
 rt.temperature              =               F
 rt.solshift                 = T
 rt.solshift.apriori         = 0.000
 rt.solshift.sigma           = 1.00
 rt.solstrnth                 = T
 rt.solstrnth.apriori         = 0.000
 rt.solstrnth.sigma           = 1.00

 # Microwindows and their parameters

 band                        =   1
 band.1.nu_start             =         785.000
 band.1.nu_stop              =         807.000
 band.1.zshift               =               F
 band.1.zshift.type = 1
 band.1.zshift.apriori       =           0.000
 band.1.zshift.sigma         =           1.000
 band.1.beam                 =               0
 band.1.calc_point_space     =       1.00E-03
 band.1.wave_factor          =           1.000
 band.1.max_opd              =          81.967
 band.1.omega                =           5.981
 band.1.apodization_code     =               0
 band.1.gasb                 = CCL4   H2O    O3     CO2

 kb = F
 kb.temperature = T
 kb.slope = T
 kb.curvature = T 
 kb.solshft = T
 kb.solstrnth = T 
 kb.phase = T
 kb.wshift = T
 kb.apod_fcn = T
 kb.phase_fcn = T 
 kb.zshift = T 
 kb.sza = T
 kb.omega = T 
 kb.line = T
 kb.line.type = 1 
 kb.line.gas = CCL4
 kb.profile = T
 kb.profile.gas = CO2 

 out.level = 1
 out.g_matrix = T
