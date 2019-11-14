 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = ../../linelist/solar/120621/solar.dat
 file.in.linelist               = 02031.021722-02039.778278.hbin
 
 # Definition for retrieval gases
 
 gas.layers                  =              41
 gas.profile.list            = CO2 O3
 gas.column.list             = H2O OCS  
 gas.profile.CO2.correlation          =               T
 gas.profile.CO2.correlation.type     =               1
 gas.profile.CO2.correlation.minalt  =               0.0
 gas.profile.CO2.correlation.maxalt  =               150.0
 gas.profile.CO2.correlation.width    =               1000.0
 gas.profile.CO2.logstate             =               F
 gas.profile.CO2.scale                =               1.0
 gas.profile.CO2.sigma                =
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
 gas.profile.O3.correlation.type     =               1
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
 gas.column.OCS.scale               = 1.0
 gas.column.OCS.sigma               = 1.0
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
 
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =               F
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
 
 # Microwindows and their parameters
 
 band                        =   1
 band.1.nu_start             =        2035.08
 band.1.nu_stop              =        2035.72
 band.1.zshift		     =	      F
 band.1.beam                 =               
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor          =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.3923
 band.1.apodization_code =               0
 band.1.gasb                 = CO2 H2O    O3     OCS    
 band.2.nu_start             =        809.16
 band.2.nu_stop              =        809.56
 band.2.zshift		     =	      F
 band.2.beam                 =               
 band.2.calc_point_space                   =       0.500E-03
 band.2.wave_factor              =           1.000
 band.2.max_opd= 180.0000
 band.2.omega= 2.3923
 band.2.apodization_code                  =               0
 band.2.gasb                 = CO2 O3 H2O
 band.3.nu_start             =        968.9
 band.3.nu_stop              =        969.3
 band.3.zshift		     =	      F
 band.3.beam                 =  1
 band.3.beam.model = PS
 band.3.beam.1.apriori       = 0.01 0.22 969.0 0.0
 band.3.beam.1.sigma         = 0.1 0.01 1.0 0.0
 band.3.calc_point_space                   =       0.500E-03
 band.3.wave_factor              =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 2.3923
 band.3.apodization_code                  =               0
 band.3.gasb                 = CO2 O3 


 out.level= 1
 out.smeas_matrix = T
 
 out.gas_spectra= T
 out.k_matrix= T
 out.ak_matrix= T
 out.g_matrix= T
 out.sa_matrix= T
 out.retprofiles= T
 out.aprprofiles= T
 file.out.ak_matrix= ak.out
 file.out.k_matrix= k.out
 file.out.g_matrix= g.out
 file.out.kb_matrix= kb.out
 file.out.sa_matrix= sa.out
 file.out.seinv_vector = seinv.out
 file.out.retprofiles= rprfs.table
 file.out.aprprofiles= aprfs.table
 file.out.summary= summary
 kb= F
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
