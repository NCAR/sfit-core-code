
 # General

 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 00985.500000-01114.500000.hbin

 # Definition for retrieval gases

 gas.layers                  =              41
 gas.profile.list            = O3 
 gas.column.list            = CO2 H2O
 gas.profile.O3.regmethod   = 'OEM'
 gas.profile.O3.regmethod.lambda   = 0.001	
 gas.profile.O3.correlation               =              F
 gas.profile.O3.scale               =              1.0
 gas.profile.O3.sigma               = 
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 
 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 
 0.5
 gas.profile.H2O.regmethod   = 'TP'
 gas.profile.H2O.regmethod.lambda   = 0.001	
 gas.profile.H2O.correlation               =              F
 gas.profile.H2O.scale               =              1.0
 gas.profile.H2O.sigma               = 
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
 1.0
 gas.column.O3.scale               =              1.0
 gas.column.O3.sigma               = 0.1 
 gas.column.CO2.scale               =             1.0
 gas.column.CO2.sigma               = 1.0 
 gas.column.H2O.scale              =          1.0
 gas.column.H2O.sigma               = 1.0 


 # Forward model parameters

 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.pressure_shift           = T
 fw.emission                 =               T
 fw.emission.T_infinity      =		   2.7
 fw.emission.object	     =  	    .e.
 fw.emission.normalized	     =	 	     F
 fw.mtckd_continuum = T
 fw.continuum = T
 fw.continuum.type = 1
 fw.continuum.strength = 0.01
 fw.continuum.order = 4
 # Retrieval parameter

 rt                          =               T
 rt.lm                       =               F
 rt.lm.gamma_start           =           1.0e2
 rt.lm.gamma_inc             =           1.0e1
 rt.lm.gamma_dec             =           1.0e1
 rt.convergence              =           0.1
 rt.max_iteration            =              17
 rt.wshift                   =               T
 rt.wshift.type              = 3
 rt.wshift.apriori           =           0.000
 rt.wshift.sigma             =           0.100
 rt.slope                    =               F
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           0.100
 rt.temperature              =               F
 rt.continuum                = T
 rt.continuum.sigma          = 0.1
 rt.temperature.sigma =
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 



 # Microwindows and their parameters

 band                        =  1   
 band.1.nu_start             =          1000.0
 band.1.nu_stop              =          1050.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =            0.1
 band.1.wave_factor          =           1.000
 band.1.max_opd                 =            1.00
 band.1.omega                =            45.0 
 band.1.apodization_code        =               0
 band.1.gasb                 = O3     CO2 H2O   
 band.2.tempretb = F   

 band.2.nu_start             =           750.0
 band.2.nu_stop              =           800.0
 band.2.zshift               =               1
 band.2.beam                 =               0
 band.2.dn                   =            0.04
 band.2.wavfac               =           1.000
 band.2.pmax                 =            1.00
 band.2.omega                =             3.86
 band.2.iap                  =               0
 band.2.snr                  =          0.0007
 band.2.gasb                 = O3     CO2 H2O



 out.level = 1
 out.gas_spectra = T
 out.gas_spectra.type = 1
 out.pbpfile = T
 out.sainv_matrix = T
