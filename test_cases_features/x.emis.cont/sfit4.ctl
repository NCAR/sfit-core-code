
 # General

 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 00395.034545-00804.965455.hbin

 # Definition for retrieval gases

 gas.layers                  =              41
 gas.profile.list            = 
 gas.column.list            = CO2  O3 H2O CH4 N2O 
 gas.profile.O3.correlation               =              F
 gas.profile.O3.scale               =              1.0
 gas.profile.O3.sigma               = 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1
 gas.column.O3.scale               =              1.0
 gas.column.O3.sigma               = 0.1 
 gas.column.CO2.scale               =             1.0
 gas.column.CO2.sigma               = 1.0 
 gas.column.H2O.scale              =          1.0
 gas.column.H2O.sigma               = 1.0 
 gas.column.CH4.scale              =          1.0
 gas.column.CH4.sigma               = 1.0 
 gas.column.N2O.scale              =          1.0
 gas.column.N2O.sigma               = 1.0 
 gas.column.CO.scale              =          1.0
 gas.column.CO.sigma               = 1.0 


 # Forward model parameters

 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.pressure_shift           =               T
 fw.emission                 =               T
 fw.emission.T_infinity      =		   2.7
 fw.emission.object	     =  	    .e.
 fw.emission.normalized	     =	 	     F
 fw.continuum                = T
 fw.continuum.type           = 2
 fw.continuum.strength       = 1
 fw.continuum.z              = 0.5
 fw.tips = F
 # Retrieval parameter

 rt                          =               F
 rt.lm                       =               T
 rt.lm.gamma_start           =           1.0e3
 rt.lm.gamma_inc             =           1.0e1
 rt.lm.gamma_dec             =           1.0e1
 rt.convergence              =           0.1
 rt.max_iteration            =              17
 rt.wshift                   =               F
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
 band.1.nu_start             =          400.0
 band.1.nu_stop              =          800.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =            0.04
 band.1.wave_factor          =           1.000
 band.1.max_opd                 =         11.0
 band.1.omega                =            26.0 
 band.1.apodization_code        =               0
 band.1.gasb                 = CO2 O3 N2O CH4 H2O
 band.1.tempretb = F   

 band.2.nu_start             =           750.0
 band.2.nu_stop              =           800.0
 band.2.zshift               =               1
 band.2.beam                 =               0
 band.2.dn                   =            0.4
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
