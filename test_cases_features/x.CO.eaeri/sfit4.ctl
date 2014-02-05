
 # General

 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 01985.500000-02214.500000.hbin

 # Definition for retrieval gases

 gas.layers                  =              41
 gas.isoflag                 =               F
 gas.profile.list            = CO
 gas.column.list            = O3 N2O H2O OCS 
 gas.profile.CO.correlation               =              F
 gas.profile.CO.scale               =              1.0
 gas.profile.CO.sigma               = 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1
 gas.column.O3.scale               =              1.0
 gas.column.O3.sigma               = 0.1 
 gas.column.N2O.scale               =             1.0
 gas.column.N2O.sigma               = 1.0 
 gas.column.H2O.scale              =          1.0
 gas.column.H2O.sigma               = 1.0 
 gas.column.OCS.scale              =          1.0
 gas.column.OCS.sigma               = 1.0 


 # Forward model parameters

 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.pressure_shift           = T
 fw.emission                 =               T
 fw.emission.T_infinity      =		   2.7
 fw.emission.object	     =  	    .e.
 fw.emission.normalized	     =	 	     F
 # Retrieval parameter

 rt                          =              T 
 rt.lm                       =               T
 rt.lm.gamma_start           =           1.0e5
 rt.lm.gamma_inc             =           1.0e1
 rt.lm.gamma_dec             =           1.0e1
 rt.convergence              =           0.1
 rt.max_iteration            =              17
 rt.wshift                   =               T
 rt.wshift.type              =               3
 rt.wshift.apriori           =           0.000
 rt.wshift.sigma             =           0.100
 rt.slope                    =               T
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           0.100
 rt.continuum                = F
 rt.continuum.type = 1
 rt.continuum.strength       = 0.0
 rt.continuum.tilt       = 0.0
 rt.temperature              =               F
 rt.temperature.sigma =
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 



 # Microwindows and their parameters

 band                        =  1   
 band.1.nu_start             =          2000.0
 band.1.nu_stop              =          2200.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =            0.04
 band.1.wave_factor          =           1.000
 band.1.max_opd                 =            1.00
 band.1.omega                =              45.0
 band.1.apodization_code        =               0
 band.1.gasb                 = CO O3 N2O H2O OCS
 band.1.tempretb = F 



 out.level = 1
 out.gas_spectra = T
 out.pbpfile = T
