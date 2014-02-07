
 # General

 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = ../../linelist/solar/120621/solar.dat
 file.in.linelist               = 01985.500000-02214.500000.hbin

 # Definition for retrieval gases

 gas.layers                  =              41
 gas.isoflag                 =               F
 gas.profile.list            = CO 
 gas.column.list            =  N2O H2O
 gas.profile.CO.correlation               =              F
 gas.profile.CO.scale               =              1.0
 gas.profile.CO.sigma               = 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1
 gas.profile.H2O.correlation               =              F
 gas.profile.H2O.scale               =              1.0
 gas.profile.H2O.sigma               = 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
 0.1
 gas.column.CO.scale               =              1.0
 gas.column.CO.sigma               = 0.1 
 gas.column.N2O.scale               =             1.0
 gas.column.N2O.sigma               = 1.0 
 gas.column.H2O.scale              =          1.0
 gas.column.H2O.sigma               = 1.0 
 gas.column.OCS.scale              =          1.0
 gas.column.OCS.sigma               = 1.0 


 # Forward model parameters

 fw.delnu                    =           0.10000
 fw.lshapemodel              =              1 
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
 rt.convergence              =           0.2
 rt.max_iteration            =              17
 rt.wshift                   =               F
 rt.slope                    =               F
 rt.continuum                = T 
 rt.continuum.order = 3 
 rt.continuum.apriori       = 0.0 
 rt.continuum.sigma       = 0.1 
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
 band.1.calc_point_space     =            0.1
 band.1.wave_factor          =           1.000
 band.1.max_opd                 =            1.00
 band.1.omega                =              45.0
 band.1.apodization_code        =               0
 band.1.gasb                 = CO N2O H2O 
 band.1.tempretb = F



 out.level = 1
 out.gas_spectra = T
 out.pbpfile = T
