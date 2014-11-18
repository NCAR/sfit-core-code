
 # General

 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 00854.834444-00905.165556.hbin

 # Definition for retrieval gases

 gas.layers                  =              41
 gas.isoflag                 =               F
 gas.profile.list            = HNO3 H2O
 gas.column.list            =  
 gas.profile.HNO3.correlation               =              F
 gas.profile.HNO3.scale               =              1.0
 gas.profile.HNO3.sigma               = 
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
 gas.column.O3.scale               =              1.0
 gas.column.O3.sigma               = 1.0 
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
 # Retrieval parameter

 rt                          =               T
 rt.lm                       =               F
 rt.lm.gamma_start           =           1.0e5
 rt.lm.gamma_inc             =           1.0e1
 rt.lm.gamma_dec             =           1.0e1
 rt.convergence              =           0.1
 rt.max_iteration            =              17
 rt.wshift                   =               T
 rt.wshift.type              = 3
 rt.wshift.apriori           =           0.0
 rt.wshift.sigma             =           0.100
 rt.slope                    =               F
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           0.100
 rt.continuum                =T 
 rt.continuum.order = 0
 rt.continuum.apriori = 0.0
 rt.continuum.sigma       = 1.0
 rt.temperature = F
 rt.temperature.sigma =
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 



 # Microwindows and their parameters

 band                        =  1   
 band.1.nu_start             =          860.0
 band.1.nu_stop              =          900.0
 band.1.zshift               =              F
 band.1.zshift.apriori       = 0.0
 band.1.zshift.sigma         = 1.0
 band.1.zshift.type = 1 
 band.1.beam                 =               0
 band.1.calc_point_space     =            0.01
 band.1.wave_factor          =           0.9999999860332405
 band.1.max_opd                 =          9.0
 band.1.omega                =            6.3636
 band.1.apodization_code        =               0
 band.1.gasb                 = HNO3    H2O
 band.1.tempretb = F




 out.level = 2
 out.gas_spectra = T
 out.gas_spectra.type = 2
 out.pbpfile = T