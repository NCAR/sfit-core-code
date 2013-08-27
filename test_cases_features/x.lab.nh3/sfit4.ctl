 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = nh3_1.dat
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 00795.893004-01204.106996.hbin

 
 # Definition for retrieval gases
 
 gas.layers                  =              1

 gas.column.list                         = NH3 O3
 gas.profile.list                         = 
 

 gas.column.NH3.logstate             =               F
 gas.column.NH3.scale             =              1.0
 gas.column.NH3.sigma               = 1.0 

 gas.column.O3.logstate             =               F
 gas.column.O3.scale             =              1.0
 gas.column.O3.sigma               = 1.0
 
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.emission.T_infinity      =		 6000.0
 fw.emission.object	     =  	    .e.
 fw.emission.normalized	     =	 	     T
 fw.isotope_separation       =               F
 fw.lab = T
 fw.lab.length = 0.01
# fw.lab.pressure = 13.3322
# fw.lab.pressure = 9.33254 
 fw.lab.pressure =  4.66627
 fw.lab.temperature = 292.84995

 # Retrieval parameter

 kb = F
 rt                          =              F
  
 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             =  800.0
# band.1.nu_stop              =  931.4
# band.1.nu_start             =  930.43
 band.1.nu_stop              =  1200.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.0500E-03
 band.1.wave_factor           =           1.000
 band.1.max_opd    = 98.11368
 band.1.omega  = 3.0
 band.1.apodization_code                  =               7
 band.1.gasb                 = NH3 O3
 
 out.level = 1
 out.gas_spectra = F
 out.gas_spectra.type = 1
