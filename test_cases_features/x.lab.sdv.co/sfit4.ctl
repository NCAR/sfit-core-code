 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = B0051.2b.asc
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 04095.941672-04204.058328.hbin
 

 
 # Definition for retrieval gases
 
 gas.layers                  =              1

 gas.column.list                         = CO
 gas.profile.list                         = 
 

 gas.column.CO.logstate             =               F
 gas.column.CO.scale             =              1.0
 gas.column.CO.sigma               = 1.0 

 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.linemixing = T
 fw.linemixing.gas = CO
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.emission.T_infinity      =		 6000.0
 fw.emission.object	     =  	    .e.
 fw.emission.normalized	     =	 	     F
 fw.isotope_separation       =               F
 fw.lab = T
 fw.lab.length = 0.4076
 fw.lab.pressure = 943.46
 fw.lab.temperature = 298.0

 # Retrieval parameter

 kb = F
 rt                          =              T
 rt.max_iteration = 15
 rt.convergence = 0.1
 rt.wshift = T
 rt.wshift.type = 3
 rt.wshift.apriori = 0.0
 rt.wshift.sigma = 0.1  
 rt.slope = T
 rt.slope.apriori = 0.0
 rt.slope.sigma = 0.1
 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             =  4100.0
 band.1.nu_stop              =  4200.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.0500E-03
 band.1.wave_factor           =           0.9999973396
 band.1.max_opd    = 180
 band.1.omega  = 1.2
 band.1.apodization_code                  =               0
 band.1.gasb                 = CO
 
 out.level = 1
 out.gas_spectra = F
 out.gas_spectra.type = 1
