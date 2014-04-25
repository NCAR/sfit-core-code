 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = ch4_1.dat
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 04226.944346-04236.055654.hbin

 
 # Definition for retrieval gases
 
 gas.layers                  =              1

 gas.column.list                         = CH4
 gas.profile.list                         = 
 

 gas.column.CH4.logstate             =               F
 gas.column.CH4.scale             =              1.0
 gas.column.CH4.sigma               = 1.0 

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
 fw.lab.length = 0.2
 fw.lab.pressure =  93.325
 fw.lab.temperature = 297.55

 # Retrieval parameter

 kb = F
 rt                          =              T
 rt.max_iteration = 20
 rt.convergence = 0.1
 rt.wshift = T
 rt.wshift.type = 3
 rt.wshift.apriori = 0.0
 rt.wshift.sigma = 0.1
  
 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             =  4231.0
 band.1.nu_stop              =  4232.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor           =           1.000
 band.1.max_opd    = 188.378
 band.1.omega  = 2.2
 band.1.apodization_code                  =               7
 band.1.gasb                 = CH4
 
 out.level = 1
 out.gas_spectra = F
 out.gas_spectra.type = 1
