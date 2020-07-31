 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.isotope                = isotope.inp
 file.in.spectrum               = B0052.2b.asc
 file.in.solarlines             = solar.dat
 file.in.linelist               = 04095.941722-04404.058278.hbin
 

 
 # Definition for retrieval gases
 
 gas.layers                  =              0

 gas.column.list             = COO1 COO2 #COO3 COO4
 gas.profile.list                         = 
 

 gas.column.CO.logstate             =               F
 gas.column.CO.scale             =              1.0
 gas.column.CO.sigma               = 1.0 

 gas.column.COO1.logstate             =               F
 gas.column.COO1.scale             =              1.0
 gas.column.COO1.sigma               = 1.0 

 gas.column.COO2.logstate             =               F
 gas.column.COO2.scale             =              1.0
 gas.column.COO2.sigma               = 1.0 

 gas.column.COO3.logstate             =               F
 gas.column.COO3.scale             =              1.0
 gas.column.COO3.sigma               = 1.0 

 gas.column.COO4.logstate             =               F
 gas.column.COO4.scale             =              1.0
 gas.column.COO4.sigma               = 1.0 

 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =   4            
 fw.lshapemodel.sdv          = T
 fw.linemixing = F
 fw.linemixing.gas = CO # COO1 COO2 COO3 COO4
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.emission.T_infinity      =		 6000.0
 fw.emission.object	     =  	    .e.
 fw.emission.normalized	     =	 	     F
 fw.isotope_separation       =               T

 cell = 1 2
 cell.1.temperature = 151.2
 cell.1.pressure = 133.1487
 cell.1.gas = COO1
 cell.1.vmr = 0.0702
 cell.1.path = 40.760
 cell.2.temperature = 151.2
 cell.2.pressure = 133.1487
 cell.2.gas = COO2
 cell.2.vmr = 0.0702
 cell.2.path = 40.760


 # Retrieval parameter

 kb = F
 rt                          =            T
 rt.max_iteration = 15
 rt.convergence = 0.1
 rt.wshift = T
 rt.wshift.type = 3
 rt.wshift.apriori = 0.0
 rt.wshift.sigma = 0.1  
 rt.slope = T
 rt.slope.apriori = 0.0
 rt.slope.sigma = 0.1
 rt.temperature = F
 rt.temperature.sigma = 0.01 0.01
 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             = 4100.0
 band.1.nu_stop             = 4400.0
# band.1.nu_start             =  4222.2
# band.1.nu_stop              =  4236.5
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor           =           1.0
 band.1.max_opd    = 180
 band.1.omega  = 1.2
 band.1.apodization_code                  =               0
 band.1.gasb                 = COO1 COO2 
 band.1.tempretb = F

 out.level = 1
 out.gas_spectra = F
 out.gas_spectra.type = 1
