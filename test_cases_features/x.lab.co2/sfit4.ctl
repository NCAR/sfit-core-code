 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = 961027R0.002.asc
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 00780.821656-00811.178344.hbin

 
 # Definition for retrieval gases
 
 gas.layers                  =              1

 gas.column.list                         = CO2 
 gas.profile.list                         = 
  
 gas.column.CO2.logstate             =               F
 gas.column.CO2.scale             =              1.0
 gas.column.CO2.sigma               = 1.0

 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               1
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               F
 fw.linemixing = F
 fw.linemixing.gas = CO2 
 fw.lab = T
 fw.lab.length = 0.24
 fw.lab.pressure = 13.3322
 fw.lab.temperature = 293.15

 # Retrieval parameter

 kb = F
 rt                          =              F
 rt.max_iteration = 10
 rt.convergence = 0.1
 rt.wshift = T
 rt.wshift.type = 3
 rt.wshift.apriori = 0.0
 rt.wshift.sigma = 0.1  
 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             =         785.000
 band.1.nu_stop              =         807.000
# band.1.nu_start             =  800.0
# band.1.nu_stop              =  931.4
# band.1.nu_start             =  930.43
# band.1.nu_stop              =  1200.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.0500E-03
 band.1.wave_factor           =           1.0001
 band.1.max_opd    = 58.868
 band.1.omega  = 3.0
 band.1.apodization_code                  =               7
 band.1.gasb                 = CO2
 
 out.level = 1
 out.gas_spectra = F
 out.gas_spectra.type = 1
 out.ak_matrix = T
