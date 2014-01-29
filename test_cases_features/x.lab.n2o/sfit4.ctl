 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = 131216_46m4_RT_1000mb2_av.ref
 file.in.spectrum               = 131216_46m4_RT_1000mb2_av
 file.in.solarlines             = /home/mathias/linelist/solar/120621/solar.dat
 file.in.linelist               = 02155.883444-02264.116556.hbin
 file.out.ak_matrix = ak.out
 
 # Definition for retrieval gases
 
 gas.layers                  =              1

 gas.column.list                         = N2O CO CO2 H2O
 gas.profile.list                         = 
 

 gas.column.N2O.logstate             =               F
 gas.column.N2O.scale             =              1.0
 gas.column.N2O.sigma               = 1.0 

 gas.column.CO.logstate             =               F
 gas.column.CO.scale             =              1.0
 gas.column.CO.sigma               = 1.0

 gas.column.H2O.logstate             =               F
 gas.column.H2O.scale             =              1.0
 gas.column.H2O.sigma               = 1.0

 gas.column.CO2.logstate             =               F
 gas.column.CO2.scale             =              1.0
 gas.column.CO2.sigma               = 1.0
 
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.linemixing = T
 fw.linemixing.gas = N2O
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.isotope_separation       =               F
 fw.lab = T
 fw.lab.length = 46.4
# fw.lab.pressure =  498.2
 fw.lab.pressure =  1000.2
 fw.lab.temperature = 293.922

 # Retrieval parameter

 kb = F
 rt                          =              T
 rt.max_iteration = 20
 rt.convergence = 0.1
 rt.wshift = F
 rt.wshift.type = 3
 rt.wshift.apriori = 0.0
 rt.wshift.sigma = 0.1
 rt.slope = T
 rt.slope.apriori = 0.0
 rt.slope.sigma = 0.1
 rt.curvature = T
 rt.curvature.apriori = 0.0
 rt.curvature.sigma = 0.1
  
 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             =  2160.0
 band.1.nu_stop              =  2260.0
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.500E-03
 band.1.wave_factor           =           1.000
 band.1.max_opd    = 90
 band.1.omega  = 2.75
 band.1.apodization_code     =  0
 band.1.gasb                 = N2O CO2 CO H2O
 
 out.level = 1
 out.gas_spectra = T
 out.gas_spectra.type = 1
 out.ak_matrix = T
