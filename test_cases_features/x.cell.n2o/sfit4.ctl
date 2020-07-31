 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = CH4_AIR.ref
 file.in.isotope                = isotope.input
 file.in.modulation_fcn         = apod.dat
 file.in.phase_fcn         = phase.dat
 file.in.spectrum               = br_n2o25_20190730.cell
 file.in.linelist               = 02162.989156-02228.755844.hbin
 # Definition for retrieval gases
 
 gas.layers                  =              0
 gas.column.list                         = N2O
 gas.profile.list                         = 
 gas.column.N2O.logstate = F
 gas.column.N2O.scale               = 1.0
 gas.column.N2O.sigma               = 0.1
 gas.column.CH42.logstate = F
 gas.column.CH42.scale               = 1.0
 gas.column.CH42.sigma               = 1.0
 gas.column.CH43.logstate = F
 gas.column.CH43.scale               = 1.0
 gas.column.CH43.sigma               = 1.0
 gas.column.HBR.logstate = F
 gas.column.HBR.scale               = 1.0
 gas.column.HBR.sigma               = 1.0
 gas.column.CO2.logstate = F
 gas.column.CO2.scale               = 1.0
 gas.column.CO2.sigma               = 1.0
 gas.column.H2O.logstate = F
 gas.column.H2O.scale               = 1.0
 gas.column.H2O.sigma               = 1.0
   
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.lshapemodel.sdv          =           T 
 fw.linemixing = T
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               T
 fw.apod_fcn.type = 2
 fw.apod_fcn.order = 9
 fw.phase_fcn                =               T
 fw.phase_fcn.type           =               2	
 fw.phase_fcn.order                =         9
 fw.emission                 =               F
 fw.isotope_separation       =              F 

 # CELL parameter
 cell = 1
 cell.1.temperature = 308.9
 cell.1.pressure = 0.894
 cell.1.gas = N2O
 cell.1.vmr = 1.0
 cell.1.path = 2.2


# Retrieval parameter
 rt                          =            T
 rt.lm                               = T
 rt.lm.gamma_start                   = 1.0e5
 rt.lm.gamma_inc                     = 10.0
 rt.lm.gamma_dec                     = 10.0
 rt.convergence                      = 0.1
 rt.max_iteration                    = 30
 rt.dwshift                          = F
 rt.wshift                           = T
 rt.wshift.type                      = 3
 rt.wshift.apriori                   = 3e-6
 rt.wshift.sigma                     = 0.100
 rt.slope                            = T
 rt.slope.apriori                    = 0.000
 rt.slope.sigma                      = 1.00
 rt.curvature                        = T
 rt.curvature.apriori                = 0.000
 rt.curvature.sigma                  = 1.00
 rt.phase                            = F
 rt.phase.apriori                    = 0.000
 rt.phase.sigma                      = 0.1
 rt.apod_fcn                         = T
 rt.apod_fcn.apriori                 = 1.000
 rt.apod_fcn.sigma                   = 0.100
 rt.phase_fcn                        = T
 rt.phase_fcn.apriori                = 1.000
 rt.phase_fcn.sigma                  = 1.0
 rt.solshift                         = F
 rt.solstrnth                        = F
 rt.temperature                      = F


 # Microwindows and their parameters
 
 band                        =   1 2 3
 band.1.nu_start             =        2167.03
 band.1.nu_stop              =        2185.25
 band.1.zshift               =              F
 band.1.beam                 = 1 2            
 band.1.beam.model = IP
 band.1.beam.1.apriori = 0.001 0.93 2173 0.0
 band.1.beam.1.sigma = 1.0 0.0 1.0 0.0	
 band.1.beam.2.apriori = 0.001 0.24 2173 0.0
 band.1.beam.2.sigma = 1.0 0.0 1.0 0.0	
 band.1.calc_point_space     =       0.100E-03
 band.1.wave_factor          =      1.0	
# band.1.wave_factor          =      1.000000455
 band.1.max_opd		     =        257.14285
 band.1.omega	             =		2.3923
 band.1.apodization_code     =               0
 band.1.gasb                 = N2O
 band.1.tempretb =F
 band.2.nu_start             =        2222.825
 band.2.nu_start             =        2222.825
 band.2.nu_stop              =        2223.019	
 band.2.zshift               =              F
 band.2.beam                 = 1 2            
 band.2.beam.model = IP
 band.2.beam.1.apriori = 0.001 0.93 2223 0.0
 band.2.beam.1.sigma = 1.0 0.0 1.0 0.0	
 band.2.beam.2.apriori = 0.001 0.24 2223 0.0
 band.2.beam.2.sigma = 1.0 0.0 1.0 0.0	
 band.2.calc_point_space     =       0.1e-3
 band.2.wave_factor          =      1.0	
# band.1.wave_factor          =      1.000000455
 band.2.max_opd		     =        257.14285
 band.2.omega	             =		2.3923
 band.2.apodization_code     =               0
 band.2.gasb                 = N2O
 band.2.tempretb =F
 band.3.nu_start             =        2224.457
 band.3.nu_stop              =        2224.715
 band.3.zshift               =              F
 band.3.beam                 = 1 2            
 band.3.beam.model = IP
 band.3.beam.1.apriori = 0.001 0.93 2173 0.0
 band.3.beam.1.sigma = 1.0 0.0 1.0 0.0	
 band.3.beam.2.apriori = 0.001 0.24 2173 0.0
 band.3.beam.2.sigma = 1.0 0.0 1.0 0.0	
 band.3.calc_point_space     =       0.100E-03
 band.3.wave_factor          =      1.0	
# band.1.wave_factor          =      1.000000455
 band.3.max_opd		     =        257.14285
 band.3.omega	             =		2.3923
 band.3.apodization_code     =               0
 band.3.gasb                 = N2O
 band.3.tempretb =F


 out.level = 1
 out.gas_spectra = T
 out.g_matrix = T