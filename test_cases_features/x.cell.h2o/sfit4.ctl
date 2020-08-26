 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = CH4_AIR.ref
 file.in.isotope                = isotope.input
 file.in.modulation_fcn         = apod.dat
 file.in.phase_fcn         = phase.dat
 file.in.spectrum               = art_spectrum
 file.in.linelist               = 00494.950000-01405.050000.hbin
 # Definition for retrieval gases
 
 gas.layers                  =              0
 gas.column.list                         = 
 gas.profile.list                         = H2O 
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
 gas.profile.H2O.logstate = F
 gas.profile.H2O.scale               = 1.0
 gas.profile.H2O.sigma               = 1.0
   
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.lshapemodel.sdv          =           T 
 fw.linemixing = T
 fw.solar_spectrum	     =               F
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.isotope_separation       =              F
 fw.tips = F
 fw.mtckd_continuum = T

 # CELL parameter
 cell = 1
 cell.1.temperature = 293.0
 cell.1.pressure = 1000.0
 cell.1.gas = H2O
 cell.1.vmr = 1.0
 cell.1.path = 10.0


# Retrieval parameter
 rt                          =            T
 rt.lm                               = T
 rt.lm.gamma_start                   = 1.0e5
 rt.lm.gamma_inc                     = 10.0
 rt.lm.gamma_dec                     = 10.0
 rt.convergence                      = 0.1
 rt.max_iteration                    = 30
 rt.dwshift                          = F
 rt.wshift                           = F
 rt.wshift.type                      = 3
 rt.wshift.apriori                   = 3e-6
 rt.wshift.sigma                     = 0.100
 rt.slope                            = F
 rt.slope.apriori                    = 0.000
 rt.slope.sigma                      = 1.00
 rt.curvature                        = F
 rt.curvature.apriori                = 0.000
 rt.curvature.sigma                  = 1.00
 rt.phase                            = F
 rt.phase.apriori                    = 0.000
 rt.phase.sigma                      = 0.1
 rt.apod_fcn                         = F
 rt.phase_fcn                        = F
 rt.phase_fcn.apriori                = 1.000
 rt.phase_fcn.sigma                  = 1.0
 rt.solshift                         = F
 rt.solstrnth                        = F
 rt.temperature                      = F


 # Microwindows and their parameters
 
 band                        =   1 
 band.1.nu_start             =        500.0
 band.1.nu_stop              =        1400.0
 band.1.zshift               =              F
 band.1.beam                 =      
 band.1.beam.model = IP
 band.1.beam.1.apriori = 0.001 0.93 2173 0.0
 band.1.beam.1.sigma = 1.0 0.0 1.0 0.0	
 band.1.beam.2.apriori = 0.001 0.24 2173 0.0
 band.1.beam.2.sigma = 1.0 0.0 1.0 0.0	
 band.1.calc_point_space     =       0.100
 band.1.wave_factor          =      1.0	
# band.1.wave_factor          =      1.000000455
 band.1.max_opd		     =        10.0
 band.1.omega	             =		2.3923
 band.1.apodization_code     =               0
 band.1.gasb                 = H2O
 band.1.tempretb =F

 out.level = 1
 out.gas_spectra = T
 out.g_matrix = T