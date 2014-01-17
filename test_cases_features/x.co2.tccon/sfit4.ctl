
 # General

 file.in.stalayers              = station.layers
 file.in.isotope                = isotope.inp
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = ../../linelist/solar/120621/solar.dat
 file.in.linelist               = 06175.842223-08009.162777.hbin

 # Definition for retrieval gases

 gas.layers                  =              71
 gas.profile.list            =
 gas.column.list             = CO2 H2O CH4 HDO O2 O2CIA
 gas.column.CO2.scale               = 1.0
 gas.column.CO2.sigma               = 1.0
 gas.column.O2.scale                = 1.0
 gas.column.O2.sigma                = 1.0
 gas.column.O2CIA.scale             = 0.00001
 gas.column.O2CIA.sigma             = 1.0
 gas.column.H2O.scale               = 1.0
 gas.column.H2O.sigma               = 1.0
 gas.column.CH4.scale               = 1.0
 gas.column.CH4.sigma               = 1.0
 gas.column.HDO.scale               = 1.0
 gas.column.HDO.sigma               = 1.0

 # Forward model parameters

 fw.delnu                    =             0.1
 fw.lshapemodel              =               0
 fw.solar_spectrum           =               T
 fw.linemixing               =               T
 fw.linemixing.gas           =             CO2
 fw.isotope_separation       =               T
 fw.pressure_shift           =               T
 # Retrieval parameter

 rt                          =               T
 rt.lm                       =               F
 rt.lm.gamma_start           =            1.e5
 rt.lm.gamma_inc             =            10.0
 rt.lm.gamma_dec             =            10.0
 rt.convergence              =             0.1
 rt.max_iteration            =              10
 rt.wshift =                                 T
 rt.wshift.type              =               3
 rt.wshift.apriori           =             0.0
 rt.wshift.sigma                =          0.1
 rt.slope                    =               T
 rt.slope.apriori            =             0.0
 rt.slope.sigma                 =          0.1
 rt.solshift                 =               T
 rt.solshift.apriori         =             0.0
 rt.solshift.sigma              =          0.1
 rt.solstrnth                 =               T
 rt.solstrnth.apriori         =             0.0
 rt.solstrnth.sigma              =          0.1
 rt.temperature              =               F

 # Kb matrices
 kb = T
 kb.profile = T
 kb.profile.gas = CO2 O2

 # Microwindows and their parameters

 band                        =   1 2 3
 band.1.nu_start             =        6180.000
 band.1.nu_stop              =        6260.000
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       5.000E-02
 band.1.wave_factor          =           1.000
 band.1.max_opd  = 64.2924
 band.1.omega= 2.3923
 band.1.apodization_code     =               0
 band.1.gasb                 = CO2 H2O CH4 HDO
 band.2.nu_start             =        6297.000
 band.2.nu_stop              =        6382.000
 band.2.zshift               =               F
 band.2.beam                 =               0
 band.2.calc_point_space     =       5.000E-02
 band.2.wave_factor          =           1.000
 band.2.max_opd = 64.2924
 band.2.omega= 2.3923
 band.2.apodization_code     =               0
 band.2.gasb                 = CO2 H2O CH4 HDO
 band.3.nu_start             =        7765.000
 band.3.nu_stop              =        8005.000
 band.3.zshift               =               F
 band.3.beam                 =               0
 band.3.calc_point_space     =       5.000E-03
 band.3.wave_factor          =           1.000
 band.3.max_opd= 64.2924
 band.3.omega= 2.3923
 band.3.apodization_code    =               0
 band.3.gasb                 = O2 O2CIA H2O


# Output Files Section

 out.level                           = 1
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.k_matrix                        = T
 out.shat_matrix                     = T
 out.retprofiles                     = T
 out.aprprofiles                     = T
 out.ab_matrix                       = F
 out.ak_matrix                       = F
 out.summary                         = T
 out.pbpfile                         = T
 out.channel                         = F
 out.parm_vectors                    = T
 out.seinv_vector                    = F
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = F
 out.raytrace                        = F
 out.raytrace.type                   = 1
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F
