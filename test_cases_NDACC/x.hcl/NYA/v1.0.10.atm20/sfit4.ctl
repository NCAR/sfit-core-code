 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.input
 file.in.solarlines             = solar.dat
 file.in.linelist               = 02723.671422-02930.058578.hbin
 

 file.out.pbpfile                    = pbpfile
 file.out.statevec                   = statevec
 file.out.k_matrix                   = k.output
 file.out.kb_matrix                  = kb.output
 file.out.sa_matrix                  = sa.complete
 file.out.retprofiles                = rprfs.table
 file.out.aprprofiles                = aprfs.table
 file.out.ab_matrix                  = ab.output
 file.out.summary                    = summary
 file.out.parm_vectors               = parm.output
 file.out.seinv_vector               = seinv.output
 file.out.sainv_matrix               = sainv.complete
 file.out.smeas_matrix               = smeas.target
 file.out.shat_matrix                = shat.complete
 file.out.ssmooth_matrix             = ssmooth.target
 file.out.ak_matrix                  = ak.target
 file.out.solarspectrum              = solar.output
 file.out.g_matrix		             = d.complete

 # Definition for retrieval gases
 
 
 gas.layers                  =              41
 gas.profile.list                 = HCL 
 gas.column.list                  = CH4 NO2 O3 N2O HDO
 gas.profile.HCL.correlation        =               F
 gas.profile.HCL.correlation.type        =               2
 gas.profile.HCL.correlation.width        =               4.0
 gas.profile.HCL.correlation.minalt        =               0.0
 gas.profile.HCL.correlation.maxalt        =               120.0
 gas.profile.HCL.logstate             =               F
 gas.profile.HCL.scale = 1.0		
 gas.profile.HCL.sigma                =
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1
 gas.column.HDO.logstate            =               F
 gas.column.HDO.scale               = 1.0
 gas.column.HDO.sigma               = 1.0
 gas.column.H2O.logstate            =               F
 gas.column.CH4.logstate            =               F
 gas.column.CH4.scale               =             1.0
 gas.column.CH4.sigma               = 1.0
 gas.column.N2O.logstate            =               F
 gas.column.N2O.scale               =             1.0
 gas.column.N2O.sigma               = 1.0
 gas.column.O3.logstate            =               F
 gas.column.O3.scale = 1.0
 gas.column.O3.sigma               = 1.0
 gas.column.NO2.logstate            =               F
 gas.column.NO2.scale               =             1.0
 gas.column.NO2.sigma               = 1.0
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum           =               T
 fw.pressure_shift          =               T
 fw.isotope_separation                    =               T
 fw.tips               =F
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =               F
 rt.lm.gamma_start	     =           1.0e3
 rt.lm.gamma_inc             =            10.0
 rt.lm.gamma_dec             =            10.0
 rt.convergence              =             0.1
 rt.max_iteration            =              17
 rt.wshift = T
 rt.wshift.type                   =               3
 rt.wshift.apriori           =          4.8316E-06
 rt.wshift.sigma                =           0.100
 rt.slope                    = T
 rt.slope.apriori            =           0.000
 rt.slope.sigma                 =           0.100
 rt.curvature                = F
 rt.curvature.apriori        =           0.000
 rt.curvature.sigma             =           0.100
 rt.dwshift                   =               F
 rt.temperature              =               F
 rt.temperature.sigma                =
 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01
 0.01
 
 # Microwindows and their parameters
 
 band                        = 1 2 3
 band.1.nu_start             =        2727.7300
 band.1.nu_stop              =        2727.8300
 band.1.zshift               =               F
 band.1.zshift.apriori       =             0.0
 band.1.zshift.sigma            =             0.2
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.900E-03
 band.1.wave_factor               =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.3932
 band.1.apodization_code     =               0
 band.1.gasb                 = HCL O3 HDO
 band.1.tempretb = F
 band.2.nu_start             =        2775.7000
 band.2.nu_stop              =        2775.8000
 band.2.zshift               =               F
 band.2.zshift.apriori       =             0.0
 band.2.zshift.sa            =             0.2
 band.2.beam                 =               0
 band.2.calc_point_space     =       0.900E-03
 band.2.wave_factor          =           1.000
 band.2.max_opd= 180.0000
 band.2.omega= 2.3932
 band.2.apodization_code     =               0
 band.2.gasb                 = HCL O3 N2O
 band.3.nu_start             =        2925.800
 band.3.nu_stop              =        2926.00
 band.3.zshift               =               F
 band.3.zshift.apriori       =             0.0
 band.3.zshift.sa            =             0.2
 band.3.beam                 =               0
 band.3.calc_point_space     =       0.900E-03
 band.3.wave_factor          =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 2.3932
 band.3.apodization_code     =               0
 band.3.gasb                 = HCL CH4 NO2 
 
 kb = T
 kb.temperature = T

 out.level                           = 1
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.g_matrix			                = T
 out.k_matrix                        = T
 out.shat_matrix                     = T
 out.retprofiles                     = T
 out.aprprofiles                     = T
 out.ab_matrix                       = F
 out.ak_matrix                       = T
 out.summary                         = T
 out.pbpfile                         = T
 out.channel                         = F
 out.parm_vectors                    = T
 out.seinv_vector                    = T
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = T
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F

