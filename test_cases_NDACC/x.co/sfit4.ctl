 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.solarlines             = /home/mathias/sfit-core-code/linelist/solar/120621/solar.dat
 file.in.linelist               = 02053.641622-02163.208378.hbin
 file.out.ak_matrix= ak.out
 
 # Definition for retrieval gases
 
 gas.layers                  =              43
 gas.isoflag                 =               F
 gas.profile.list            = CO O3 N2O
 gas.column.list             = CO2 OCS H2O
 
 gas.profile.CO.correlation  =               F
 gas.profile.CO.logstate     =               F
 gas.profile.CO.scale        =             1.0
 gas.profile.CO.sigma        =
 0.12680   0.11340   0.12250   0.13420   0.15000
 0.16040   0.17320   0.18970   0.19780   0.20700
 0.21760   0.23010   0.24490   0.26310   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390
 
 gas.profile.O3.correlation  =               F
 gas.profile.O3.logstate     =               F
 gas.profile.O3.scale        =             1.0
 gas.profile.O3.sigma        =
 0.12680   0.11340   0.12250   0.13420   0.15000
 0.16040   0.17320   0.18970   0.19780   0.20700
 0.21760   0.23010   0.24490   0.26310   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27390   0.27390   0.27390   0.27390   0.27390
 0.27380   0.27380   0.27380
 
 gas.profile.H2O.correlation  =               F
 gas.profile.H2O.logstate     =               F
 gas.profile.H2O.scale        =             1.0
 gas.profile.H2O.sigma        =
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0
 
 gas.profile.N2O.correlation  =               F
 gas.profile.N2O.logstate     =               F
 gas.profile.N2O.scale        =             1.0
 gas.profile.N2O.sigma        =
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0 1.0 1.0  1.0 1.0 1.0 1.0 1.0
 1.0 1.0 1.0
 
 gas.column.H2O.logstate     =                F
 gas.column.H2O.scale        =              1.0
 gas.column.H2O.sigma        =              1.0
 
 gas.column.N2O.logstate     =                F
 gas.column.N2O.scale        =              1.0
 gas.column.N2O.sigma        =              1.0
 
 gas.column.CO2.logstate     =                F
 gas.column.CO2.scale        =              1.0
 gas.column.CO2.sigma        =              1.0
 
 gas.column.OCS.logstate     =                F
 gas.column.OCS.scale        =              1.0
 gas.column.OCS.sigma        =              1.0
 
 
 # Forward model parameters
 
 fw.delnu                    =           0.10000
 fw.lshapemodel              =               0
 fw.solar_spectrum           =               T
 fw.pressure_shift           =               T
 fw.apod_fcn                 =               F
 fw.phase_fcn                =               F
 fw.emission                 =               F
 fw.isotope_separation       =               F
 
 
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =               F
 rt.convergence              =             0.1
 rt.max_iteration            =              17
 rt.wshift                   =               T
 rt.wshift.type              =               3
 rt.wshift.apriori           =           0.000
 rt.wshift.sigma             =           0.100
 rt.slope                    =               T
 rt.slope.apriori            =           0.000
 rt.slope.sigma              =           0.100
 rt.curvature                =               F
 rt.solshift                 =               T
 rt.solshift.apriori         =           0.000
 rt.solshift.sigma           =           0.100
 rt.phase                    = T
 rt.phase.apriori            = 0.0
 rt.phase.sigma              = 0.1
 rt.dwshift                  =              T
 rt.temperature              =               F
 
 # Microwindows and their parameters
 
 band                        =   1 2 3
 band.1.nu_start             =        2057.700
 band.1.nu_stop              =        2058.000
 band.1.zshift               =               F
 band.1.beam                 =               0
 band.1.calc_point_space     =       0.200E-03
 band.1.wave_factor          =           1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.3923
 band.1.apodization_code     =               0
 band.1.gasb                 = CO     O3    CO2    OCS
 
 band.2.nu_start             =        2069.560
 band.2.nu_stop              =        2069.760
 band.2.zshift               =               F
 band.2.beam                 =               0
 band.2.calc_point_space     =       0.200E-03
 band.2.wave_factor          =           1.000
 band.2.max_opd= 180.0000
 band.2.omega= 2.3923
 band.2.apodization_code     =               0
 band.2.gasb                 = CO     O3    CO2    OCS
 
 band.3.nu_start             =        2157.500
 band.3.nu_stop              =        2159.150
 band.3.zshift               =               F
 band.3.beam                 =               0
 band.3.calc_point_space     =       0.200E-03
 band.3.wave_factor          =           1.000
 band.3.max_opd= 180.0000
 band.3.omega= 2.3923
 band.3.apodization_code     =               0
 band.3.gasb                 = CO     O3    CO2  N2O   H2O
 
 
 out.level= 1
 out.gas_spectra= F
 out.k_matrix= T
 out.ak_matrix= T
 out.g_matrix= T
 out.sa_matrix= T
 out.shat_matrix= T
 out.seinv_vector= T
 out.retprofiles= T
 out.pbpfile= T
 out.statevec= T
 file.out.pbpfile= pbpfile
 file.out.statevec= statevec
 file.out.ak_matrix= ak.out
 file.out.k_matrix= k.out
 file.out.g_matrix= g.out
 file.out.kb_matrix= kb.out
 file.out.sa_matrix= sa.out
 file.out.shat_matrix= shat.out
 file.out.seinv_vector= seinv.out
 file.out.retprofiles= rprfs.table
 file.out.aprprofiles= aprfs.table
 file.out.summary= summary
 kb= T
 kb.slope= T
 kb.curvature= T
 kb.solshft= T
 kb.solstrnth= T
 kb.temperature= T
 kb.phase= T
 kb.omega= T
 kb.max_opd= T
 kb.zshift= T
 kb.sza= T
 out.aprprofiles = T
