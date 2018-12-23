 
 # General
 
 file.in.stalayers              = station.layers
 file.in.refprofile             = reference.prf
 file.in.spectrum               = spectrum
 file.in.isotope                = isotope.inp
 file.in.solarlines= solar.dat
 file.in.linelist= 02910.549356-02928.965644.hbin
 
 # Definition for retrieval gases
 
 gas.layers                   =              48
 gas.profile.list             = NO2
 gas.column.list              = CH4 CH3D H2O HDO H218O #OCS C2H6
 gas.profile.NO2.correlation  = T
 gas.profile.NO2.correlation.type  = 2
 gas.profile.NO2.correlation.width  =           4.000
 gas.profile.NO2.correlation.minalt  =           0.000
 gas.profile.NO2.correlation.maxalt  =         100.000
 gas.profile.NO2.logstate     =               F
 gas.profile.NO2.scale        = 1.0
 gas.profile.NO2.sigma        =
 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2
 
 gas.column.CH4.logstate            =               F
 gas.column.CH4.scale               = 1.0
 gas.column.CH4.sigma               = 1.0
 
 gas.column.CH3D.logstate               =              F
 gas.column.CH3D.scale               = 1.0
 gas.column.CH3D.sigma               = 1.0
 
 gas.column.H2O.logstate            =               F
 gas.column.H2O.scale               = 1.0
 gas.column.H2O.sigma               = 1.0
 
 gas.column.H218O.logstate            =               F
 gas.column.H218O.scale               = 1.0
 gas.column.H218O.sigma               = 1.0
 
 gas.column.H217O.logstate            =               F
 gas.column.H217O.scale               = 1.0
 gas.column.H217O.sigma               = 1.0
 
 gas.column.HDO.logstate            =               F
 gas.column.HDO.scale               = 1.0
 gas.column.HDO.sigma               = 1.0
 
 gas.column.OCS.logstate            =               F
 gas.column.OCS.scale               = 1.0
 gas.column.OCS.sigma               = 1.0
 
 gas.column.C2H6.logstate            =               F
 gas.column.C2H6.scale               = 1.0
 gas.column.C2H6.sigma               = 1.0
 
 # Forward model parameters
 
 fw.delnu                    = 0.1
 fw.lshapemodel              = 0
 fw.solar_spectrum           = T
 fw.pressure_shift           = F
 fw.apod_fcn                 = F
 fw.phase_fcn                = F
 fw.emission                 =               F
 fw.isotope_separation       =               T
 
 # Retrieval parameter
 
 rt                          =               T
 rt.lm                       =              F
 rt.lm.gamma_start           =      100000.0
 rt.lm.gamma_dec             =          10.0
 rt.lm.gamma_inc             =          10.0
 rt.convergence              =           0.1
 rt.max_iteration            =              20
 rt.wshift                   =  T
 rt.wshift.type              =  3
 rt.wshift.apriori           =  0.000
 rt.wshift.sigma             =  0.100
 rt.slope                    =	T
 rt.slope.apriori            =  0.000
 rt.slope.sigma              =  0.100
 rt.curvature                    =	T
 rt.curvature.apriori            =  0.000
 rt.curvature.sigma              =  0.100
 rt.phase                    = T
 rt.phase.apriori            = 0.000
 rt.phase.sigma              = 0.100
 rt.solshift		     =		    T
 rt.solshift.apriori	     =             0.0
 rt.solshift.sigma           =             1.0
 rt.solstrnth		     =		    T
 rt.solstrnth.apriori	     =             0.0
 rt.solstrnth.sigma           =             1.0
 
 # Microwindows and their parameters
 
 band                        = 1 2 3 5
 band.1.nu_start             =         2914.590
 band.1.nu_stop              =         2914.707
 band.1.zshift               =         F
 band.1.calc_point_space    =       0.9e-3
 band.1.wave_factor= 0.999998990000
 band.1.max_opd= 180.0000
 band.1.omega= 2.7512
 band.1.apodization_code     =               0
 band.1.gasb                 = NO2 CH4 H2O #C2H6
 
 band.2.nu_start             =         2918.100
 band.2.nu_stop              =         2918.350
 band.2.zshift               =              F
 band.2.calc_point_space     =       0.9e-3
 band.2.wave_factor= 0.999998990000
 band.2.max_opd= 180.0000
 band.2.omega= 2.7512
 band.2.apodization_code     =               0
 band.2.gasb                 = NO2 H2O CH4 #C2H6
 
 band.3.nu_start             =         2919.400
 band.3.nu_stop              =         2919.650
 band.3.zshift               =              F
 band.3.calc_point_space     =       0.9e-3
 band.3.wave_factor= 0.999998990000
 band.3.max_opd= 180.0000
 band.3.omega= 2.7512
 band.3.apodization_code     =               0
 band.3.gasb                 = NO2 CH4 H2O HDO H218O #OCS C2H6
 
 
 band.4.nu_start             =         2922.360
 band.4.nu_stop              =         2922.750
 band.4.zshift               =              F
 band.4.calc_point_space     =       0.9e-3
 band.4.wave_factor         =           1.000005
 band.4.max_opd = 180.0000
 band.4.omega= 2.3923
 band.4.apodization_code    =               0
 band.4.gasb                 = NO2 CH4 HDO CH3D H218O
 
 
 band.5.nu_start             =         2924.750
 band.5.nu_stop              =         2924.925
 band.5.zshift               =         F
 band.5.calc_point_space     =       0.9e-3
 band.5.wave_factor= 0.999998990000
 band.5.max_opd= 180.0000
 band.5.omega= 2.7512
 band.5.apodization_code     =               0
 band.5.gasb                 = NO2 CH4 H2O CH3D HDO #OCS
 
 
 
 
 
 out.level= 1
 out.gas_spectra= F
 out.k_matrix= T
 out.ak_matrix= T
 out.g_matrix= T
 out.sa_matrix= T
 out.shat_matrix= T
 out.seinv_vector= T
 out.retprofiles= T
 out.aprprofiles= T
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
 kb.line = T
 kb.line.type = 1
 kb.line.gas = NO2
 kb.slope= T
 kb.curvature= T
 kb.solshft= T
 kb.solstrnth= T
 kb.temperature= T
 kb.phase = T
 kb.omega= T
 kb.zshift= T
 kb.sza= T