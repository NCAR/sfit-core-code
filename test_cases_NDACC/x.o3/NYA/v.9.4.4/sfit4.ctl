# SFIT4 V0.9.4.4 input control file

 # Input filenames

 file.in.stalayers           = ./station.layers
 file.in.refprofile          = ./reference.prf
 file.in.spectrum            = ./spectrum
 file.in.isotope             = ./isotope.input
 file.in.solarlines          = ./solar.dat
 file.in.linelist            = ./00778.501722-01008.558278.hbin

 # Output filenames

 file.out.ak_matrix          = ak.out
 file.out.k_matrix           = k.out
 file.out.g_matrix           = g.out
 file.out.kb_matrix          = kb.out
 file.out.sa_matrix          = sa.out
 file.out.retprofiles        = rprfs.table
 file.out.aprprofiles        = aprfs.table
 file.out.summary            = summary
 file.out.seinv_vector       = seinv.vector
 file.out.pbpfile            = pbpfile
 file.out.retprofiles        = rprfs.table
 file.out.aprprofiles        = aprfs.table
 file.out.statevec           = statevec
 file.out.shat_matrix        = shat.out


 # Definition for Retrieval Gases

 gas.layers                  = 41
 gas.profile.list            = O3 H2O O3668 O3686
 gas.column.list             = CO2 C2H4

 gas.profile.O3.correlation         = T
 gas.profile.O3.correlation.type    = 2
 gas.profile.O3.correlation.width   = 6.0
 gas.profile.O3.correlation.minalt  = 0.000
 gas.profile.O3.correlation.maxalt  = 120.000

 gas.profile.O3.logstate     = F
 gas.profile.O3.scale        = 1.0
 gas.profile.O3.sigma        =
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.20000   0.20000   0.20000   0.20000   0.20000
 0.2

 gas.profile.H2O.correlation = F
 gas.profile.H2O.logstate    = F
 gas.profile.H2O.scale       = 1.0
 gas.profile.H2O.sigma       =
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1   0.1   0.1   0.1   0.1
 0.1

 gas.profile.O3668.correlation         = T
 gas.profile.O3668.correlation.type    = 2
 gas.profile.O3668.correlation.width   = 4.0
 gas.profile.O3668.correlation.minalt  = 0.000
 gas.profile.O3668.correlation.maxalt  = 120.000
 gas.profile.O3668.scale     = 1.0
 gas.profile.O3668.sigma     =
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2

 gas.profile.O3686.correlation         = T
 gas.profile.O3686.correlation.type    = 2
 gas.profile.O3686.correlation.width   = 4.0
 gas.profile.O3686.correlation.minalt  = 0.000
 gas.profile.O3686.correlation.maxalt  = 120.000
 gas.profile.O3686.logstate  = F
 gas.profile.O3686.scale     = 1.0
 gas.profile.O3686.sigma     =
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2   0.2   0.2   0.2   0.2
 0.2

 gas.column.CO2.scale        = 1.0
 gas.column.CO2.sigma        = 0.01

 gas.column.C2H4.scale       = 1.0
 gas.column.C2H4.sigma       = 0.05


 # Forward model parameters

 fw.delnu                    = 0.10000
 fw.lshapemodel              = 1
 fw.solar_spectrum           = T
 fw.pressure_shift           = T
 fw.apod_fcn                 = F
 fw.phase_fcn                = F
 fw.isotope_separation       = T


 # Retrieval Parameters

 rt                          = T
 rt.lm                       = T
 rt.lm.gamma_start           = 1.0e5
 rt.lm.gamma_inc             = 10.0
 rt.lm.gamma_dec             = 10.0
 rt.convergence              = 0.1
 rt.max_iteration            = 17
 rt.wshift                   = T
 rt.wshift.type              = 2
 rt.wshift.apriori           = 0.000
 rt.wshift.sigma             = 0.100
 rt.slope                    = T
 rt.slope.apriori            = 0.000
 rt.slope.sigma              = 0.01
 rt.curvature                = T
 rt.curvature.apriori        = 0.000
 rt.curvature.sigma          = 0.01
 rt.phase                    = F
 rt.phase.apriori            = 0.000
 rt.phase.sigma              = 0.0100
 rt.dwshift                  = T
 rt.solshift                 = T
 rt.solshift.apriori         = 0.00
 rt.solshift.sigma           = 0.01
 rt.temperature              = F
 rt.temperature.sigma        =
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001
 0.001

 # Microwindows and their Specific Parameters

 band                        = 1 2 3 4

 band.1.nu_start             = 782.56
 band.1.nu_stop              = 782.86
 band.1.zshift		           = T
 band.1.zshift.type          = 2
 band.1.beam                 =
 band.1.beam.model           = PS
 band.1.beam.1.apriori       = 0.01 0.02 782.56 1.0
 band.1.beam.1.sigma         = 1.0  0.0 1.0 0.0
 band.1.calc_point_space     = 0.500E-03
 band.1.wave_factor          = 1.00D0
 band.1.max_opd              = 180.0000
 band.1.omega                = 2.3923
 band.1.apodization_code     = 0
 band.1.gasb                 = O3 H2O CO2 O3668 O3686
 band.1.tempretb             = F

 band.2.nu_start             = 788.85
 band.2.nu_stop              = 789.37
 band.2.zshift		           = T
 band.2.zshift.type          = 2
 band.2.beam                 =
 band.2.calc_point_space     = 0.500E-03
 band.2.wave_factor          = 1.000
 band.2.max_opd              = 180.0000
 band.2.omega                = 2.3923
 band.2.apodization_code     = 0
 band.2.gasb                 = O3 H2O CO2 O3668 O3686

 band.3.nu_start             = 993.3
 band.3.nu_stop              = 993.8
 band.3.zshift		           = T
 band.3.zshift.type          = 2
 band.3.beam                 = 1
 band.3.beam.model           = PS
 band.3.beam.1.apriori       = 0.05 0.2 993.52 1.0
 band.3.beam.1.sigma         = 1.0  0.0 1.0    0.0
 band.3.calc_point_space     = 0.500E-03
 band.3.wave_factor          = 1.000
 band.3.max_opd              = 180.0000
 band.3.omega                = 2.3923
 band.3.apodization_code     = 0
 band.3.gasb                 = O3 H2O CO2 C2H4 O3668 O3686

 band.4.nu_start             = 1000.000
 band.4.nu_stop              = 1004.500
 band.4.zshift		           = T
 band.4.zshift.type          = 1
 band.4.zshift.apriori       = 0.000
 band.4.zshift.sigma         = 0.0200
 band.4.beam                 =
 band.4.calc_point_space     = 0.500E-03
 band.4.wave_factor          = 1.000
 band.4.max_opd              = 180.0000
 band.4.omega                = 2.3923
 band.4.apodization_code     = 0
 band.4.gasb                 = O3 H2O CO2 O3668 O3686


# Overwrite Input SNR for Specific Regions

 sp.snr                      =
 sp.snr.1.nu_start           = 1000.9
 sp.snr.1.nu_stop            = 1001.25
 sp.snr.1.snr                = 1.0
 sp.snr.2.nu_start           = 1003.4
 sp.snr.2.nu_stop            = 1004.6
 sp.snr.2.snr                = 1.0


# Output Files

 out.level                   = 1
 out.smeas_matrix            = T
 out.gas_spectra             = T
 out.k_matrix                = T
 out.ak_matrix               = T
 out.g_matrix                = T
 out.sa_matrix               = T
 out.retprofiles             = T
 out.aprprofiles             = T
 out.seinv_vector            = T
 out.raytrace                = T
 out.raytrace.type           = 1


# Derivatives for Uncertainty Calculations

 kb                          = F
 kb.slope                    = T
 kb.curvature                = T
 kb.solshft                  = T
 kb.solstrnth                = T
 kb.temperature              = T
 kb.phase                    = T
 kb.omega                    = T
 kb.zshift                   = T
 kb.sza                      = T
 kb.line                     = T
 kb.line.type                = 1
 kb.line.gas                 = retrieval
