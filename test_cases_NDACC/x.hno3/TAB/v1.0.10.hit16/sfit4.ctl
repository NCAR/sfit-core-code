# Filenames Section

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = t15asc.4
 file.in.modulation_fcn              = ils_poly2.mf
 file.in.phase_fcn                   = ils_poly2.mf
 file.in.sa_matrix                   = sa.input
 file.in.isotope                     = isotope.input
 file.in.solarlines                  = solar.dat
 file.in.linelist                    = 00863.009054-00874.040946.hbin

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

 gas.layers                          = 47
 gas.profile.list                    = HNO3 H2O

 gas.profile.HNO3.correlation          = T
 gas.profile.HNO3.correlation.type     = 2 	
 gas.profile.HNO3.correlation.width    = 2.0
 gas.profile.HNO3.correlation.minalt   = 0.0
 gas.profile.HNO3.correlation.maxalt   = 120.0
 gas.profile.HNO3.logstate             = F
 gas.profile.HNO3.scale                = 1.0
 gas.profile.HNO3.sigma                =
    0.133631    0.147442    0.163956    0.174608    0.185952    0.198107
    0.210912    0.225189    0.241121    0.258544    0.277778    0.298807
    0.320750    0.344214    0.366618    0.385758    0.401610    0.410997
    0.416667    0.421076    0.427179    0.431934    0.438529    0.445435
    0.450835    0.458349    0.464238    0.472456    0.481125    0.487950
    0.497519    0.505076    0.515437    0.526462    0.534828    0.546848
    0.556243    0.570173    0.585206    0.596338    0.613601    0.627456
    0.647117    0.668750    0.686156    0.712832    0.742065

#0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3
#0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3
#0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3
#0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 

 gas.profile.H2O.correlation          = T
 gas.profile.H2O.correlation.type     = 2 	
 gas.profile.H2O.correlation.width    = 4.0
 gas.profile.H2O.correlation.minalt   = 0.0
 gas.profile.H2O.correlation.maxalt   = 120.0
 gas.profile.H2O.logstate             = F
 gas.profile.H2O.scale                = 1.0
 gas.profile.H2O.sigma                =
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
    0.100000    0.100000    0.100000    0.100000    0.100000

 gas.column.list                     = OCS NH3 CO2
 gas.column.H2O.scale                = 1.0
 gas.column.H2O.sigma                = 0.2
 gas.column.OCS.scale                = 1.0
 gas.column.OCS.sigma                = 0.2
 gas.column.NH3.scale                = 1.0
 gas.column.NH3.sigma                = 0.2
 gas.column.CO2.scale                = 1.0
 gas.column.CO2.sigma                = 0.2

# Forward model parameters

 fw.raytonly                         = F
 fw.isotope_separation               = F
 fw.delnu                            = 0.1000 	
 fw.lshapemodel                      = 0
 fw.linemixing                       = T 	# new 
 fw.linemixing.gas                   = CO2
 fw.solar_spectrum                   = T
 fw.pressure_shift                   = T
 fw.apod_fcn                         = F
 fw.apod_fcn.type                    = 4
 fw.apod_fcn.order                   = 2
 fw.phase_fcn                        = F
 fw.phase_fcn.type                   = 4
 fw.phase_fcn.order                  = 2
 fw.emission                         = F
 fw.tips                             = F

# Retrieval parameter

rt                      = T 
rt.lm                   = F 
rt.lm.gamma_start       = 1.0e5 
rt.lm.gamma_inc         = 10.0 
rt.lm.gamma_dec         = 10.0 
rt.convergence          = 0.01 
rt.tolerance            = 0.1 
rt.max_iteration        = 25 
rt.dwshift              = F 
rt.wshift               = T 
rt.wshift.type          = 2 
rt.wshift.apriori       = 0.000 
rt.wshift.sigma         = 0.100 
rt.slope                = T 
rt.slope.apriori        = 0.000 
rt.slope.sigma          = 0.100 
rt.curvature            = F 
rt.curvature.apriori    = 0.000 
rt.curvature.sigma      = 0.100 
rt.phase                = T 
rt.phase.apriori        = 0.000 
rt.phase.sigma          = 0.010 
rt.apod_fcn             = F 
rt.apod_fcn.apriori     = 1.000 
rt.apod_fcn.sigma       = 0.200 
rt.phase_fcn            = F 
rt.phase_fcn.apriori    = 1.000 
rt.phase_fcn.sigma      = 0.200 
rt.solshift             = T 
rt.solshift.apriori     = 0.000 
rt.solshift.sigma       = 0.010 
rt.solstrnth            = T 
rt.solstrnth.apriori    = 0.000 
rt.solstrnth.sigma      = 0.010 
rt.temperature          = F 
rt.ifcalc_se            = T


# Kb derivative calculations

 kb                                  = T 
 kb.temperature                      = T
 kb.slope                            = F
 kb.curvature                        = F
 kb.solshft                          = F
 kb.solstrnth                        = F
 kb.phase                            = F
 kb.wshift                           = F
 kb.apod_fcn                         = T
 kb.phase_fcn                        = T
 kb.zshift                           = F
 kb.sza                              = T
 kb.omega                            = F
 kb.line                             = T
 kb.line.type                        = 1
 kb.line.gas                         = target


# Microwindows and their parameters

 band                                = 1
 band.1.nu_start                     = 867.050
 band.1.nu_stop                      = 870.000
 band.1.zshift                       = F	
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.200
 band.1.calc_point_space             = 0.0005	
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 257.00
 band.1.omega                        = 4.067
 band.1.apodization_code             = 0
 band.1.gasb                         = HNO3 H2O OCS NH3 CO2 
 band.1.tempretb                     = F
 band.1.beam                         = 1 2
 band.1.beam.model                   = PS
 band.1.beam.1.apriori               = -0.01 2.00 868.00 0.0
 band.1.beam.1.sigma                 =  1.0  0.0    1.0  0.0
 band.1.beam.2.apriori               =  0.01 1.0  868.25 0.0
 band.1.beam.2.sigma                 =  1.0  0.0    1.0  0.0



# Spectrum parameters

 #sp.snr                              = 1 
 #sp.snr.1.nu_start                   = 4038.86
 #sp.snr.1.nu_stop                    = 4039.05
 #sp.snr.1.snr                        = 150.0



# Output Files Section

 out.level                           = 1
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.g_matrix			             = T
 out.k_matrix                        = T
 out.shat_matrix                     = F
 out.retprofiles                     = T
 out.aprprofiles                     = T
 out.ab_matrix                       = F
 out.ak_matrix                       = T
 out.summary                         = T
 out.pbpfile                         = T
 out.channel                         = T
 out.parm_vectors                    = F
 out.seinv_vector                    = T
 out.sainv_matrix                    = F
 out.smeas_matrix                    = T
 out.ssmooth_matrix                  = T
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F
