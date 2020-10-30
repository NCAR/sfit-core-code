# Filenames Section

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = t15asc.4
 file.in.modulation_fcn              = ils_poly2.mf
 file.in.phase_fcn                   = ils_poly2.phs
 file.in.sa_matrix                   = sa.input
 file.in.isotope                     = isotope.input
 file.in.solarlines                  = solar.dat
 file.in.linelist                    = 03263.999054-03303.640946.hbin

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
 file.out.g_matrix   	             = d.complete

# Definition for retrieval gases

 gas.layers                              = 47
 gas.profile.list                        = HCN
 gas.profile.HCN.correlation             = T
 gas.profile.HCN.correlation.type        = 2
 gas.profile.HCN.correlation.width       = 6.0
 gas.profile.HCN.correlation.minalt      = 0.0
 gas.profile.HCN.correlation.maxalt      = 120.0
 gas.profile.HCN.logstate                = F
 gas.profile.HCN.scale                   = 1.0
 gas.profile.HCN.sigma                   =
  0.1336  0.1474  0.1640  0.1746  0.1860  0.1981  0.2109  0.2252  0.2411  0.2585
  0.2778  0.2988  0.3208  0.3442  0.3666  0.3858  0.4016  0.4110  0.4167  0.4211
  0.4272  0.4319  0.4385  0.4454  0.4508  0.4583  0.4642  0.4725  0.4811  0.4880
  0.4975  0.5051  0.5160  0.5276  0.5376  0.5508  0.5618  0.5774  0.5938  0.6081
  0.6275  0.6439  0.6670  0.6940  0.7157  0.7487  0.7866

 gas.column.list                     = H2O C2H2 CO2 O3
 gas.column.H2O.scale                = 1.0
 gas.column.H2O.sigma                = 1.0
 gas.column.C2H2.scale               = 1.0
 gas.column.C2H2.sigma               = 1.0
 gas.column.CO2.scale                = 1.0
 gas.column.CO2.sigma                = 1.0
 gas.column.O3.scale                 = 1.0
 gas.column.O3.sigma                 = 1.0


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
 kb.line.gas                         = retrieval


# Microwindows and their parameters

 band                                = 1 2 3

 band.1.nu_start                     = 3268.0400
 band.1.nu_stop                      = 3268.4000
 band.1.zshift                       = F
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.200
 band.1.calc_point_space             = 0.0005
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 257.00
 band.1.omega                        = 2.3923
 band.1.apodization_code             = 0
 band.1.gasb                         = HCN H2O C2H2 O3
 band.1.tempretb                     = F

 band.2.nu_start                     = 3287.1000
 band.2.nu_stop                      = 3287.3500
 band.2.zshift                       = F
 band.2.zshift.type                  = 1
 band.2.zshift.apriori               = 0.000
 band.2.zshift.sigma                 = 0.200
 band.2.calc_point_space             = 0.0005
 band.2.wave_factor                  = 1.000
 band.2.max_opd                      = 257.00
 band.2.omega                        = 2.3923
 band.2.apodization_code             = 0
 band.2.gasb                         = HCN H2O C2H2 CO2
 band.2.tempretb                     = F

 band.3.nu_start                     = 3299.4000
 band.3.nu_stop                      = 3299.6000
 band.3.zshift                       = F
 band.3.zshift.type                  = 1
 band.3.zshift.apriori               = 0.000
 band.3.zshift.sigma                 = 0.200
 band.3.calc_point_space             = 0.0005
 band.3.wave_factor                  = 1.000
 band.3.max_opd                      = 257.00
 band.3.omega                        = 2.3923
 band.3.apodization_code             = 0
 band.3.gasb                         = HCN H2O
 band.3.tempretb                     = F

 sp.snr                              = 
 sp.snr.1.nu_start                   = 3268.0000
 sp.snr.1.nu_stop                    = 3268.3800
 sp.snr.1.snr                        = 100.0 #	changed	from 200 to check is oscillations are reduced

 sp.snr.2.nu_start                   = 3287.0000
 sp.snr.2.nu_stop                    = 3287.4800
 sp.snr.2.snr                        = 100.0  # changed from 200 to check is oscillations are reduced

 sp.snr.3.nu_start                   = 3200.0000
 sp.snr.3.nu_stop                    = 3400.0000
 sp.snr.3.snr                        = 120.0 #	changed	from 200 to check is oscillations are reduced


# Output Files Section

 out.level                           = 1
 out.gas_spectra                     = T
 out.gas_spectra.type                = 1
 out.sa_matrix                       = T
 out.statevec                        = T
 out.g_matrix			     = T
 out.k_matrix                        = T
 out.shat_matrix                     = F
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
