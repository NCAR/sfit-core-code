# Filenames Section

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = t15asc.4
 file.in.modulation_fcn              = ils19990502b.dat
 
 file.in.phase_fcn                   = ils19990502b.dat
 
 file.in.sa_matrix                   = sa.input
 file.in.isotope                     = isotope.input
 file.in.solarlines                  = solar.dat
 file.in.linelist                    = 04034.769054-04043.110946.hbin

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
 file.out.g_matrix		     = d.complete

# Definition for retrieval gases

 gas.layers                          = 47
 gas.profile.list                    = HF H2O

 gas.profile.HF.correlation          = T
 gas.profile.HF.correlation.type     = 2 	
 gas.profile.HF.correlation.width    = 4.0
 gas.profile.HF.correlation.minalt   = 0.0
 gas.profile.HF.correlation.maxalt   = 120.0
 gas.profile.HF.logstate             = F
 gas.profile.HF.scale                = 1.0
 gas.profile.HF.sigma                =
    0.026726    0.029488    0.032791    0.034922    0.037190    0.039621
    0.042182    0.045038    0.048224    0.051709    0.055556    0.059761
    0.064150    0.068843    0.073324    0.077152    0.080322    0.082199
    0.083333    0.084215    0.085436    0.086387    0.087706    0.089087
    0.090167    0.091670    0.092848    0.094491    0.096225    0.097590
    0.099504    0.101015    0.103098    0.105263    0.106960    0.109383
    0.111290    0.114035    0.116985    0.119327    0.122720    0.125432
    0.129391    0.133762    0.137283    0.142538    0.148446

#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000    0.100000
#    0.100000    0.100000    0.100000    0.100000    0.100000

 gas.profile.H2O.correlation          = T
 gas.profile.H2O.correlation.type     = 2 	
 gas.profile.H2O.correlation.width    = 4.0
 gas.profile.H2O.correlation.minalt   = 0.0
 gas.profile.H2O.correlation.maxalt   = 120.0
 gas.profile.H2O.logstate             = F
 gas.profile.H2O.scale                = 1.0
 gas.profile.H2O.sigma                =
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000    0.500000
    0.500000    0.500000    0.500000    0.500000    0.500000


 gas.column.list                     = CH4 HDO
 gas.column.H2O.scale                = 1.0
 gas.column.H2O.sigma                = 0.1
 gas.column.CH4.scale                = 1.0
 gas.column.CH4.sigma                = 1.0
 gas.column.HDO.scale                = 1.0
 gas.column.HDO.sigma                = 0.1


# Forward model parameters

 fw.raytonly                         = F
 fw.isotope_separation               = T
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

 band                                = 1
 band.1.nu_start                     = 4038.81
 band.1.nu_stop                      = 4039.07
 band.1.zshift                       = F	
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.200
 band.1.calc_point_space             = 0.0005	
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 257.00
 band.1.omega                        = 2.3923
 band.1.apodization_code             = 0
 band.1.gasb                         = HF H2O CH4 HDO
 band.1.tempretb                     = F

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
 out.channel                         = F
 out.parm_vectors                    = F
 out.seinv_vector                    = F
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = T
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F
