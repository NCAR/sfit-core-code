# Filenames Section

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = t15asc.4
 file.in.modulation_fcn              = ils19990502b.dat
 file.in.phase_fcn                   = ils19990502b.dat
 file.in.sa_matrix                   = sa.input
 file.in.isotope                     = isotope.input
 file.in.solarlines                  = solar.dat
 file.in.linelist                    = 00775.968750-00785.431250.hbin

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
 gas.profile.list                    = CLONO2 O3 H2O

 gas.profile.CLONO2.correlation          = T
 gas.profile.CLONO2.correlation.type     = 2 	
 gas.profile.CLONO2.correlation.width    = 10.0
 gas.profile.CLONO2.correlation.minalt   = 0.0
 gas.profile.CLONO2.correlation.maxalt   = 120.0
 gas.profile.CLONO2.logstate             = F
 gas.profile.CLONO2.scale                = 1.0
 gas.profile.CLONO2.sigma                =
    0.133631    0.147442    0.163956    0.174608    0.185952    0.198107
    0.210912    0.225189    0.241121    0.258544    0.277778    0.298807
    0.320750    0.344214    0.366618    0.385758    0.401610    0.410997
    0.416667    0.421076    0.427179    0.431934    0.438529    0.445435
    0.450835    0.458349    0.464238    0.472456    0.481125    0.487950
    0.497519    0.505076    0.515491    0.526316    0.534798    0.546914
    0.556449    0.570173    0.584925    0.596635    0.613601    0.627160
    0.646955    0.668810    0.686414    0.712688    0.742229

#   1.67741   8.25106   8.79695   2.57177   2.04584   1.90738
#   1.96225   2.13191   1.86314   1.44480   1.14718   0.93399
#   0.75008   0.57258   0.41802   0.30693   0.24546   0.22902
#   0.23930   0.25677   0.27082   0.28831   0.32295   0.38278
#   0.47130   0.58013   0.69487   0.79841   0.85670   0.85100
#   0.82698   0.86343   0.99708   1.20216   1.39458   1.47518
#   1.37532   1.09181   0.74833   0.65863   0.69033   0.70980
#   0.71871   0.72384   0.73640   0.75367   0.77007


 gas.profile.O3.correlation          = T
 gas.profile.O3.correlation.type     = 2 	
 gas.profile.O3.correlation.width    = 5.0
 gas.profile.O3.correlation.minalt   = 0.0
 gas.profile.O3.correlation.maxalt   = 120.0
 gas.profile.O3.logstate             = F
 gas.profile.O3.scale                = 1.0
 gas.profile.O3.sigma                =
   1.10977   0.85691   0.72938   0.60894   0.78347   0.63890
   0.33472   0.14597   0.13249   0.13566   0.12104   0.10213
   0.09148   0.09670   0.11460   0.13507   0.15191   0.16455
   0.17319   0.17969   0.18350   0.18182   0.18002   0.18791
   0.20933   0.24030   0.27507   0.30406   0.31312   0.30290
   0.29240   0.29878   0.32205   0.33246   0.28702   0.20408
   0.13832   0.10347   0.08387   0.07151   0.06941   0.08321
   0.11214   0.14897   0.18890   0.23492   0.28632
 
 gas.profile.H2O.correlation          = T
 gas.profile.H2O.correlation.type     = 2 	
 gas.profile.H2O.correlation.width    = 5.0
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


 gas.column.list                     = CO2 
 gas.column.CO2.scale                = 1.0
 gas.column.CO2.sigma                = 1.0
 gas.column.H2O.scale                = 1.0
 gas.column.H2O.sigma                = 0.1
 gas.column.C2H2.scale               = 1.0
 gas.column.C2H2.sigma               = 1.0

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

 band                                = 1 2
 band.1.nu_start                     = 780.100
 band.1.nu_stop                      = 780.350
 band.1.zshift                       = F	
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.200
 band.1.calc_point_space             = 0.0005	
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 80.00
 band.1.omega                        = 4.067
 
 band.1.gasb                         = CLONO2 H2O
 band.1.tempretb                     = F

 band.2.nu_start                     = 780.700
 band.2.nu_stop                      = 781.300
 band.2.zshift                       = F	
 band.2.zshift.type                  = 1
 band.2.zshift.apriori               = 0.000
 band.2.zshift.sigma                 = 0.200
 band.2.calc_point_space             = 0.0005	
 band.2.wave_factor                  = 1.000
 band.2.max_opd                      = 80.00
 band.2.omega                        = 4.067
 
 band.2.gasb                         = CO2 O3 H2O CLONO2
 band.2.tempretb                     = F


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
 out.g_matrix			                   = T
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
 out.seinv_vector                    = T
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = T
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F
