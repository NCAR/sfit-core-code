# General

 file.in.stalayers                   = station.layers
 file.in.refprofile                  = reference.prf
 file.in.spectrum                    = t15asc.4
 file.in.modulation_fcn              = ils_poly.mf
 file.in.phase_fcn                   = ils_poly.phs
 file.in.sa_matrix                   = sa.input
 file.in.solarlines                  = /Users/jamesw/FDP/sfit/400/linelist.all/solar/120621/solar.dat
 file.in.isotope                     = isotope.input.co2
# file.in.linelist                    = ./04877.895100-04892.104900.hbin
 file.in.linelist                    = ./04867.895100-04892.104900.hbin

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

# Definition for retrieval gases

 gas.layers                          = 44
 gas.profile.list                    = CO2

 gas.profile.CO2.scale               = 1.0D0
 gas.profile.CO2.correlation         = T
 gas.profile.CO2.correlation.type    = 2
 gas.profile.CO2.correlation.width   = 500.0
 gas.profile.CO2.correlation.minalt  = 0.000
 gas.profile.CO2.correlation.maxalt  = 120.000
 gas.profile.CO2.logstate            = F
 gas.profile.CO2.sigma               =
    0.00534 0.00590 0.00656 0.00698 0.00744
    0.00792 0.00844 0.00901 0.00965 0.01034
    0.01111 0.01195 0.01283 0.01377 0.01466
    0.01543 0.01606 0.01644 0.01667 0.01684
    0.01709 0.01728 0.01754 0.01782 0.01803
    0.01833 0.01857 0.01890 0.01925 0.01952
    0.01990 0.02020 0.02062 0.02102 0.02136
    0.02182 0.02218 0.02270 0.02327 0.02372
    0.02434 0.02486 0.02559 0.02642

 gas.column.list                     = O13CO CO18O HDO H2O H218O

 gas.column.O13CO.scale                     = 1.0
 gas.column.O13CO.sigma                     = 0.1

 gas.column.CO18O.scale                     = 1.0
 gas.column.CO18O.sigma                     = 0.1

 gas.column.HDO.scale                       = 1.0
 gas.column.HDO.sigma                       = 0.1

 gas.column.H2O.scale                       = 1.0
 gas.column.H2O.sigma                       = 0.5

 gas.column.H218O.scale                     = 1.0
 gas.column.H218O.sigma                     = 0.1

 gas.column.H217O.scale                     = 1.0
 gas.column.H217O.sigma                     = 0.1



# Forward model parameters
 fw.raytonly                         = F
 fw.isotope_separation               = T
 fw.delnu                            = 0.1D0
 fw.lshapemodel                      = 0
 fw.linemixing                       = T
 fw.linemixing.gas                   = CO2
 fw.solar_spectrum                   = T
 fw.pressure_shift                   = T
 fw.apod_fcn                         = T
 fw.apod_fcn.type                    = 2
 fw.apod_fcn.order                   = 5
 fw.phase_fcn                        = T
 fw.phase_fcn.type                   = 2
 fw.phase_fcn.order                  = 5
 fw.emission                         = F
# fw.emission.T_infinity = 6000.0
# fw.emission.object = E
# fw.emission.normalized = T


# Retrieval parameter

 rt                                  = T
 rt.lm                               = T
 rt.lm.gamma_start                   = 10000.
 rt.lm.gamma_dec                     = 10.0
 rt.lm.gamma_inc                     = 10.0
 rt.convergence                      = 0.05
 rt.tolerance                        = 0.05
 rt.max_iteration                    = 21
 rt.wshift                           = T
 rt.wshift.type                      = 1
 rt.wshift.apriori                   = 0.000
 rt.wshift.sigma                     = 0.0100
 rt.dwshift                          = T
 rt.slope                            = T
 rt.slope.apriori                    = 0.000
 rt.slope.sigma                      = 0.100
 rt.curvature                        = F
 rt.curvature.apriori                = 0.000
 rt.curvature.sigma                  = 0.100
 rt.phase                            = F
 rt.phase.apriori                    = 0.000
 rt.phase.sigma                      = 0.000
 rt.apod_fcn                         = T
 rt.apod_fcn.apriori                 = 0.000
 rt.apod_fcn.sigma                   = 0.100
 rt.phase_fcn                        = T
 rt.phase_fcn.apriori                = 0.000
 rt.phase_fcn.sigma                  = 0.100
 rt.solshift                         = T
 rt.solshift.apriori                 = 0.00
 rt.solshift.sigma                   = 0.100
 rt.solstrnth                        = F
 rt.solstrnth.apriori                = 0.000
 rt.solstrnth.sa                     = 0.010
 rt.ifcalc_se                        = F
 rt.temperature                      = F
 rt.temperature.sigma                =
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    1.0 1.0 1.0 1.0 1.0 1.0 1.0



# Microwindows and their parameters

 band                                = 1
 band.1.nu_start                     = 4872.000
 band.1.nu_stop                      = 4888.000
 band.1.zshift                       = T
 band.1.zshift.type                  = 1
 band.1.zshift.apriori               = 0.000
 band.1.zshift.sigma                 = 0.0100
 band.1.calc_point_space             = 0.900E-03
 band.1.wave_factor                  = 1.000
 band.1.max_opd                      = 100.000
 band.1.omega                        = 1.914
 band.1.apodization_code             = 0
 band.1.snr                          = 800.0
 band.1.beam                         = 1
 band.1.beam.model                   = IP
 band.1.beam.1.apriori               = 0.01161 1.97530 4882.98981 0.00000
 band.1.beam.1.sa                    = 0.01 0.01 0.01 0.000
 band.1.gasb                         = CO2 O13CO CO18O HDO H2O H218O
 band.1.tempretb                     = F

# spectrum parameters

 sp.snr                              =
 sp.snr.1.nu_start                   = 4882.0
 sp.snr.1.nu_stop                    = 4888.0
 sp.snr.1.snr                        = 200.0
 sp.snr.2.nu_start                   =
 sp.snr.2.nu_stop                    =
 sp.snr.2.snr                        = 100.0


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
 out.seinv_vector                    = T
 out.sainv_matrix                    = F
 out.smeas_matrix                    = F
 out.ssmooth_matrix                  = F
 out.raytrace                        = T
 out.raytrace.type                   = 3
 out.solarspectrum                   = F
 out.levmardet                       = F
 out.xscdetail                       = F

# Kb derivative calculations

 kb                                  = F
 kb.temperature                      = F
 kb.slope                            = T
 kb.curvature                        = T
 kb.solshft                          = T
 kb.solstrnth                        = T
 kb.phase                            = T
 kb.wshift                           = T
 kb.apod_fcn                         = T
 kb.phase_fcn                        = T
 kb.zshift                           = T
 kb.sza                              = T
 kb.omega                            = T
 kb.line                             = T
 kb.line.type                        = 1
 kb.line.gas                         = CO2

