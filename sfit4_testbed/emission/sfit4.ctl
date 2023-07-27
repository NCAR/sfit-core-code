file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = spectrum 
file.in.solarlines = solar.dat 
file.in.linelist = 00935.856908-01114.143092.hbin 
gas.layers = 39 
gas.isoflag = T 
gas.profile.list = O3 H2O 
gas.column.list = CO2 CH4 N2O 
gas.profile.O3.correlation = T 
gas.profile.O3.correlation.type = 6 
gas.profile.O3.correlation.lambda = 100 
gas.profile.O3.logstate = F 
gas.profile.O3.scale = 1.0 
gas.profile.O3.sigma = 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
gas.profile.H2O.logstate = F 
gas.profile.H2O.correlation = T 
gas.profile.H2O.correlation.type = 6 
gas.profile.H2O.correlation.lambda = 100 
gas.profile.H2O.scale = 1.0 
gas.profile.H2O.sigma = 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 
gas.column.CO2.scale = 1.0 
gas.column.CO2.sigma = 1.0 
gas.column.N2O.scale = 1.0 
gas.column.N2O.sigma = 1.0 
gas.column.CH4.scale = 1.0 
gas.column.CH4.sigma = 1.0 
fw.isotope_separation = F 
fw.delnu = 1.0 
fw.lshapemodel = 1 
fw.pressure_shift = T 
fw.apod_fcn = F 
fw.emission = T 
fw.emission.T_infinity = 2.7 
fw.emission.object = .e. 
fw.emission.normalized = F 
fw.mtckd_continuum = T 
rt = T
rt.lm = T 
rt.lm.gamma_start = 1.0e5 
rt.lm.gamma_inc = 1.0e1 
rt.lm.gamma_dec = 1.0e1 
rt.convergence = 0.5 
rt.max_iteration = 20 
rt.dwshift = F 
rt.wshift = T 
rt.wshift.type = 2 
rt.wshift.apriori = 0.000 
rt.wshift.sigma = 0.100 
rt.slope = F 
rt.slope.apriori = 1.000 
rt.slope.sigma = 0.100 
rt.curvature = F 
rt.curvature.apriori = 0.000 
rt.curvature.sigma = 0.100 
rt.solshift = F 
rt.solshift.apriori = 0.000 
rt.solshift.sigma = 0.100 
rt.phase = F 
rt.phase.apriori = 1.000 
rt.phase.sigma = 0.500 
rt.apod_fcn = F 
rt.apod_fcn.type = 0 
rt.apod_fcn.apriori = 0.200 
rt.apod_fcn.sigma = 0.200 
rt.phase_fcn = F 
rt.phase_fcn.type = 0 
rt.phase_fcn.apriori = 1.000 
rt.phase_fcn.sigma = 0.200 
rt.temperature = F 
rt.temperature.sigma = 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 0.001 
kb = T 
kb.slope = F 
kb.curvature = F 
kb.solshft = F 
kb.solstrnth = F 
kb.temperature = T 
kb.phase = F 
kb.omega = F 
kb.zshift = F 
kb.sza = F 
kb.line = T 
kb.line.type = 1 
kb.line.gas = O3 
band = 1 2
band.1.nu_start = 950.0 
band.1.nu_stop = 1050.0 
band.1.zshift = F 
band.1.beam = 0 
band.1.calc_point_space = 0.01 
band.1.wave_factor = 1.000 
band.1.max_opd = 1.035 
band.1.omega = 46.0 
band.1.apodization_code = 0 
band.1.gasb = O3 H2O CO2 CH4 N2O 
band.1.tempretb = F 
band.2.nu_start = 1051.0 
band.2.nu_stop = 1100.0 
band.2.zshift = F 
band.2.beam = 0 
band.2.calc_point_space = 0.01 
band.2.wave_factor = 1.000 
band.2.max_opd = 1.035 
band.2.omega = 46.0 
band.2.apodization_code = 0 
band.2.gasb = O3 H2O CO2 CH4 N2O 
band.2.tempretb = F 
sp.snr.1.nu_start = 950.0 
sp.snr.1.nu_stop = 1100.0 
sp.snr.1.snr = 0.0002 
out.level = 1 
out.gas_spectra = T 
out.gas_spectra.type = 1 
out.pbpfile = T 
out.smeas_matrix = F 
out.k_matrix = T 
out.ak_matrix = T 
out.g_matrix = T 
out.sa_matrix = T 
out.shat_matrix = F 
out.statevec = T 
out.retprofiles = F 
out.aprprofiles = F 
out.parm_vectors = F 
out.ab_matrix = F 
out.summary = T 
out.seinv_vector = T 
file.out.ak_matrix = ak.out 
file.out.pbpfile = pbpfile.out 
file.out.k_matrix = k.out 
file.out.g_matrix = g.out 
file.out.kb_matrix = kb.out 
file.out.sa_matrix = sa.out 
file.out.statevec = statevec 
file.out.retprofiles = rprfs.table 
file.out.aprprofiles = aprfs.table 
file.out.summary = summary 
file.out.seinv_vector = seinv.vector 
fw.tips = F 
