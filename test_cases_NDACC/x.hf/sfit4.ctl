file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = spectrum 
file.in.isotope = isotope.input 
file.in.solarlines = ../../linelist/solar/120621/solar.dat 
file.in.linelist = 03996.841422-04114.258578.hbin
file.out.pbpfile = pbpfile 
file.out.statevec = statevec 
file.out.k_matrix = k.output 
file.out.kb_matrix = kb.out 
file.out.sa_matrix = sa.complete 
file.out.retprofiles = rprfs.table 
file.out.aprprofiles = aprfs.table 
file.out.ab_matrix = ab.output 
file.out.summary = summary 
file.out.parm_vectors = parm.output 
file.out.seinv_vector = seinv.output 
file.out.sainv_matrix = sainv.complete 
file.out.smeas_matrix = smeas.target 
file.out.shat_matrix = shat.complete 
file.out.ssmooth_matrix = ssmooth.target 
file.out.ak_matrix = ak.target 
file.out.solarspectrum = solar.output 
gas.layers = 41 
gas.profile.list = HF 
gas.profile.HF.correlation = F 
gas.profile.HF.logstate = F 
gas.profile.HF.scale = 1.0d0 
gas.profile.HF.sigma = 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000 
gas.column.list = H2O CH4 HDO 
gas.column.H2O.scale = 0.2d0 
gas.column.H2O.sigma = 0.1d0 
gas.column.HDO.scale = 0.2d0 
gas.column.HDO.sigma = 0.01d0 
gas.column.CH4.scale = 1.0d0 
gas.column.CH4.sigma = 0.1d0 
fw.isotope_separation = T 
fw.delnu = 0.100 
fw.lshapemodel = 0
fw.linemixing = F
fw.linemixing.gas = CO2 
fw.solar_spectrum = T 
fw.pressure_shift = T 
fw.apod_fcn = F 
fw.phase_fcn = F 
fw.emission = F 
rt = F
rt.lm = F 
rt.convergence = 0.1 
rt.max_iteration = 30 
rt.dwshift = F 
rt.wshift = F
rt.wshift.type = 3 
rt.wshift.apriori = 0.000 
rt.wshift.sigma = 0.100 
rt.slope = T 
rt.slope.apriori = 0.000 
rt.slope.sigma = 0.100 
rt.curvature = T 
rt.curvature.apriori = 0.000 
rt.curvature.sigma = 0.100 
rt.phase = F 
rt.phase.apriori = 0.000 
rt.phase.sigma = 0.100 
rt.apod_fcn = F 
rt.apod_fcn.apriori = 0.000 
rt.apod_fcn.sigma = 0.000 
rt.phase_fcn = F 
rt.phase_fcn.apriori = 0.000 
rt.phase_fcn.sigma = 0.000 
rt.solshift = F 
rt.solshift.apriori = 0.000 
rt.solshift.sigma = 0.500 
rt.solstrnth = F 
rt.solstrnth.apriori = 0.000 
rt.solstrnth.sigma = 0.500 
rt.temperature = F 
kb = F
kb.temperature = T 
kb.slope = T 
kb.curvature = T 
kb.solshft = T 
kb.solstrnth = T 
kb.phase = T 
kb.dwshift = F 
kb.wshift = F 
kb.apod_fcn = F 
kb.phase_fcn = F 
kb.zshift = T 
kb.sza = T 
kb.omega = T 
kb.max_opd = T 
kb.line = T 
kb.line.type = 1 
kb.line.gas = HF 
band = 3 1 2
band.1.nu_start = 4038.7 
band.1.nu_stop = 4039.2 
band.1.zshift = F 
band.1.calc_point_space = 0.900E-03 
band.1.wave_factor = 1.000 
band.1.max_opd = 180.0 
band.1.omega = 2.3923 
band.1.apodization_code = 0 
band.1.beam = 0 
band.1.gasb = HF H2O CH4 HDO 
band.1.tempretb = F 
band.2.nu_start = 4109.4 
band.2.nu_stop = 4110.2 
band.2.zshift = F 
band.2.calc_point_space = 0.900E-03 
band.2.wave_factor = 1.000 
band.2.max_opd = 180.0 
band.2.omega = 2.3923 
band.2.apodization_code = 0 
band.2.beam = 0 
band.2.gasb = HF H2O HDO CH4 
band.2.tempretb = F 
band.3.nu_start = 4000.9 
band.3.nu_stop = 4001.45 
band.3.zshift = F 
band.3.calc_point_space = 0.900E-03 
band.3.wave_factor = 1.000 
band.3.max_opd = 180.0 
band.3.omega = 2.3923 
band.3.apodization_code = 0 
band.3.beam = 0 
band.3.gasb = HF H2O HDO CH4 
band.3.tempretb = F 
out.level = 1 
out.gas_spectra = T 
out.gas_spectra.type = 1 
out.sa_matrix = T 
out.statevec = T 
out.k_matrix = T 
out.shat_matrix = T 
out.retprofiles = T 
out.aprprofiles = T 
out.ab_matrix = F 
out.ak_matrix = F 
out.summary = T 
out.pbpfile = T 
out.channel = F 
out.parm_vectors = T 
out.seinv_vector = F 
out.sainv_matrix = T 
out.smeas_matrix = F 
out.ssmooth_matrix = F 
out.raytrace = T 
out.raytrace.type = 1 
out.solarspectrum = F 
out.levmardet = F 
out.xscdetail = F 
