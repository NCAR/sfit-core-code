file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = spectrum 
file.in.isotope = isotope.input 
file.in.solarlines = solar.dat 
file.in.linelist = 00775.791722-00786.928278.hbin 
gas.layers = 48 
gas.profile.list = CLONO2 O3 
gas.column.list = HNO3 H2O CO2 
gas.profile.CLONO2.correlation = T 
gas.profile.CLONO2.correlation.type = 2 
gas.profile.CLONO2.correlation.width = 8.0 
gas.profile.CLONO2.correlation.minalt = 0.0 
gas.profile.CLONO2.correlation.maxalt = 120.0 
gas.profile.CLONO2.logstate = F 
gas.profile.CLONO2.scale = 1.0 
gas.profile.CLONO2.sigma = 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0 
gas.profile.O3.correlation = F 
gas.profile.O3.scale = 1.0 
gas.profile.O3.sigma = 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7   0.7   0.7 0.7   0.7   0.7 
gas.column.HNO3.scale = 1.0 
gas.column.HNO3.sigma = 1.0 
gas.column.H2O.scale = 1.0 
gas.column.H2O.sigma = 1.0 
gas.column.CO2.scale = 1.0 
gas.column.CO2.sigma = 1.0 
fw.delnu = 0.10000 
fw.lshapemodel = 0 
fw.solar_spectrum = F 
fw.pressure_shift = T 
fw.apod_fcn = F 
fw.phase_fcn = F 
fw.isotope_separation = F 
rt = T 
rt.lm = T 
rt.lm.gamma_start = 1e5 
rt.lm.gamma_inc = 10 
rt.lm.gamma_dec = 10 
rt.convergence = 0.1 
rt.max_iteration = 15 
rt.wshift = T 
rt.wshift.type = 3 
rt.wshift.apriori = 0.000 
rt.wshift.sigma = 0.100 
rt.slope = T 
rt.slope.apriori = 0.000 
rt.slope.sigma = 0.100 
rt.curvature = F 
rt.curvature.apriori = 0.000 
rt.curvature.sigma = 0.100 
rt.phase = T 
rt.phase.apriori = 0.000 
rt.phase.sigma = 0.200 
rt.dwshift = F 
band = 1 2 
band.1.nu_start = 779.85 
band.1.nu_stop = 780.45 
band.1.zshift = F 
band.1.beam = 0 
band.1.calc_point_space = 0.500E-03 
band.1.wave_factor = 0.999998990000 
band.1.max_opd = 180.0000 
band.1.omega = 2.3923 
band.1.apodization_code = 0 
band.1.beam = 1 2 
band.1.beam.model = PS 
band.1.beam.1.apriori = 0.01 0.239 780.0 1.0 
band.1.beam.1.sigma = 0.1 0.01 1.0 0.0 
band.1.beam.2.apriori = 0.01 0.18 780.0 1.0 
band.1.beam.2.sigma = 0.1 0.01 1.0 0.0 
band.1.gasb = CLONO2 CO2 O3 HNO3 
band.2.nu_start = 782.55 
band.2.nu_stop = 782.87 
band.2.zshift = F 
band.2.beam = 1 
band.2.beam.model = PS 
band.2.beam.1.apriori = 0.01 0.239 780.0 1.0 
band.2.beam.1.sigma = 0.1 0.01 0.1 0.0 
band.2.calc_point_space = 0.500E-03 
band.2.wave_factor = 0.999998990000 
band.2.max_opd = 180.0000 
band.2.omega = 2.3923 
band.2.apodization_code = 0 
band.2.gasb = O3 CO2 H2O HNO3 
kb.line = T 
kb.line.type = 1 
kb.line.gas = CLONO2 
kb = T 
kb.slope = T 
kb.curvature = T 
kb.temperature = T 
kb.phase = T 
kb.omega = T 
kb.zshift = T 
kb.sza = T 
out.level = 1 
fw.tips = F 
