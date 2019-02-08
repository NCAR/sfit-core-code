file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = spectrum 
file.in.solarlines = solar.dat 
file.in.linelist = 02971.941722-02982.058278.hbin 
gas.layers = 48 
gas.profile.list = C2H6 H2O CH4 
gas.profile.C2H6.correlation = F 
gas.profile.C2H6.logstate = F 
gas.profile.C2H6.scale = 1.0 
gas.profile.C2H6.sigma = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
gas.profile.H2O.correlation = F 
gas.profile.H2O.logstate = F 
gas.profile.H2O.scale = 1.0 
gas.profile.H2O.sigma = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
gas.profile.O3.correlation = F 
gas.profile.O3.logstate = F 
gas.profile.O3.scale = 1.0 
gas.profile.O3.sigma = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
gas.profile.CH4.correlation = F 
gas.profile.CH4.logstate = F 
gas.profile.CH4.scale = 1.0 
gas.profile.CH4.sigma = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
gas.column.list = O3 
gas.column.H2O.scale = 1.0 
gas.column.H2O.sigma = 1.0 
gas.column.C2H2.scale = 1.0 
gas.column.C2H2.sigma = 1.0 
gas.column.O3.scale = 1.0 
gas.column.O3.sigma = 1.0 
gas.column.CH4.scale = 1.0 
gas.column.CH4.sigma = 1.0 
fw.isotope_separation = F 
fw.delnu = 0.1000 
fw.lshapemodel = 0 
fw.linemixing = F 
fw.solar_spectrum = T 
fw.pressure_shift = T 
fw.emission = F 
rt = T 
rt.lm = F 
rt.convergence = 0.1 
rt.max_iteration = 17 
rt.dwshift = F 
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
rt.phase.sigma = 0.100 
rt.solshift = F 
rt.solshift.apriori = 0.000 
rt.solshift.sigma = 0.500 
rt.solstrnth = F 
rt.solstrnth.apriori = 1.000 
rt.solstrnth.sigma = 0.500 
rt.temperature = F 
kb = F 
kb.temperature = T 
kb.slope = T 
kb.curvature = T 
kb.solshft = T 
kb.solstrnth = T 
kb.phase = T 
kb.wshift = F 
kb.apod_fcn = F 
kb.phase_fcn = F 
kb.zshift = T 
kb.sza = T 
kb.omega = T 
kb.line = T 
kb.line.type = 1 
kb.line.gas = C2H6 
band = 1 
band.1.nu_start = 2976.0 
band.1.nu_stop = 2978.0 
band.1.zshift = F 
band.1.zshift.type = 1 
band.1.zshift.apriori = 0.000 
band.1.zshift.sigma = 0.200 
band.1.calc_point_space = 0.0005 
band.1.wave_factor = 1.000 
band.1.max_opd = 180.0000 
band.1.omega = 1.1962 
band.1.apodization_code = 0 
band.1.gasb = C2H6 H2O O3 CH4 
band.1.tempretb = F 
band.2.nu_start = 2982.6000 
band.2.nu_stop = 2984.500 
band.2.zshift = T 
band.2.zshift.type = 2 
band.2.zshift.apriori = 0.000 
band.2.zshift.sigma = 0.200 
band.2.calc_point_space = 0.0005 
band.2.wave_factor = 1.000 
band.2.max_opd = 257.1429 
band.2.omega = 2.2727 
band.2.apodization_code = 0 
band.2.gasb = C2H6 H2O O3 CH4 
band.2.tempretb = F 
out.level = 1 
out.gas_spectra = F 
out.k_matrix = T 
out.ak_matrix = T 
out.g_matrix = T 
out.sa_matrix = T 
out.shat_matrix = T 
out.seinv_vector = T 
out.retprofiles = T 
out.aprprofiles = T 
out.pbpfile = T 
out.statevec = T 
file.out.pbpfile = pbpfile 
file.out.statevec = statevec 
file.out.ak_matrix = ak.out 
file.out.k_matrix = k.out 
file.out.g_matrix = g.out 
file.out.kb_matrix = kb.out 
file.out.sa_matrix = sa.out 
file.out.shat_matrix = shat.out 
file.out.seinv_vector = seinv.out 
file.out.retprofiles = rprfs.table 
file.out.aprprofiles = aprfs.table 
file.out.summary = summary 
