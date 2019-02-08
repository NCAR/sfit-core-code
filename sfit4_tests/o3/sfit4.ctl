file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = spectrum 
file.in.isotope = isotope.input 
file.in.modulation_fcn = ils.dat 
file.in.phase_fcn = ils.dat 
file.in.solarlines = solar.dat 
file.in.linelist = 00778.501722-01008.558278.hbin 
gas.layers = 48 
gas.profile.list = O3 H2O 
gas.column.list = CO2 O3668 O3686 C2H4 
gas.profile.O3.correlation = F 
gas.profile.O3.logstate = F 
gas.profile.O3.scale = 1.0 
gas.profile.O3.sigma = 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 0.1   0.1   0.1   0.1   0.1 
gas.profile.H2O.correlation = F 
gas.profile.H2O.logstate = F 
gas.profile.H2O.scale = 1.0 
gas.profile.H2O.sigma = 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 1.0   1.0   1.0   1.0   1.0 
gas.column.O3668.scale = 1.0 
gas.column.O3668.sigma = 1.0 
gas.column.O3686.scale = 1.0 
gas.column.O3686.sigma = 1.0 
gas.column.CO2.scale = 1.0 
gas.column.CO2.sigma = 1.0 
gas.column.C2H4.scale = 1.0 
gas.column.C2H4.sigma = 1.0 
fw.delnu = 0.1 
fw.lshapemodel = 0 
fw.solar_spectrum = T 
fw.pressure_shift = T 
fw.apod_fcn = F 
fw.phase_fcn = F 
fw.isotope_separation = T 
rt = T 
rt.lm = F 
rt.convergence = 0.1 
rt.max_iteration = 15 
rt.wshift = T 
rt.wshift.type = 3 
rt.wshift.apriori = 0.000 
rt.wshift.sigma = 0.100 
rt.slope = T 
rt.slope.apriori = 0.000 
rt.slope.sigma = 0.100 
rt.curvature = T 
rt.curvature.apriori = 0.000 
rt.curvature.sigma = 0.100 
rt.phase = T 
rt.phase.apriori = 0.000 
rt.phase.sigma = 0.200 
rt.apod_fcn = F 
rt.apod_fcn.apriori = 0.0 
rt.apod_fcn.sigma = 0.1 
rt.phase_fcn = F 
rt.dwshift = F 
band = 1 2 3 4 
band.1.nu_start = 782.56 
band.1.nu_stop = 782.86 
band.1.zshift = F 
band.1.beam = 0 
band.1.calc_point_space = 0.500E-03 
band.1.wave_factor = 0.999998990000 
band.1.max_opd = 180.0000 
band.1.omega = 3.5885 
band.1.apodization_code = 0 
band.1.gasb = O3  H2O    CO2 O3668 O3686 
band.2.nu_start = 788.85 
band.2.nu_stop = 789.37 
band.2.zshift = F 
band.2.beam = 0 
band.2.calc_point_space = 0.500E-03 
band.2.wave_factor = 0.999998990000 
band.2.max_opd = 180.0000 
band.2.omega = 3.5885 
band.2.apodization_code = 0 
band.2.gasb = O3     H2O    CO2 O3668 O3686 
band.3.nu_start = 993.3 
band.3.nu_stop = 993.8 
band.3.zshift = F 
band.3.beam = 0 
band.3.calc_point_space = 0.500E-03 
band.3.wave_factor = 0.999998990000 
band.3.max_opd = 180.0000 
band.3.omega = 3.5885 
band.3.apodization_code = 0 
band.3.gasb = O3     H2O    CO2 C2H4 O3668 O3686 
band.4.nu_start = 1000.000 
band.4.nu_stop = 1004.500 
band.4.zshift = T 
band.4.zshift.type = 1 
band.4.zshift.apriori = 0.000 
band.4.zshift.sigma = 0.200 
band.4.beam = 0 
band.4.calc_point_space = 0.500E-03 
band.4.wave_factor = 0.999998990000 
band.4.max_opd = 180.0000 
band.4.omega = 3.5885 
band.4.apodization_code = 0 
band.4.gasb = O3     H2O    CO2 O3668 O3686 
sp.snr = 1 2 3 
sp.snr.1.nu_start = 1000.00 
sp.snr.1.nu_stop = 1000.8 
sp.snr.1.snr = 1.0 
sp.snr.2.nu_start = 1001.0 
sp.snr.2.nu_stop = 1001.3 
sp.snr.2.snr = 1.0 
sp.snr.3.nu_start = 1003.16 
sp.snr.3.nu_stop = 1004.5 
sp.snr.3.snr = 1.0 
kb = F 
kb.line = T 
kb.line.gas = target 
out.level = 1 
out.smeas_matrix = F 
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
kb.solshft = T 
kb.temperature = T 
kb.omega = T 
kb.zshift = T 
