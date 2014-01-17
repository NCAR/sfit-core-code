file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = spectrum
file.in.modulation_fcn = apod_poly.dat
file.in.phase_fcn = phase_poly.dat
file.in.sa_matrix = sainv.input 
file.in.isotope = isotope.input 
file.in.solarlines = ../../linelist/solar/120621/solar.dat 
file.in.linelist = 02609.659356-02908.070644.hbin
file.out.pbpfile = pbpfile 
file.out.statevec = statevec 
file.out.k_matrix = k.output 
file.out.g_matrix = g.out
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
gas.layers = 47 
gas.profile.list = CH4 
gas.profile.CH4.correlation = T 
gas.profile.CH4.correlation.type = 5 
gas.profile.CH4.logstate = F 
gas.profile.CH4.scale = 1.0d0 
gas.profile.CH4.sigma = 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000    1.000000 1.000000    1.000000    1.000000    1.000000    1.000000 
gas.column.list = CO2 HDO NO2
gas.column.H2O.scale = 0.2d0 
gas.column.H2O.sigma = 0.1d0 
gas.column.CO2.scale = 1.0d0 
gas.column.CO2.sigma = 0.1d0 
gas.column.HDO.scale = 0.2d0 
gas.column.HDO.sigma = 0.01d0 
gas.column.NO2.scale = 1.0d0 
gas.column.NO2.sigma = 0.1d0 
fw.raytonly = F 
fw.isotope_separation = T 
fw.delnu = 0.100 
fw.lshapemodel = 0 
fw.linemixing = T 
fw.linemixing.gas = CO2 
fw.solar_spectrum = T 
fw.pressure_shift = T 
fw.apod_fcn =  F
fw.apod_fcn.type = 2
fw.apod_fcn.order = 2 
fw.phase_fcn = F 
fw.phase_fcn.type = 2 
fw.phase_fcn.order = 2 
fw.emission = F 
rt = T 
rt.lm = F 
rt.lm.gamma_start = 1.0e5 
rt.lm.gamma_inc = 10.0 
rt.lm.gamma_dec = 10.0 
rt.convergence = 0.01 
rt.tolerance = 0.050 
rt.max_iteration = 30 
rt.dwshift = F 
rt.wshift = T 
rt.wshift.type = 3 
rt.wshift.apriori = 0.000 
rt.wshift.sigma = 0.100 
rt.slope = T
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
rt.phase = F 
rt.phase.apriori = 0.000 
rt.phase.sigma = 0.100 
rt.apod_fcn = F
rt.apod_fcn.apriori = 0.0 
rt.apod_fcn.sigma = 0.1 
rt.phase_fcn = F 
rt.phase_fcn.apriori = 0.0 
rt.phase_fcn.sigma = 0.1 
rt.solshift = F 
rt.solshift.apriori = 0.000 
rt.solshift.sigma = 0.500 
rt.solstrnth = F 
rt.solstrnth.apriori = 0.000 
rt.solstrnth.sigma = 0.500 
rt.temperature = F 
kb = T 
kb.temperature = T
kb.slope = T
kb.curvature = T 
kb.solshft = T
kb.solstrnth = T 
<<<<<<< HEAD
kb.phase = T
kb.dwshift = F 
kb.wshift = F 
kb.apod_fcn = F 
kb.phase_fcn = F 
kb.zshift = T
kb.sza = T
kb.omega = T 
kb.max_opd = T 
kb.line = F 
kb.slope = T
kb.curvature = T 
kb.solshft = T
kb.solstrnth = T 
kb.phase = T
=======
kb.phase = F
>>>>>>> BugFix_Mathias
kb.wshift = T
kb.dwshift = T
kb.apod_fcn = T
kb.phase_fcn = T 
kb.zshift = F 
kb.sza = T
kb.omega = T 
kb.line = T
kb.line.type = 1 
<<<<<<< HEAD
kb.line.gas = CH4 
=======
kb.line.gas =  retrieval
>>>>>>> BugFix_Mathias
kb.profile = T
kb.profile.gas = CO2 
band = 1 2 3 4 
band.1.nu_start = 2613.7000 
band.1.nu_stop = 2615.4000 
band.1.zshift = F 
band.1.calc_point_space = 0.900E-03 
band.1.wave_factor = 1.000 
band.1.max_opd = 257.143 
band.1.omega = 2.2727 
band.1.apodization_code = 0 
band.1.beam = 0 
band.1.gasb = CH4 HDO CO2 
band.1.tempretb = F 
band.2.nu_start = 2650.6000 
band.2.nu_stop = 2651.3000 
band.2.zshift = F 
band.2.calc_point_space = 0.900E-03 
band.2.wave_factor = 1.000 
band.2.max_opd = 257.143 
band.2.omega = 2.2727 
band.2.apodization_code = 0 
band.2.beam = 0 
band.2.gasb = CH4 HDO CO2 
band.2.tempretb = F 
band.3.nu_start = 2835.5000 
band.3.nu_stop = 2835.8000 
band.3.zshift = F 
band.3.calc_point_space = 0.900E-03 
band.3.wave_factor = 1.000 
band.3.max_opd = 257.143 
band.3.omega = 2.2727 
band.3.apodization_code = 0 
band.3.beam = 0 
band.3.gasb = CH4 
band.3.tempretb = F 
band.4.nu_start = 2903.600 
band.4.nu_stop = 2904.0300 
band.4.zshift = F 
band.4.calc_point_space = 0.900E-03 
band.4.wave_factor = 1.000 
band.4.max_opd = 257.143 
band.4.omega = 2.2727 
band.4.apodization_code = 0 
band.4.beam = 0 
band.4.gasb = CH4 NO2 
band.4.tempretb = F 
sp.snr.1.nu_start = 2600. 
sp.snr.1.nu_stop = 2650.0 
sp.snr.1.snr = 200.0 
sp.snr.2.nu_start = 0 
sp.snr.2.nu_stop = 0 
sp.snr.2.snr = 200.0 
out.level = 1 
out.gas_spectra = F 
