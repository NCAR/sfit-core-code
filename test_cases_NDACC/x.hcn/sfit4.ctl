file.in.stalayers = station.layers 
file.in.refprofile = reference.prf 
file.in.spectrum = t15asc.4 
file.in.modulation_fcn = Bruker_EAP_simple.dat 
file.in.phase_fcn = Bruker_EAP_simple.dat 
file.in.sa_matrix = sa.input 
file.in.isotope = isotope.input 
file.in.solarlines = ../../linelist/solar/120621/solar.dat 
file.in.linelist = 03263.959056-03291.520944.hbin 
file.out.pbpfile = pbpfile 
file.out.statevec = statevec 
file.out.k_matrix = k.output 
file.out.kb_matrix = kb.output 
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
gas.profile.list = HCN O3
gas.profile.HCN.correlation = F
gas.profile.HCN.logstate = F 
gas.profile.HCN.scale = 1.0 
gas.profile.HCN.sigma = 0.1336  0.1474  0.1640  0.1746  0.1860  0.1981  0.2109  0.2252  0.2411  0.2585 0.2778  0.2988  0.3208  0.3442  0.3666  0.3858  0.4016  0.4110  0.4167  0.4211 0.4272  0.4319  0.4385  0.4454  0.4508  0.4583  0.4642  0.4725  0.4811  0.4880 0.4975  0.5051  0.5160  0.5276  0.5376  0.5508  0.5618  0.5774  0.5938  0.6081 0.6275  0.6439  0.6670  0.6940  0.7157  0.7487  0.7866 
gas.profile.O3.correlation = F 
gas.profile.O3.logstate = F 
gas.profile.O3.scale = 1.0 
gas.profile.O3.sigma = 0.1336  0.1474  0.1640  0.1746  0.1860  0.1981  0.2109  0.2252  0.2411  0.2585 0.2778  0.2988  0.3208  0.3442  0.3666  0.3858  0.4016  0.4110  0.4167  0.4211 0.4272  0.4319  0.4385  0.4454  0.4508  0.4583  0.4642  0.4725  0.4811  0.4880 0.4975  0.5051  0.5160  0.5276  0.5376  0.5508  0.5618  0.5774  0.5938  0.6081 0.6275  0.6439  0.6670  0.6940  0.7157  0.7487  0.7866 
gas.column.list = H2O C2H2 CO2 CH4
gas.column.H2O.scale = 1.0 
gas.column.H2O.sigma = 1.0 
gas.column.CH4.scale = 1.0 
gas.column.CH4.sigma = 1.0 
gas.column.C2H2.scale = 1.0 
gas.column.C2H2.sigma = 1.0 
gas.column.CO2.scale = 1.0 
gas.column.CO2.sigma = 1.0 
fw.raytonly = F 
fw.isotope_separation = F 
fw.delnu = 0.1000 
fw.lshapemodel = 0 
fw.linemixing = T 
fw.linemixing.gas = CO2 
fw.solar_spectrum = T 
fw.pressure_shift = T
fw.apod_fcn = F 
fw.apod_fcn.type = 2 
fw.apod_fcn.order = 1 
fw.phase_fcn = F 
fw.phase_fcn.type = 4 
fw.phase_fcn.order = 0 
fw.emission = F 
rt = T 
rt.lm = F 
rt.lm.gamma_start = 1.0e5 
rt.lm.gamma_inc = 10.0 
rt.lm.gamma_dec = 10.0 
rt.convergence = 0.01 
rt.tolerance = 0.1 
rt.max_iteration = 17 
rt.dwshift = F 
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
rt.phase = F
rt.phase.apriori = 0.000 
rt.phase.sigma = 0.200 
rt.apod_fcn = F 
rt.apod_fcn.apriori = 1.000 
rt.apod_fcn.sigma = 0.200 
rt.phase_fcn = F 
rt.phase_fcn.apriori = 1.000 
rt.phase_fcn.sigma = 0.200 
rt.solshift = F
rt.solshift.apriori = 1.000 
rt.solshift.sigma = 0.500 
rt.solstrnth = F 
rt.solstrnth.apriori = 1.000 
rt.solstrnth.sigma = 0.500 
rt.temperature = F 
kb = T 
kb.temperature = T 
kb.slope = F 
kb.curvature = F 
kb.solshft = F 
kb.solstrnth = F 
kb.phase = F 
kb.wshift = F 
kb.apod_fcn = F 
kb.phase_fcn = F 
kb.zshift = F 
kb.sza = T 
kb.omega = F 
kb.max_opd = F 
kb.line = F 
kb.line.type = 1 
kb.line.gas = C2H6 
band = 1 2
band.1.nu_start = 3268.0000 
band.1.nu_stop = 3268.3800 
band.1.zshift = F 
band.1.zshift.type = 1 
band.1.zshift.apriori = 0.000 
band.1.zshift.sigma = 0.200 
band.1.calc_point_space = 0.0005 
band.1.wave_factor = 1.000 
band.1.max_opd = 257.14 
band.1.omega = 2.75 
band.1.apodization_code = 0 
band.1.gasb = HCN H2O C2H2 O3 CO2 CH4
band.1.tempretb = F 
band.2.nu_start = 3287.0000 
band.2.nu_stop = 3287.4800 
band.2.zshift = F 
band.2.zshift.type = 1 
band.2.zshift.apriori = 0.000 
band.2.zshift.sigma = 0.200 
band.2.calc_point_space = 0.0005 
band.2.wave_factor = 1.000 
band.2.max_opd = 257.14 
band.2.omega = 2.75 
band.2.apodization_code = 0 
band.2.gasb = HCN H2O C2H2 CO2 O3
band.2.tempretb = F 
sp.snr = 
sp.snr.1.nu_start = 3268.0000 
sp.snr.1.nu_stop = 3268.3800 
sp.snr.1.snr = 200.0 
sp.snr.2.nu_start = 3287.0000 
sp.snr.2.nu_stop = 3287.4800 
sp.snr.2.snr = 200.0 
out.level = 1 
out.gas_spectra = T
out.gas_spectra.type = 1 
out.sa_matrix = T 
out.statevec = T 
out.g_matrix = T 
out.k_matrix = T 
out.shat_matrix = F 
out.retprofiles = T 
out.aprprofiles = T 
out.ab_matrix = F 
out.ak_matrix = T 
out.summary = T 
out.pbpfile = T 
out.channel = F 
out.parm_vectors = F 
out.seinv_vector = F 
out.sainv_matrix = F 
out.smeas_matrix = F 
out.ssmooth_matrix = F 
out.raytrace = F 
out.raytrace.type = 1 
out.solarspectrum = F 
out.levmardet = F 
out.xscdetail = F 
