 file.in.stalayers = station.layers
 file.in.refprofile = reference.prf
 file.in.spectrum = spectrum
 file.in.solarlines= solar.dat
 file.in.linelist= 03263.959056-03291.520944.hbin
 
 gas.layers = 48
 gas.profile.list = HCN H2O
 gas.profile.HCN.correlation = F
 gas.profile.HCN.logstate = F
 gas.profile.HCN.scale = 1.0
 gas.profile.HCN.sigma =
 0.2109  0.2252  0.2411  0.2585 0.2778  0.2988
 0.3208  0.3442  0.3666  0.3858 0.4016  0.4110
 0.4167  0.4211  0.4272  0.4319 0.4385  0.4454
 0.4508  0.4583  0.4642  0.4725 0.4811  0.4880
 0.4975  0.5051  0.5160  0.5276 0.5376  0.5508
 0.5618  0.5774  0.5938  0.6081 0.6275  0.6439
 0.6670  0.6940  0.7157  0.7487 0.7866  0.7866
 0.7866  0.7866  0.7866  0.7866 0.7866  0.7866
 gas.profile.H2O.correlation = F
 gas.profile.H2O.logstate = F
 gas.profile.H2O.scale = 1.0
 gas.profile.H2O.sigma =
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1
 gas.column.list =  O3 C2H2 CH4
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
 rt.curvature = T
 rt.curvature.apriori = 0.000
 rt.curvature.sigma = 0.100
 rt.phase = F
 rt.solshift = T
 rt.solshift.apriori = 0.000
 rt.solshift.sigma = 0.500
 rt.solstrnth = F
 rt.solstrnth.apriori = 1.000
 rt.solstrnth.sigma = 0.500
 rt.temperature = F
 
 kb= T
 kb.temperature= T
 kb.slope= F
 kb.curvature= F
 kb.solshft= T
 kb.solstrnth= T
 kb.phase= T
 kb.wshift = F
 kb.apod_fcn = F
 kb.phase_fcn = F
 kb.zshift= T
 kb.sza= T
 kb.omega= T
 kb.line = T
 kb.line.type = 1
 kb.line.gas = HCN
 
 band = 1 2
 band.1.nu_start = 3268.0000
 band.1.nu_stop = 3268.3800
 band.1.zshift = F
 band.1.zshift.type = 1
 band.1.zshift.apriori = 0.000
 band.1.zshift.sigma = 0.200
 band.1.calc_point_space = 0.0005
 band.1.wave_factor = 1.000
 band.1.max_opd= 180.0000
 band.1.omega= 2.7512
 band.1.apodization_code = 0
 band.1.gasb = HCN H2O C2H2 O3
 band.1.tempretb = F
 
 band.2.nu_start = 3287.0000
 band.2.nu_stop = 3287.4800
 band.2.zshift = F
 band.2.zshift.type = 1
 band.2.zshift.apriori = 0.000
 band.2.zshift.sigma = 0.200
 band.2.calc_point_space = 0.0005
 band.2.wave_factor = 1.000
 band.2.max_opd= 180.0000
 band.2.omega= 2.7512
 band.2.apodization_code = 0
 band.2.gasb = HCN H2O C2H2 O3 CH4
 band.2.tempretb = F
 
 out.level= 1
 out.gas_spectra= F
 out.k_matrix= T
 out.ak_matrix= T
 out.g_matrix= T
 out.sa_matrix= T
 out.shat_matrix= T
 out.seinv_vector= T
 out.retprofiles= T
 out.aprprofiles= T
 out.pbpfile= T
 out.statevec= T
 file.out.pbpfile= pbpfile
 file.out.statevec= statevec
 file.out.ak_matrix= ak.out
 file.out.k_matrix= k.out
 file.out.g_matrix= g.out
 file.out.kb_matrix= kb.out
 file.out.sa_matrix= sa.out
 file.out.shat_matrix= shat.out
 file.out.seinv_vector= seinv.out
 file.out.retprofiles= rprfs.table
 file.out.aprprofiles= aprfs.table
 file.out.summary= summary
