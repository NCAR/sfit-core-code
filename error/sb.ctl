# Sb values for error analysis calculations

 sb.temperature.random                      =   
 9  9  9  9  9  9  9  7  7  7  7  6  6  5  5  2  2  2  2  2  2  2  2  2
 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

# sb.temperature.systematic                      =
# 9  9  9  9  9  9  9  7  7  7  7  6  6  5  5  2  2  2  2  2  2  2  2 2
# 2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2
# natural units here but should be relative, converted in error_cal.py
# sb.slope.random                     = 0.1
# sb.curvature.random                 = 0.1
# sb.solshft.random                   = 0.5
# sb.solstrnth.random                 = 0.1
 
# sb.phase.random                     = 0.0005 # natural units

# sb.wshift.random                    = 0.1

# sb.apod_fcn.random                  = 0.2
# sb.apod_fcn.systematic              = 0.2
# sb.phase_fcn.random                 = 0.2
# sb.phase_fcn.systematic             = 0.2

# sb.band.1.zshift.random             = 0.01
# sb.band.2.zshift.random             = 0.01
# sb.band.3.zshift.random       	     = 0.01
# sb.band.4.zshift.random       	     = 0.01
# sb.band.1.zshift.systematic         = 0.01
# sb.band.2.zshift.systematic   	     = 0.01
# sb.band.3.zshift.systematic   	     = 0.01
# sb.band.4.zshift.systematic   	     = 0.01

 sb.sza.random                       = 0.0056#0.125 absolute
# sb.sza.systematic                   = 0.001
# sb.omega.random                     = 0.001
# sb.omega.systematic                 = 0.001
# sb.max_opd.random                   = 0.005
# sb.max_opd.systematic               = 0.005
sb.lineInt.systematic               = 0.1 # relative
sb.lineTAir.systematic              = 0.1 # relative
sb.linePAir.systematic              = 0.1 # relative

sb.lineInt.HCN.systematic               = 0.1 # relative
sb.lineTAir.HCN.systematic              = 0.1 # relative
sb.linePAir.HCN.systematic              = 0.1 # relative

sb.lineInt.H2O.systematic               = 0.1 # relative
sb.lineTAir.H2O.systematic              = 0.1 # relative
sb.linePAir.H2O.systematic              = 0.1 # relative


# Output

 out.srandom                           = T
 out.srandom.all                       = F
 out.ssystematic.interferer.all        = T      
 out.ssystematic                       = T
 out.ssystematic.all                   = T
 out.srandom.interferer.vmr            = T
 out.ssystematic.interferer.vmr        = T
 out.srandom.vmr                       = T
# out.srandom.vmr.all                   = T
 out.ssystematic.vmr                   = T
# out.ssystematic.vmr.all               = F

# file names

 file.out.srandom                      = srandom.output
 file.out.srandom.all                  = srandom.all.output
 file.out.ssytematic                   = ssystematic.output
 file.out.ssytematic.all               = ssystematic.all.output
