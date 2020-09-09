import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
sys.path.append('/home/mathias/sfit-processing-environment/Lib_MP')
from read_result_sfit4 import pbp

#s_jn_wi = np.loadtxt('comp_justus/2-mit-10-cm-opd',skiprows=3)
#s_jn_wi = np.loadtxt('comp_justus/3-mit-10-cm-opd-und-1-mrad',skiprows=3)
s_jn_wo = np.loadtxt('comp_justus/4-ohne-cont',skiprows=3)
s_jn_wi = np.loadtxt('comp_justus/4-mit-10-cm-opd-1-mrad-alle-linien-aber-5-gase',skiprows=3)
#s_jn_wi = np.loadtxt('comp_justus/5-mit-allen-linien-und allen gasen-und-opd-und-fov',skiprows=3)
#s_jn_wi = np.loadtxt('comp_justus/table71-mit-cont-0.1', skiprows=3)
#s_jn_wo = np.loadtxt('comp_justus/table71-ohne-cont-0.1', skiprows=3)

#lbl1 = nc.Dataset('comp_lblrtm/clear_sky_mit_continuum.nc')
#lbl2 = nc.Dataset('comp_lblrtm/clear_sky_ohne_continuum.nc')
#s1_h2o = pbp('pbpfile_h2o_wi_mtckd')
#s2_h2o = pbp('pbpfile_h2o_wo_mtckd')

s1 = pbp('pbpfile_wi_mtckd')
s2 = pbp('pbpfile_wo_mtckd')

lbl1 = nc.Dataset('comp_lblrtm/spectral_radiance_all_continua.nc')
lbl2 = nc.Dataset('comp_lblrtm/spectral_radiance_no_continuum.nc')
s1_lbl = pbp('pbpfile_lbl_wi_mtckd')
s2_lbl = pbp('pbpfile_lbl_wo_mtckd')

f = plt.figure(1)

plt.subplot(141)
plt.plot(s_jn_wo[:,0], s_jn_wo[:,-1], label='jn ohne')
plt.plot(s_jn_wi[:,0], s_jn_wi[:,-1], label='jn mit')
plt.legend()
plt.subplot(142)
plt.plot(lbl2['wavenumber'], lbl2['spectral_radiance'])
plt.plot(lbl1['wavenumber'], lbl1['spectral_radiance'])
plt.legend()
plt.subplot(143)
plt.plot(s2.nu[0], s2.clc[0])
plt.plot(s1.nu[0], s1.clc[0])
plt.legend()
plt.subplot(144)
plt.plot(s_jn_wi[:,0], s_jn_wi[:,-1], label='jn mit')
plt.plot(s1.nu[0], s1.clc[0],label='sfit4 mit')
plt.legend()
f.show()

f1=plt.figure(2)
plt.subplot(121)
plt.plot(s_jn_wi[:,0], s_jn_wi[:,-1], label='jn mit')
plt.plot(s1.nu[0], s1.clc[0],label='sfit4 mit')
plt.legend()

plt.subplot(122)
plt.plot(s_jn_wo[:,0], s_jn_wo[:,-1], label='jn ohne')
plt.plot(s2.nu[0], s2.clc[0],label='sfit4 ohne')
plt.legend()
f1.show()

f2=plt.figure(3)
plt.subplot(121)
plt.plot(lbl1['wavenumber'], lbl1['spectral_radiance'][:]/1000, label='lbl mit')
plt.plot(s1_lbl.nu[0], s1_lbl.clc[0],label='sfit4 mit')
plt.legend()

plt.subplot(122)
plt.plot(lbl2['wavenumber'], lbl2['spectral_radiance'][:]/1000, label='lbl ohne')
plt.plot(s2_lbl.nu[0], s2_lbl.clc[0],label='sfit4 ohne')
plt.legend()

f2.show()
input()

