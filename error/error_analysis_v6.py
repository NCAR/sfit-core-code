"""
error_analysis_v6.py

Python script which calculates systematic and random uncertainty covariance matrix 
using output from sfit4 and sb values from ctl file

Outputs:
error_analysis.log: summary of calculated erros in %
error_analysis.summary: summary of calculated erros in molec / cm^2
outputs for covariance matrices are controlled by sfit4.ctl file
options:
  total systematic uncertainty covariance matrix
  total random uncertainty covariance matrix
  all systematic uncertainty covariance matrices
  all random uncertainty covariance matrices

Updates 09/07/2013

1) Made changes to reconcile changes to version 2 of file read
  a) n_profile and n_column are now calculated in file_read
2) Added vmr covariance matrix outpur
3) Changed error summary to error log
4) Added error summary of uncertainties in mol cm^-2
5) Fixed size of AK_int2 and Se_int2 for case of interfering species retrieved as profile
6) Improved outputs, allows choice of untis of covariance matrix between vmr and molec /cm 2

"""
import numpy as np
import file_read_v4 as rout
import read_prfs as prfs
import sfit4_ctl_simple as sfit4

# Read in control files
ctl = sfit4.sfit4_ctl()
ctl.read_ctl_file('sfit4.ctl')
b = sfit4.sfit4_ctl()
b.read_ctl_file('sb.ctl')


# Read in outfiles files from sfit4: 

# k.output to calculate averaging kernel
# d.complete to calculate averaging kernel and measurement error
# ssmooth.target if present or sa.complete to calculate smoothing error
# summary file for the SNR to calculate the measurement error
# kb.output for non-retrieved parameter error calculation
# rprfs.table and aprfs.table for the airmass, apriori and retrieved profiles 

a = rout.read_output_for_errorcalc(1)

# If all files that are required do not exist, STOP

if (a.filemissing):
  exit()

# Primary retrieved gas is assumed to be the first gas in the profile gas list.  
# If no gases are retrieved as a profile, the primary gas is assumed to be the first gas in the column gas list.

if (a.n_profile > 0):
  primgas = a.gas_profile[0]
else:
  primgas = a.gas_column[0]

r = prfs.prfs(primgas)

# Gain matrix for the retrieved profile of the gas in questions only,
Dx = a.D[a.x_start:a.x_stop,:]
# Calculate the unscaled Averaging Kernel: AK = D*K
AK = np.dot(a.D,a.K)

# Unscaled Averaging Kernel for the retrieved profile of the gas in questions only,
AKx = AK[a.x_start:a.x_stop,a.x_start:a.x_stop]

# Scaled Averaging Kernel
AKvmr = np.dot(np.dot(np.diag(1/r.aprf),AKx),np.diag(r.aprf))

# Retrieved total column
retdenscol = np.sum(r.rprf*r.airmass)

# A priori density profile
aprdensprf = r.aprf*r.airmass

# DOFS for column
col_dofs = np.trace(AKx)

# List of all parameters 
Kb_labels = ['temperature','solshft','solstrnth','phase','wshift','dwshift','sza','lineInt','lineTAir','linePAir','slope','curvature','apod_fcn','phase_fcn','omega','max_opd','zshift']

# List of all calculated random error covariance matrices, percent errors and labels
S_ran = []

# List of all calculated systematic error covariance matrices, percent errors and labels 
S_sys = [] 

try:
  filename = b.value['file.out.error.log']
except KeyError, e:
  filename = 'error_analysis.log'

f = open(filename,'w')

# Calculate Smoothing error

if (a.Ssmooth_file=='green'):
  # Use Ssmooth matrix output by sfit4 if available.  
  # Ss_err = r.aprf density profile * Ssmooth * transpose(r.aprf density profile)
  S_tmp = a.Ssmooth
  S_err_tmp = np.sqrt(np.dot(np.dot(aprdensprf,a.Ssmooth),aprdensprf.T))
  S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),a.Ssmooth),np.diag(aprdensprf))
else:
  # Otherwise calculate from Rogers: Ss=(A-I)*Sa*(A-I)T
  S_tmp = np.dot(np.dot((AKvmr - np.identity(AKvmr.shape[0])),a.sa[a.x_start:a.x_stop,a.x_start:a.x_stop]),(AKvmr - np.identity(AKvmr.shape[0])).T)
  S_err_tmp = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp),aprdensprf.T))
  S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))

S_sys.append((S_tmp,S_tmp2,S_err_tmp,'smoothing'))

# Calculate Measurement error using SNR calculated from the spectrum

# Calculate from Rogers: Sm = Dx * Sm * Dx.T
S_tmp = np.dot(np.dot(Dx,a.se),Dx.T)
S_err_tmp = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp),aprdensprf.T))
S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))

S_ran.append((S_tmp,S_tmp2,S_err_tmp,'measurement'))
f.write('%s     %.4E\n'%('Column amount = ', retdenscol))
f.write('%s     %f\n'%('DOFS (total column) = ', col_dofs))
f.write('%s     %f\n'%('Sm (%) = ', (S_ran[0][2]/retdenscol)*100))
f.write('%s     %f\n'%('Ss (%) = ', (S_sys[0][2]/retdenscol)*100))

# Interference error 1 - retrieval parameters
AK_int1 = AK[0:a.x_start,a.x_start:a.x_stop].T
Se_int1 =  a.sa[0:a.x_start,0:a.x_start]
  
S_tmp = np.dot(np.dot(AK_int1,Se_int1),AK_int1.T)
S_err_tmp = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp),aprdensprf.T)) 
S_tmp2 =  np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))

f.write('%s     %f\n'%('Sint1_random (retrieval params) (%) = ', S_err_tmp/retdenscol*100))

S_ran.append((S_tmp,S_tmp2,S_err_tmp,'retrieval parameters'))    

# Interference error 2 - interfering species

n_int2 = a.n_profile + a.n_column - 1
n_int2_column = (a.n_profile-1)*a.n_layer + a.n_column

AK_int2 = AK[a.x_stop:a.x_stop+n_int2_column,a.x_start:a.x_stop].T

Se_int2 = a.sa[a.x_stop:a.x_stop+n_int2_column,a.x_stop:a.x_stop+n_int2_column]

S_tmp = np.dot(np.dot(AK_int2,Se_int2),AK_int2.T)
S_err_tmp = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp),aprdensprf.T))
S_tmp_2 =  np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))

S_ran.append((S_tmp,S_tmp_2,S_err_tmp,'interfering species'))

f.write('%s     %f\n'%('Sint2_random (intf. spec.) (%) = ', (S_err_tmp/retdenscol)*100))

# Errors from parameters not retrieved by Sfit4

zerolev_band_b = [];
for k in ctl.value['band'].strip().split():
  if (ctl.value['band.'+k+'.zshift'].strip() == 'F'):  # only include bands where zshift is retrieved
    zerolev_band_b.append(k)

for s in Kb_labels:
  Kb_exists = True
  try:
    DK = np.dot(Dx,a.Kb[s])
  except KeyError, e:
    Kb_exists = False
  
  if (Kb_exists):
    for k in ['random', 'systematic']:
#      print s
      try:
        if (s == 'zshift'):    
          Sb = np.zeros((len(zerolev_band_b),len(zerolev_band_b)))
          for i in range(0,len(zerolev_band_b)):
            Sb[i,i] = float( b.value['sb.band.'+zerolev_band_b[i]+'.zshift.'+k])**2
        elif (DK.shape[1] == 1):
          Sb = np.zeros((1,1))
          Sb[0,0] = float(b.value['sb.'+s+'.'+k])**2
        elif (DK.shape[1] == a.n_window):
          Sb = np.dot(float(b.value['sb.'+s+'.'+k])**2,np.identity(a.n_window))
        elif (DK.shape[1] == a.n_layer):  # Only the temperature error fits this requirement
          temp_sb = b.value['sb.temperature.'+k].strip().split()
          Sb = np.zeros((a.n_layer,a.n_layer))
          for j in range(0, len(temp_sb)):
            Sb[j,j] = (float(temp_sb[j])/r.T[j])**2  # convert degrees to relative units
          del temp_sb
        else:
           Sb = 0
      except KeyError, e:
        Sb = 0
      if np.sum(Sb) == 0:
        if (s == 'zshift'):
          f.write('%s\n'%('sb.band.x.'+s+'.'+k+' for all bands where zshift is not retrieved is 0 or is not specified, error covariance matrix not calculated'))
        else:
          f.write('%s\n'%('sb.'+s+'.'+k+' for '+s+' is 0 or is not specified, error covariance matrix not calculated'))
      else:
        S_tmp = np.dot(np.dot(DK,Sb),DK.T)
        S_err_tmp = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp),aprdensprf.T))
        S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))
        f.write('%s     %f\n'%('S.'+s+'.'+k+' (%) =', (S_err_tmp/retdenscol)*100))

        if (k == 'random'):
          S_ran.append((S_tmp,S_tmp2,S_err_tmp,s))

        if (k == 'systematic'):
          S_sys.append((S_tmp,S_tmp2,S_err_tmp,s))

S_random_err = 0
S_random = np.zeros((a.n_layer,a.n_layer))
S_random_vmr = np.zeros((a.n_layer,a.n_layer))

for k in range(0,len(S_ran)):
  S_random_err += S_ran[k][2]**2
  S_random_vmr += S_ran[k][0]
  S_random += S_ran[k][1]
S_random_err = np.sqrt(S_random_err)

f.write('%s     %f\n'%('Total random error (%) = ', (S_random_err/retdenscol)*100))

S_systematic_err = 0
S_systematic =  np.zeros((a.n_layer,a.n_layer))
S_systematic_vmr =  np.zeros((a.n_layer,a.n_layer))

for k in range(1,len(S_sys)):
  S_systematic_err += S_sys[k][2]**2
  S_systematic_vmr += S_sys[k][0]
  S_systematic += S_sys[k][1]
S_systematic_err = np.sqrt(S_systematic_err)
f.write('%s     %f\n'%('Total systematic error (%) = ', (S_systematic_err/retdenscol)*100))

f.write('%s     %f\n'%('Total error exclusive of smoothing error (%) = ',np.sqrt(S_systematic_err**2 + S_random_err**2)/retdenscol*100))
f.close()

try:
  filename = b.value['file.out.error.summary']
except KeyError, e:
  filename = 'error_analysis.summary'

with open(filename, 'w') as fout:
  fout.write(a.header + 'ERROR SUMMARY\n\n')
  fout.write('  %s  %.4e\n'%('Total random uncertainty in mol cm^-2 :', S_random_err))
  fout.write('  %s  %.4e\n'%('Total systematic uncertainty in mol cm^-2 :', S_systematic_err))
  for k in range(0, len(S_ran)):
    fout.write('  %s%s %s  %.4e\n'%(S_ran[k][3][0].upper(),S_ran[k][3][1:],'random uncertainty in mol cm^-2 :', S_ran[k][2]))
  for k in range(0, len(S_sys)):
    fout.write('  %s%s %s  %.4e\n'%(S_sys[k][3][0].upper(),S_sys[k][3][1:],'systematic uncertainty in mol cm^-2 :', S_sys[k][2]))

# Output Covariance matrices in mol/cm^2

try:
  output = b.value['out.ssystematic'].strip()
except KeyError, e:
  output = False

if (output == 'T'):
  try:
    filename = b.value['file.out.ssystematic']
  except KeyError, e:
    filename = 'ssystematic.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
    fout.write('          ' + str(S_systematic.shape[0]) + '          ' + str(S_systematic.shape[1]) +'\n')
    for row in  S_systematic:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom'].strip()
except KeyError, e:
  output = False
if (output == 'T'):

  try:
    filename = b.value['file.out.srandom']
  except KeyError, e:
    filename = 'srandom.output'

  with open(filename, 'w') as fout: 
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
    fout.write('          ' + str(S_random.shape[0]) + '          ' + str(S_random.shape[1])+'\n')
    for row in S_random:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.all'].strip()
except KeyError, e:
  output = False

if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.all']
  except KeyError, e:
    filename = 'srandom.all.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRICES IN MOL CM^-2\n')
    fout.write('          ' + str(S_ran[0][1].shape[0]) + '          ' + str(S_ran[0][1].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_ran)):
      fout.write('  '+S_ran[k][3].upper()+' ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
      for row in  S_ran[k][1]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

try:
  output = b.value['out.ssystematic.all'].strip()
except KeyError, e:
  output = False

if (output == 'T'): 
  try:
    filename = b.value['file.out.ssytematic.all'].strip()
  except KeyError, e:
    filename = 'ssystematic.all.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRICES IN MOL CM^-2\n')
    fout.write('          ' + str(S_sys[0][1].shape[0]) + '          ' + str(S_sys[0][1].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_sys)):
      fout.write('  '+S_sys[k][3].upper()+' ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
      for row in  S_sys[k][1]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

# Output Covariance matrices in vmr

try:
  output = b.value['out.ssystematic.vmr'].strip()
except KeyError, e:
  output = False

if (output == 'T'):
  try:
    filename = b.value['file.out.ssystematic.vmr']
  except KeyError, e:
    filename = 'ssystematic.vmr.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRIX IN VMR UNITS\n')
    fout.write('          ' + str(S_systematic_vmr.shape[0]) + '          ' + str(S_systematic_vmr.shape[1]) +'\n')
    for row in  S_systematic_vmr:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.vmr'].strip()
except KeyError, e:
  output = False
if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.vmr']
  except KeyError, e:
    filename = 'srandom.vmr.output'

  with open(filename, 'w') as fout: 
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRIX IN VMR UNITS \n')
    fout.write('          ' + str(S_random_vmr.shape[0]) + '          ' + str(S_random_vmr.shape[1])+'\n')
    for row in S_random_vmr:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.vmr.all'].strip()
except KeyError, e:
  output = False
if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.vmr.all']
  except KeyError, e:
    filename = 'srandom.vmr.all.output'
  with open(filename, 'w') as fout:
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRICES IN VMR UNITS \n')
    fout.write('          ' + str(S_ran[0][0].shape[0]) + '          ' + str(S_ran[0][0].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_ran)):
      fout.write('  '+S_ran[k][3].upper()+' ERROR COVARIANCE MATRIX IN VMR UNITS \n')
      for row in S_ran[k][0]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

try:
  output = b.value['out.ssystematic.vmr.all'].strip()
except KeyError, e:
  output = False

if (output == 'T'): 
  try:
    filename = b.value['file.out.ssytematic.vmr.all'].strip()
  except KeyError, e:
    filename = 'ssystematic.vmr.all.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRICES IN VMR UNITS \n')
    fout.write('          ' + str(S_sys[0][0].shape[0]) + '          ' + str(S_sys[0][0].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_sys_vmr)):
      fout.write('  '+S_sys[k][3].upper()+' ERROR COVARIANCE MATRIX IN VMR UNITS \n')
      for row in  S_sys[k][0]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

