"""
error_analysis_with_interferer.py

Python script that calculates the systematic and random uncertainty covariance matrices 
for the main retrieval gas and one interfering species, if the interfering species if 
retrieved as a profile using output from sfit4 and sb values from ctl file

Outputs:
error_analysis.log: summary of calculated erros in %
error_analysis.summary: summary of calculated erros in molec / cm^2
outputs for covariance matrices are controlled by sfit4.ctl file
options:
  total systematic uncertainty covariance matrix
  total systematic uncertainty of interfering gas covariance matrix
  total random uncertainty covariance matrix
  total random uncertainty of interfering gas covariance matrix
  all systematic uncertainty covariance matrices
  all systematic uncertainty of interfering gas covariance matrices
  all random uncertainty covariance matrices
  all random uncertainty of interfering gas covariance matrices

Based on error_analysis_v6.py

Improvements that need to be made:
1) allow for multiple interfering species

"""
import numpy as np
import file_read_v4 as rout
import read_prfs as prfs
import sfit4_ctl_simple as sfit4

# Primary gas and one interferer: (use variable in case we ant to eventually allow for multiple inteferers)
ngas = 2

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

a = rout.read_output_for_errorcalc(ngas)

# If all files that are required do not exist, STOP

if (a.filemissing):
  exit()

# Primary retrieved gas is assumed to be the first gas in the profile gas list. 
# The interfering gas is assumed to be the second gas in the profile gas list.  

if (a.n_profile < 2):
 print 'Insufficient gases retrieved as profiles for error analysis of main and interfering species'
 exit()

primgas = a.gas_profile[0]
intgas = a.gas_profile[1]
r = prfs.prfs(primgas)
s = prfs.prfs(intgas)

# Gain matrix for the retrieved profile of the gas in questions only,
Dx = a.D[a.x_start:a.x_start+2*a.n_layer,:]
# Calculate the unscaled Averaging Kernel: AK = D*K
AK = np.dot(a.D,a.K)

# Unscaled Averaging Kernel for the retrieved profile of the gas in questions only,
AKx = AK[a.x_start:a.x_start+2*a.n_layer,a.x_start:a.x_start+2*a.n_layer]

# Retrieved total column
retdenscol = np.sum(r.rprf*r.airmass)
retdenscol_int = np.sum(s.rprf*s.airmass)

# A priori density profile
aprdensprf = r.aprf*r.airmass
aprdensprf_int = s.aprf*s.airmass

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

# Retrieval species: Calculate Smoothing error

# Scaled Averaging Kernel of primary gas
AKvmr = np.dot(np.dot(np.diag(1/r.aprf),AK[a.x_start:a.x_stop,a.x_start:a.x_stop]),np.diag(r.aprf))

# Calculate from Rogers: Ss=(A-I)*Sa*(A-I)T
S_tmp_1 = np.dot(np.dot((AKvmr - np.identity(AKvmr.shape[0])),a.sa[a.x_start:a.x_stop,a.x_start:a.x_stop]),(AKvmr - np.identity(AKvmr.shape[0])).T)
S_err_tmp_1 = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp_1),aprdensprf.T))
#S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))

# Calculate Smoothing error

# Scaled Averaging Kernel of Interferer
AKvmr = np.dot(np.dot(np.diag(1/s.aprf),AK[a.x_stop:a.x_stop+a.n_layer,a.x_stop:a.x_stop+a.n_layer]),np.diag(s.aprf))

# Calculate from Rogers: Ss=(A-I)*Sa*(A-I)T
S_tmp_2 = np.dot(np.dot((AKvmr - np.identity(AKvmr.shape[0])),a.sa[a.x_stop:a.x_stop+a.n_layer,a.x_stop:a.x_stop+a.n_layer]),(AKvmr - np.identity(AKvmr.shape[0])).T)
S_err_tmp_2 = np.sqrt(np.dot(np.dot(aprdensprf_int,S_tmp_2),aprdensprf_int.T))
#S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))

S_sys.append((S_tmp_1,S_tmp_2,S_err_tmp_1,S_err_tmp_2,'smoothing'))

# Calculate measurement error using SNR calculated from the spectrum for both primary and interfering species

# Calculate from Rogers: Sm = Dx * Sm * Dx.T

S_tmp = np.dot(np.dot(Dx,a.se),Dx.T)
S_err_tmp_1 = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp[:a.n_layer,:a.n_layer]),aprdensprf.T))
S_err_tmp_2 = np.sqrt(np.dot(np.dot(aprdensprf_int,S_tmp[a.n_layer:,a.n_layer:]),aprdensprf_int.T))

S_ran.append((S_tmp[:a.n_layer,:a.n_layer],S_tmp[a.n_layer:,a.n_layer:],S_err_tmp_1,S_err_tmp_2,'measurement'))

f.write('%s     %.4E\n'%('Column amount = ', retdenscol))
f.write('%s     %f\n'%('DOFS (total column) = ', np.trace(AK[a.x_start:a.x_stop,a.x_start:a.x_stop])))
f.write('%s     %f\n'%('Sm (%) = ', (S_ran[0][2]/retdenscol)*100))
f.write('%s     %f\n'%('Ss (%) = ', (S_sys[0][2]/retdenscol)*100))
f.write('%s     %.4E\n'%('Column amount of interferer= ', retdenscol_int))
f.write('%s     %f\n'%('DOFS of interferer (total column) = ', np.trace(AK[a.x_stop:a.x_stop+a.n_layer,a.x_stop:a.x_stop+a.n_layer])))
f.write('%s     %f\n'%('Sm_int (%) = ', (S_ran[0][3]/retdenscol_int)*100))
f.write('%s     %f\n'%('Ss_int (%) = ', (S_sys[0][3]/retdenscol_int)*100))

# Interference error 1 - retrieval parameters
AK_int1 = AK[:a.x_start,a.x_start:a.x_start+2*a.n_layer].T
Se_int1 =  a.sa[:a.x_start,:a.x_start]
  
S_tmp = np.dot(np.dot(AK_int1,Se_int1),AK_int1.T)
S_err_tmp_1 = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp[:a.n_layer,:a.n_layer]),aprdensprf.T))
S_err_tmp_2 = np.sqrt(np.dot(np.dot(aprdensprf_int,S_tmp[a.n_layer:,a.n_layer:]),aprdensprf_int.T))

f.write('%s     %f\n'%('Sint1_random (retrieval params) (%) = ', S_err_tmp_1/retdenscol*100))
f.write('%s     %f\n'%('Sint1_random_int (retrieval params) (%) = ', S_err_tmp_2/retdenscol_int*100))
S_ran.append((S_tmp[a.n_layer:,a.n_layer:],S_tmp[:a.n_layer,:a.n_layer],S_err_tmp_1,S_err_tmp_2,'retrieval parameters'))    

# Retrieval Species: Interference error 2 - interfering species 

n_int2 = a.n_profile + a.n_column - 1
n_int2_column = (a.n_profile-1)*a.n_layer + a.n_column

AK_int2 = AK[a.x_stop:a.x_stop+n_int2_column,a.x_start:a.x_stop].T

Se_int2 = a.sa[a.x_stop:a.x_stop+n_int2_column,a.x_stop:a.x_stop+n_int2_column]

S_tmp_1 = np.dot(np.dot(AK_int2,Se_int2),AK_int2.T)
S_err_tmp_1 = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp_1),aprdensprf.T))

f.write('%s     %f\n'%('Sint2_random (intf. spec.) (%) = ', (S_err_tmp_1/retdenscol)*100))


# Interfering species: Interference error 2 - interfering species 

tmp = AK
for x in range(0,a.n_layer):
  tmp = np.delete(tmp,a.x_stop,axis=0)
AK_int2 = tmp[a.x_start:a.x_start+n_int2_column,a.x_stop:a.x_stop+a.n_layer].T   
for x in range(0,a.n_layer):
  tmp = np.delete(a.sa,a.x_stop,axis=0)
for x in range(0,a.n_layer):
  tmp = np.delete(tmp,a.x_stop,axis=1)
Se_int2 = tmp[a.x_start:a.x_start+n_int2_column,a.x_start:a.x_start+n_int2_column]

S_tmp_2 = np.dot(np.dot(AK_int2,Se_int2),AK_int2.T)
S_err_tmp_2 = np.sqrt(np.dot(np.dot(aprdensprf_int,S_tmp_2),aprdensprf_int.T))

f.write('%s     %f\n'%('Sint2_random_int (intf. spec.) (%) = ', (S_err_tmp_2/retdenscol_int)*100))

S_ran.append((S_tmp_1,S_tmp_2,S_err_tmp_1,S_err_tmp_2,'interfering species'))

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
      try:
        if (s == 'zshift'):    
          Sb = np.zeros((len(zerolev_band_b),len(zerolev_band_b)))
          for i in range(0,len(zerolev_band_b)):
            Sb[i,i] = float( b.value['sb.band.'+zerolev_band_b[i]+'.zshift.'+k])**2
        elif (s == 'lineInt' or s == 'lineTAir' or s == 'linePAir'):
          Sblist = []
          Sb = np.zeros((1,1))
          Sb[0,0] = float( b.value['sb.'+s+'.'+primgas+'.'+k])**2
          Sblist.append(Sb)
          Sb[0,0] = float( b.value['sb.'+s+'.'+intgas+'.'+k])**2
          Sblist.append(Sb)
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
        if (s == 'lineInt' or s == 'lineTAir' or s == 'linePAir'):
          S_tmp = np.array(np.dot(np.dot(np.matrix(DK[:,0]).T,Sblist[0]),np.matrix(DK[:,0]))) # matrix conversion transposes the array, so transposition is reversed
          S_err_tmp_1 = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp[:a.n_layer,:a.n_layer]),aprdensprf.T))
          S_tmp = np.array(np.dot(np.dot(np.matrix(DK[:,1]).T,Sblist[1]),np.matrix(DK[:,1]))) # matrix conversion transposes the array, so transposition is reversed
          S_err_tmp_2 = np.sqrt(np.dot(np.dot(aprdensprf_int,S_tmp[a.n_layer:,a.n_layer:]),aprdensprf_int.T)) 
        else:
          S_tmp = np.dot(np.dot(DK,Sb),DK.T)
          S_err_tmp_1 = np.sqrt(np.dot(np.dot(aprdensprf,S_tmp[:a.n_layer,:a.n_layer]),aprdensprf.T))
          S_err_tmp_2 = np.sqrt(np.dot(np.dot(aprdensprf_int,S_tmp[a.n_layer:,a.n_layer:]),aprdensprf_int.T))
#        S_tmp2 = np.dot(np.dot(np.diag(aprdensprf),S_tmp),np.diag(aprdensprf))
        f.write('%s     %f\n'%('S.'+s+'.'+k+' (%) =', (S_err_tmp_1/retdenscol)*100))
        f.write('%s     %f\n'%('S.'+s+'.'+k+'_int (%) =', (S_err_tmp_2/retdenscol_int)*100))

        if (k == 'random'):      
          S_ran.append((S_tmp[:a.n_layer,:a.n_layer],S_tmp[a.n_layer:,a.n_layer:],S_err_tmp_1,S_err_tmp_2,s))

        if (k == 'systematic'):
          S_sys.append((S_tmp[:a.n_layer,:a.n_layer],S_tmp[a.n_layer:,a.n_layer:],S_err_tmp_1,S_err_tmp_2,s))

S_random_err = 0
S_random_err_int = 0
S_random = np.zeros((a.n_layer,a.n_layer))
S_random_int = np.zeros((a.n_layer,a.n_layer))

for k in range(0,len(S_ran)):
  S_random_err += S_ran[k][2]**2
  S_random_err_int += S_ran[k][3]**2
  S_random += S_ran[k][0]
  S_random_int += S_ran[k][1]
S_random_err = np.sqrt(S_random_err)
S_random_err_int = np.sqrt(S_random_err_int)

f.write('%s     %f\n'%('Total random error (%) = ', (S_random_err/retdenscol)*100))
f.write('%s     %f\n'%('Total random error of the inteferer (%) = ', (S_random_err_int/retdenscol_int)*100))

S_systematic_err = 0
S_systematic_err_int = 0
S_systematic =  np.zeros((a.n_layer,a.n_layer))
S_systematic_int =  np.zeros((a.n_layer,a.n_layer))

for k in range(1,len(S_sys)):
  S_systematic_err += S_sys[k][2]**2
  S_systematic_err_int += S_sys[k][3]**2
  S_systematic += S_sys[k][0]
  S_systematic_int += S_sys[k][1]
S_systematic_err = np.sqrt(S_systematic_err)
S_systematic_err_int = np.sqrt(S_systematic_err_int)

f.write('%s     %f\n'%('Total systematic error (%) = ', (S_systematic_err/retdenscol)*100))
f.write('%s     %f\n'%('Total systematic error of the interferer (%) = ', (S_systematic_err_int/retdenscol_int)*100))

f.write('%s     %f\n'%('Total error exclusive of smoothing error (%) = ',np.sqrt(S_systematic_err**2 + S_random_err**2)/retdenscol*100))
f.write('%s     %f\n'%('Total error exclusive of smoothing error of the interferer (%) = ',np.sqrt(S_systematic_err_int**2 + S_random_err_int**2)/retdenscol_int*100))
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
    fout.write('  %s%s %s  %.4e\n'%(S_ran[k][4][0].upper(),S_ran[k][4][1:],'random uncertainty in mol cm^-2 :', S_ran[k][2]))
  for k in range(0, len(S_sys)):
    fout.write('  %s%s %s  %.4e\n'%(S_sys[k][4][0].upper(),S_sys[k][4][1:],'systematic uncertainty in mol cm^-2 :', S_sys[k][2]))
  fout.write('\n\n')
  fout.write('  %s  %.4e\n'%('Total random uncertainty of the interferer in mol cm^-2 :', S_random_err_int))
  fout.write('  %s  %.4e\n'%('Total systematic uncertainty of the interferer in mol cm^-2 :', S_systematic_err_int))
  for k in range(0, len(S_ran)):
    fout.write('  %s%s %s  %.4e\n'%(S_ran[k][4][0].upper(),S_ran[k][4][1:],'random uncertainty of the interferer in mol cm^-2 :', S_ran[k][3]))
  for k in range(0, len(S_sys)):
    fout.write('  %s%s %s  %.4e\n'%(S_sys[k][4][0].upper(),S_sys[k][4][1:],'systematic uncertainty of the interferer in mol cm^-2 :', S_sys[k][3]))

# Output Covariance matrices in mol/cm^2 for main gas

try:
  output = b.value['out.ssystematic'].strip()
except KeyError, e:
  output = False

if (output == 'T'):
  try:
    filename = b.value['file.out.ssystematic']
  except KeyError, e:
    filename = 'ssystematic.output'
  tmp = np.dot(np.dot(np.diag(aprdensprf),S_systematic),np.diag(aprdensprf))
  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
    fout.write('          ' + str(S_systematic.shape[0]) + '          ' + str(S_systematic.shape[1]) +'\n')
    for row in tmp:
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
  tmp = np.dot(np.dot(np.diag(aprdensprf),S_random),np.diag(aprdensprf))
  with open(filename, 'w') as fout: 
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
    fout.write('          ' + str(S_random.shape[0]) + '          ' + str(S_random.shape[1])+'\n')
    for row in tmp:
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
    fout.write('          ' + str(S_ran[0][0].shape[0]) + '          ' + str(S_ran[0][0].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_ran)):
      fout.write('  '+S_ran[k][4].upper()+' ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
      tmp = np.dot(np.dot(np.diag(aprdensprf),S_ran[k][0]),np.diag(aprdensprf))
      for row in tmp:
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
      fout.write('  '+S_sys[k][4].upper()+' ERROR COVARIANCE MATRIX IN MOL CM^-2\n')
      tmp = np.dot(np.dot(np.diag(aprdensprf),S_sys[k][0]),np.diag(aprdensprf))
      for row in tmp:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

# Output Covariance matrices in mol/cm^2 for interferer

try:
  output = b.value['out.ssystematic.interferer'].strip()
except KeyError, e:
  output = False

if (output == 'T'):
  try:
    filename = b.value['file.out.ssystematic.interferer']
  except KeyError, e:
    filename = 'ssystematic.output.interferer'
  tmp = np.dot(np.dot(np.diag(aprdensprf_int),S_systematic_int),np.diag(aprdensprf_int))
  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRIX OF INTERFERER IN MOL CM^-2\n')
    fout.write('          ' + str(S_systematic_int.shape[0]) + '          ' + str(S_systematic_int.shape[1]) +'\n')
    for row in tmp:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.interferer'].strip()
except KeyError, e:
  output = False
if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.interferer']
  except KeyError, e:
    filename = 'srandom.output.interferer'
  tmp = np.dot(np.dot(np.diag(aprdensprf_int),S_random_int),np.diag(aprdensprf_int))
  with open(filename, 'w') as fout: 
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRIX OF INTERFERER IN MOL CM^-2\n')
    fout.write('          ' + str(S_random_int.shape[0]) + '          ' + str(S_random_int.shape[1])+'\n')
    for row in tmp:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.interferer.all'].strip()
except KeyError, e:
  output = False

if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.interferer.all']
  except KeyError, e:
    filename = 'srandom.interferer.all.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRICES IN MOL CM^-2\n')
    fout.write('          ' + str(S_ran[0][1].shape[0]) + '          ' + str(S_ran[0][1].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_ran)):
      fout.write('  '+S_ran[k][4].upper()+' ERROR COVARIANCE MATRIX OF INTERFERER IN MOL CM^-2\n')
      tmp = np.dot(np.dot(np.diag(aprdensprf),S_ran[k][1]),np.diag(aprdensprf))
      for row in tmp:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

try:
  output = b.value['out.ssystematic.interferer.all'].strip()
except KeyError, e:
  output = False

if (output == 'T'): 
  try:
    filename = b.value['file.out.ssytematic.interferer.all'].strip()
  except KeyError, e:
    filename = 'ssystematic.interferer.all.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRICES OF INTERFERER IN MOL CM^-2\n')
    fout.write('          ' + str(S_sys[0][1].shape[0]) + '          ' + str(S_sys[0][1].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_sys)):
      fout.write('  '+S_sys[k][4].upper()+' ERROR COVARIANCE MATRIX OF INTERFERER IN MOL CM^-2\n')
      tmp = np.dot(np.dot(np.diag(aprdensprf),S_sys[k][1]),np.diag(aprdensprf))
      for row in tmp:
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
    fout.write('          ' + str(S_systematic.shape[0]) + '          ' + str(S_systematic.shape[1]) +'\n')
    for row in  S_systematic:
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
    fout.write('          ' + str(S_random.shape[0]) + '          ' + str(S_random.shape[1])+'\n')
    for row in S_random:
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
      fout.write('  '+S_ran[k][4].upper()+' ERROR COVARIANCE MATRIX IN VMR UNITS \n')
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
      fout.write('  '+S_sys[k][4].upper()+' ERROR COVARIANCE MATRIX IN VMR UNITS \n')
      for row in  S_sys[k][0]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

# Output Covariance matrices in vmr of interferer

try:
  output = b.value['out.ssystematic.interferer.vmr'].strip()
except KeyError, e:
  output = False

if (output == 'T'):
  try:
    filename = b.value['file.out.ssystematic.interferer.vmr']
  except KeyError, e:
    filename = 'ssystematic.vmr.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRIX OF INTERFERER IN VMR UNITS\n')
    fout.write('          ' + str(S_systematic_int.shape[0]) + '          ' + str(S_systematic_int.shape[1]) +'\n')
    for row in  S_systematic_int:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.interferer.vmr'].strip()
except KeyError, e:
  output = False
if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.vmr']
  except KeyError, e:
    filename = 'srandom.vmr.interferer.output'

  with open(filename, 'w') as fout: 
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRIX OF INTERFERER IN VMR UNITS \n')
    fout.write('          ' + str(S_random_int.shape[0]) + '          ' + str(S_random_int.shape[1])+'\n')
    for row in S_random_int:
      fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))

try:
  output = b.value['out.srandom.interferer.vmr.all'].strip()
except KeyError, e:
  output = False
if (output == 'T'):

  try:
    filename = b.value['file.out.srandom.interferer.vmr.all']
  except KeyError, e:
    filename = 'srandom.interferer.vmr.all.output'
  with open(filename, 'w') as fout:
    fout.write(a.header + 'RANDOM ERROR COVARIANCE MATRICES OF INTERFERER IN VMR UNITS \n')
    fout.write('          ' + str(S_ran[0][1].shape[0]) + '          ' + str(S_ran[0][1].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_ran)):
      fout.write('  '+S_ran[k][4].upper()+' ERROR COVARIANCE MATRIX OF INTERFERER IN VMR UNITS \n')
      for row in S_ran[k][1]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')

try:
  output = b.value['out.ssystematic.interferer.vmr.all'].strip()
except KeyError, e:
  output = False

if (output == 'T'): 
  try:
    filename = b.value['file.out.ssytematic.interferer.vmr.all'].strip()
  except KeyError, e:
    filename = 'ssystematic.interferer.vmr.all.output'

  with open(filename, 'w') as fout:
    fout.write(a.header + 'SYSTEMATIC ERROR COVARIANCE MATRICES OF INTERFERER IN VMR UNITS \n')
    fout.write('          ' + str(S_sys[0][1].shape[0]) + '          ' + str(S_sys[0][1].shape[1])+'\n')
    fout.write('\n')
    for k in range(0, len(S_sys)):
      fout.write('  '+S_sys[k][4].upper()+' ERROR COVARIANCE MATRIX OF INTERFERER IN VMR UNITS \n')
      for row in  S_sys[k][1]:
        fout.write(' %s \n' % '  '.join('% .18E' % i for i in row))
      fout.write('\n')
      fout.write('\n')
