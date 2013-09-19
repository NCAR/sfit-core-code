"""
file_read.v4.py

This class reads in all files necessary for error_analysis.py

Updates from version 3:
1) Made the read in of multiple LineInt, LineTAir and LinePair arrays for the error analysis of interfering species.
2) Added input variable for number of gases for LineInt, LineTAir and LinePair so the same code can be used for 
error_analysis.py and error_analysis_interferer.py
2) Minor fixes

Updates from version 2:
1) Declare Kb and AB matrices first.  This forces Kb matrices for single column parameters (n_wav,1) instead of (n_wav,).  
This corrects issues with the vmr error covariance matrices. 

Updates from version 1:
1) Fixed issues with Kb matrix and AB matrix read in:
  a) accounted for situation when kb for FOV and Omega and not both calculated, thus no need to skip columns
  b) accounted for situation when kb for BckGrndCur and BckGrndSlp and not both calculated, thus no need to skip columns
  c) corrected number of columns for apd_fcn and phase_fcn
  d) added dwshift
  e) allowed for phase kb to be n_window x n_wav

2) Included self.n_profile, self.n_column, self.gas_profile and self.gas_column here so they can be used within and in the error calculation

Updates from version 0:
1) Reads in Kb and AB matrixes as a dictionnary to reduce number of variables and to make looping in errorcalc possible.
2) Changed num_layer to n_layer and num_window to n_window
"""

import numpy as np
import sfit4_ctl_simple as sfit4

class read_output_for_errorcalc:

  def __init__(self, ngas):
    ctl = sfit4.sfit4_ctl()
    ctl.read_ctl_file('sfit4.ctl')
    self.filemissing = False
    self.n_window = len(ctl.value['band'].strip().split())

    # Get retrieved gas list and number of retrieved gases (profiles and columns) retrieved.  
    try:
      self.gas_column = ctl.value['gas.column.list'].strip().split() 
    except KeyError, e:
      self.gas_column = []
    self.n_column = len(self.gas_column)

    try:
      self.gas_profile = ctl.value['gas.profile.list'].strip().split()
    except KeyError, e:
      self.gas_profile = []
    self.n_profile = len(self.gas_profile)
 

    # Read in Sa matrix
    nSkips = 3
    try:
      with open(ctl.value['file.out.sa_matrix'].rstrip()) as fopen:
        for i in range(0, nSkips):
          next(fopen)
        self.sa = np.array([[float(x) for x in row.split()] for row in fopen])
      Sa_file = True
    except IOError as errmsg:
#      print "Sa matrix file does not exist"
      Sa_file = False

    # Read in Ssmooth matrix
    nSkips = 2
    try:
      with open(ctl.value['file.out.ssmooth_matrix'].rstrip()) as fopen:
        for i in range(0, nSkips):
          next(fopen)
        self.Ssmooth = [[float(x) for x in row.split()] for row in fopen]
      self.Ssmooth_file = True
    except IOError as errmsg:
 #     print "Ssmooth matrix file does not exist"
      self.Ssmooth_file = False

    # Check is at least one of the above exists
    if (Sa_file == False and self.Ssmooth_file == False):
      print "Neither Sa file or Ssmooth file exists, one of which is required for error calculation, STOP"
      self.filemissing = True

    # Read in SaInv matrix - not necessary for error analysis
#    nSkips = 3
#    try:
#      with open(ctl.value['file.out.sainv_matrix'].rstrip()) as fopen:
#        for i in range(0, nSkips):
#          next(fopen)
#        self.sainv = [[float(x) for x in row.split()] for row in fopen]
#    except IOError as errmsg:
#      print "SaInv matrix file does not exist" 

    # Read in Shat matrix
#    nSkips = 3
#    try:
#      with open(ctl.value['file.out.shat_matrix'].rstrip()) as fopen:
#        for i in range(0, nSkips):
#          next(fopen)
#        self.shat = [[float(x) for x in row.split()] for row in fopen]
#    except IOError as errmsg:
#      print "Shat matrix file does not exist"


    # Create Se matrix using SNR from spectrum
    nSkips = 2
    try:
      fname = ctl.value['file.out.summary'].rstrip()
    except KeyError, e:
      fname = 'summary'
    try:
      snr = []
      npts = []
      with open(fname) as fopen:
        lines = fopen.readlines()
        for ln in range(len(lines)):
          line = lines[ln].strip().split()
          if len(line) > 0:
            if line[0] == 'IBAND':
              addline = 1
              tnpts = 0
              for j in range(0,self.n_window):
                line = lines[ln+addline].strip().split()
                npts.append(int(line[4]))  
                tnpts += int(line[4])
                addline += 1
                if len(line) < 10:
                  snr.append(float(lines[ln+addline].strip().split()[10-len(line)]))
                  addline += 1
                else:
                  snr.append(float(line[10]))
        self.se =  np.zeros((tnpts,tnpts))
        placer = 0 
        for i in range(0,self.n_window):
          for j in range(placer,placer+npts[i]):
            self.se[j,j] = snr[i]**-2
          placer += npts[i]
        del placer            

    except IOError as errmsg:
      print "Summary file does not exist and is required for error calculation, STOP"
      self.filemissing = True

    # Read in Smeas matrix
#    nSkips = 2
#    try:
#      with open(ctl.value['file.out.smeas_matrix'].rstrip()) as fopen:
#        for i in range(0, nSkips):
#          next(fopen)
#        self.Smeas = [[float(x) for x in row.split()] for row in fopen]
#      self.Smeas_file = True
#    except IOError as errmsg:
#      print "Smeas matrix file does not exist"
#      self.Smeas_file = False


#    # Read in AKx matrix
#    nSkips = 2
#    try:
#      with open(ctl.value['file.out.ak_matrix'].rstrip()) as fopen:
#        for i in range(0, nSkips):
#          next(fopen)
#        self.AKx = [[float(x) for x in row.split()] for row in fopen]
#      self.AKx_file = True
#    except IOError as errmsg:
#      print "AKx matrix file does not exist"
#      self.AKx_file = False

    # Read in K matrix
    try:
      with open(ctl.value['file.out.k_matrix'].rstrip()) as fopen:
        self.header = next(fopen)
        self.header = self.header[0:52]
        dimens = next(fopen).strip().split()
        self.n_wav = int(dimens[0])
        self.x_start = int(dimens[2])
        self.n_layer = int(dimens[3])
        self.x_stop = self.x_start + self.n_layer
        self.param_int = next(fopen).strip().split()
        self.K = [[float(x) for x in row.split()] for row in fopen]
    except IOError as errmsg:
      print "K matrix file does not exist and is required for error calculation, STOP"
      exit()
      self.filemissing = True

   # Read in Gain matrix
    nSkips = 3
    try:
      with open(ctl.value['file.out.gain_matrix'].rstrip()) as fopen:
        for i in range(0, nSkips):
          next(fopen)
        self.D = np.array([[float(x) for x in row.split()] for row in fopen])
    except IOError as errmsg:
      print "Gain matrix file does not exist and is required for error calculation, STOP"
      self.filemissing = True

    # Read in Ab matrix
    nSkips = 2
    try:
      with open(ctl.value['file.out.ab_matrix'].rstrip()) as fopen:
        for i in range(0, nSkips):
          next(fopen)
        self.param = next(fopen).split()
        Ab = np.matrix([[float(x) for x in row.split()] for row in fopen])
      read_array = [False, False, False, False, False, False, False, False, False, False,False, False,False,False]
      self.Ab = {}
      # Sort Ab matrix values

      for n in range(Ab.shape[1]):
        if self.param[n] == 'TEMPERAT' and read_array[0] == False:
          self.Ab['temperature'] = np.zeros((self.n_layer,self.n_layer))
          self.Ab['temperature'][:,:] = Ab[:,n:n+self.n_layer]
          read_array[0] = True
        if self.param[n] == 'BckGrdSlp' and read_array[1] == False:
          if (self.param[n+1] == 'BckGrdCur'):
            self.Ab['slope'] = np.zeros((self.n_layer,self.n_window))
            for i in range(0,self.n_window):
              self.Ab['slope'][:,i] = Ab[:,n+i*2]
          else:
            self.Ab['slope'][:,:] = Ab[:,n:n+self.n_window]
          read_array[1] = True
        if self.param[n] == 'BckGrdCur' and read_array[2] == False:
          if (self.param[n+1] == 'BckGrdSlp'):
            self.Ab['curvature'] = np.zeros((self.n_layer,self.n_window))
            for i in range(0,self.n_window):
              self.Ab['curvature'][:,i] = Ab[:,n+i*2]
          else:
            self.Ab['curvature'][:,:] = Ab[:,n:n+self.n_window]
          read_array[2] = True
        if self.param[n] == 'SolLnShft':
          self.Ab['solshft'] = np.zeros((self.n_layer,1)) 
          self.Ab['solshft'][:,:] = Ab[:,n]
        if self.param[n] == 'SolLnStrn':
          self.Ab['solstrnth'] = np.zeros((self.n_layer,1))
          self.Ab['solstrnth'][:,:] = Ab[:,n]
        if self.param[n] == 'SPhsErr' and read_array[10] == False:
          self.Ab['phase'] = np.zeros((self.n_layer,self.n_window))
          for i in range(0,self.n_window):
            self.Ab['phase'][:,:] = Ab[:,n:n+self.n.window]
            read_array[10] = True
        if self.param[n] == 'IWNumShft' and read_array[8] == False:
          self.Ab['wshift'] = np.zeros((self.n_layer,self.n_window))
          self.Ab['wshift'][:,:] = Ab[:,n:n+self.n_window]
        if self.param[n] == 'DWNumShft' and read_array[9] == False:
          self.Ab['dwshift'] = np.zeros((self.n_layer,self.n_profile+self.n_column-1))
          self.Ab['dwshift'][:,:] = Ab[:,n:n+self.n_profile+self.n_column-1]
        if self.param[n] == 'EmpApdFcn' and read_array[5] == False:
          tmp = ctl.value['fw.apod_fcn.order'].rstrip()      
          self.Ab['apod_fcn'] = np.zeros((self.n_layer,int(tmp)))    
          self.Ab['apod_fcn'][:,:] = Ab[:,n:n+int(tmp)]
          read_array[5] = True
        if self.param[n] == 'EmpPhsFnc' and read_array[6] == False:
          tmp = ctl.value['fw.phase_fcn.order'].rstrip()
          self.Ab['phase_fcn'] = np.zeros((self.n_layer,int(tmp))) 
          self.Ab['phase_fcn'][:,:] = Ab[:,n:n+int(tmp)]
          read_array[6] = True
        if self.param[n] == 'SZA':
          self.Ab['sza'] = np.zeros((self.n_layer,1))
          self.Ab['sza'][:,:] = Ab[:,n] 
        if self.param[n] == 'FOV' and read_array[3] == False:
          self.Ab['omega'] = np.zeros((self.n_layer,self.n_window))
          if (len(self.param) > n+1 and self.param[n+1] == 'OPD'):
            for i in range(0,self.n_window):
              self.Ab['omega'][:,i] = Ab[:,n+i*2]
          else:
            self.Ab['omega'][:,:] = Ab[:,n:n+self.n_window]
          read_array[3] = True
        if self.param[n] == 'OPD' and read_array[4] == False:
          self.Ab['max_opd'] = np.zeros((self.n_layer,self.n_window))
          if (len(self.param) > n+1 and self.param[n+1] == 'FOV'):
            for i in range(0,self.n_window):
              self.Ab['max_opd'][:,i] = Ab[:,n+i*2]
          else:
            self.Ab['max_opd'][:,:] = Ab[:,n:n+self.n_window]
          read_array[4] = True
        if self.param[n] == 'LineInt' and read_array[11] == False:
          self.Ab['lineInt'] = np.zeros((self.n_wav,ngas))
          self.Ab['lineInt'][:,:] = Ab[:,n+ngas]
          read_array[11] == True
        if self.param[n] == 'LineTAir' and read_array[12] == False:
          self.Ab['lineTAir'] = np.zeros((self.n_wav,ngas))
          self.Ab['lineTAir'][:,:] = Ab[:,n:n+ngas]
          read_array[12] == True
        if self.param[n] == 'LinePAir' and read_array[13] == False:
          self.Ab['linePAir'] = np.zeros((self.n_wav,ngas))
          self.Ab['linePAir'][:,:] = Ab[:,n:n+ngas]
          read_array[13] == True
        if self.param[n] == 'ZeroLev' and read_array[7] == False:
          self.Ab['zshift'] = np.zeros((self.n_layer,self.n_window))
          self.Ab['zshift'][:,:] = Ab[:,n:n+self.n_window]
          read_array[7] = True

    except IOError as errmsg:
#      print "Ab matrix file does not exist and is required for error calculation, STOP"
#      self.filemissing = True
      pass

    # Read in Kb matrix
    nSkips = 2
    try:
      with open(ctl.value['file.out.kb_matrix'].rstrip()) as fopen:
        for i in range(0, nSkips):
          next(fopen)
        self.param = next(fopen).split()
        Kb = np.matrix([[float(x) for x in row.split()] for row in fopen])
      read_array = [False, False, False, False, False, False, False, False, False, False,False,False,False,False]
      self.Kb = {}

      # Sort Kb matrix values

      for n in range(Kb.shape[1]):
        if self.param[n] == 'TEMPERAT' and read_array[0] == False:
          self.Kb['temperature'] = np.zeros((self.n_wav,self.n_layer))
          self.Kb['temperature'][:,:] = Kb[:,n:n+self.n_layer]
          read_array[0] = True
        if self.param[n] == 'BckGrdSlp' and read_array[1] == False:
          if (self.param[n+1] == 'BckGrdCur'):
            self.Kb['slope'] = np.zeros((self.n_wav,self.n_window))
            for i in range(0,self.n_window):
              self.Kb['slope'][:,i] = Kb[:,n+i*2]
          else:
            self.Kb['slope'][:,:] = Kb[:,n:n+self.n_window]
          read_array[1] = True
        if self.param[n] == 'BckGrdCur' and read_array[2] == False:
          if (self.param[n+1] == 'BckGrdSlp'):
            self.Kb['curvature'] = np.zeros((self.n_wav,self.n_window))
            for i in range(0,self.n_window):
              self.Kb['curvature'][:,i] = Kb[:,n+i*2]
          else:
            self.Kb['curvature'][:,:] = Kb[:,n:n+self.n_window]
          read_array[2] = True
        if self.param[n] == 'SolLnShft':
          self.Kb['solshft'] = np.zeros((self.n_wav,1)) 
          self.Kb['solshft'][:,:] = Kb[:,n]
        if self.param[n] == 'SolLnStrn':
          self.Kb['solstrnth'] = np.zeros((self.n_wav,1))
          self.Kb['solstrnth'][:,:] = Kb[:,n]
        if self.param[n] == 'SPhsErr' and read_array[10] == False:
          self.Kb['phase'] = np.zeros((self.n_wav,self.n_window))
          self.Kb['phase'][:,:] = Kb[:,n+self.n_window] 
          read_array[10] = True
        if self.param[n] == 'IWNumShft' and read_array[8] == False:
          self.Kb['wshift'] = np.zeros((self.n_wav,self.n_window))
          self.Kb['wshift'][:,:] = Kb[:,n:n+self.n_window]
        if self.param[n] == 'DWNumShft' and read_array[9] == False:
          self.Kb['dwshift'] = np.zeros((self.n_wav,self.n_profile+self.n_column-1))
          self.Kb['dwshift'][:,:] = Kb[:,n:n+self.n_profile+self.n_column-1]
        if self.param[n] == 'EmpApdFcn' and read_array[5] == False:
          tmp = ctl.value['fw.apod_fcn.order'].rstrip()      
          self.Kb['apod_fcn'] = np.zeros((self.n_wav,int(tmp)))    
          self.Kb['apod_fcn'][:,:] = Kb[:,n:n+int(tmp)]
          read_array[5] = True
        if self.param[n] == 'EmpPhsFnc' and read_array[6] == False:
          tmp = ctl.value['fw.phase_fcn.order'].rstrip()
          self.Kb['phase_fcn'] = np.zeros((self.n_wav,int(tmp))) 
          self.Kb['phase_fcn'][:,:] = Kb[:,n:n+int(tmp)]
          read_array[6] = True
        if self.param[n] == 'SZA':
          self.Kb['sza'] = np.zeros((self.n_wav,1))
          self.Kb['sza'][:,:] = Kb[:,n] 
        if self.param[n] == 'FOV' and read_array[3] == False:
          self.Kb['omega'] = np.zeros((self.n_wav,self.n_window))
          if (len(self.param) > n+1 and self.param[n+1] == 'OPD'):
            for i in range(0,self.n_window):
              self.Kb['omega'][:,i] = Kb[:,n+i*2]
          else:
            self.Kb['omega'][:,:] = Kb[:,n:n+self.n_window]
          read_array[3] = True
        if self.param[n] == 'OPD' and read_array[4] == False:
          self.Kb['max_opd'] = np.zeros((self.n_wav,self.n_window))
          if (len(self.param) > n+1 and self.param[n+1] == 'FOV'):
            for i in range(0,self.n_window):
              self.Kb['max_opd'][:,i] = Kb[:,n+i*2]
          else:
            self.Kb['max_opd'][:,:] = Kb[:,n:n+self.n_window]
          read_array[4] = True
        if self.param[n] == 'LineInt' and read_array[11] == False:
          self.Kb['lineInt'] = np.zeros((self.n_wav,ngas))
          self.Kb['lineInt'][:,:] = Kb[:,n+ngas]
          read_array[11] == True
        if self.param[n] == 'LineTAir' and read_array[12] == False:
          self.Kb['lineTAir'] = np.zeros((self.n_wav,ngas))
          self.Kb['lineTAir'][:,:] = Kb[:,n:n+ngas]
          read_array[12] == True
        if self.param[n] == 'LinePAir' and read_array[13] == False:
          self.Kb['linePAir'] = np.zeros((self.n_wav,ngas))
          self.Kb['linePAir'][:,:] = Kb[:,n:n+ngas]
          read_array[13] == True
        if self.param[n] == 'ZeroLev' and read_array[7] == False:
          self.Kb['zshift'] = np.zeros((self.n_wav,self.n_window))
          self.Kb['zshift'][:,:] = Kb[:,n:n+self.n_window]
          read_array[7] = True

    except IOError as errmsg:
      print "Kb matrix file does not exist and is required for error calculation, STOP"
      self.filemissing = True
