import numpy as np
import sfit4_ctl_simple as sfit4


class prfs:

  def __init__(self, gas):
    # Read in profiles matrix
    ctl = sfit4.sfit4_ctl()
    ctl.read_ctl_file('sfit4.ctl')
    nSkips = 3
    try:
      fname = ctl.value['file.out.retprofiles'].strip()
    except KeyError, e:
      fname = 'rprfs.table'
 
    with open(fname) as fopen:
      for i in range(0, nSkips):
        next(fopen)
      column_names = next(fopen).strip().split()
      full_matrix = np.array([[float(x) for x in row.split()] for row in fopen])
#    self.z = full_matrix[:,0]
    self.zbar = full_matrix[:,1]
    self.T = full_matrix[:,2]
#    self.P = full_matrix[:,3]
    self.airmass = full_matrix[:,4]
    
    for i in range(0,len(column_names)):
      if column_names[i] == gas:
        self.rprf = full_matrix[:,i]
        column_num = i
    try:
      fname = ctl.value['file.out.aprprofiles'].strip()
    except KeyError, e:
      fname = 'aprfs.table'

    with open(fname) as fopen:
      for i in range(0, nSkips+1):
        next(fopen)
      full_matrix = np.array([[float(x) for x in row.split()] for row in fopen])
      self.aprf = full_matrix[:,column_num]



      
    
