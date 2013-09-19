import sys, os 
import numpy as np
import sfit4_ctl as sfit4

gas = 'CH4_SNR300'

# TAO
lat = 43.6 # depends on site
long = 79.36 # depends on site
alt = 174 # depends on site
rearth = 6367.44 #find radius at eureka - site specific
n_layer = 48 # find a way not to have to define this

vmrfile = '../vmr/01-43_WACCM6_NO.ref' #

# folder conatain the bnr folder 
base = '/net/glacier/data/1/fts/tao-da8/'

# folder for output folder
base2 = '/net/glacier/data/1/sconway/tao-da8-sfit4/'

def PrepareRef(zptfile,vmrfile):
    """
    Prepare input file for the retrieval by copying the
    content of file zpt file and the vmr file to the reference.prf file
    """
    # Read the original files

    with open(zptfile) as fin:
      lines = fin.readlines()
    with open(vmrfile) as fin:
      lines2 = fin.readlines()

    # Open the output file
    foutname = 'reference.prf'
    with open(foutname, 'w') as fout:
    # Loop over lines and write them to file
      for ln in range(len(lines)):
         line = lines[ln]
         if (ln == 0):
              #strip the leading and trailing spaces
              lst = line.rstrip().lstrip().split(' ')

              #get rid of any extra zeros in the list
              lst = [el for el in lst if el != '']

              #make new line
              line = '   %s   %s   %1.0f\n' % (lst[0], lst[1], 99)
         fout.write(line)
#      fout.write(' \n')
      for ln in range(len(lines2)):
          line = lines2[ln]
          fout.write(line)

def PrepareSpec(bnrfile):
    """
    Prepare input file for the prepspec by inputing the lat, long, 
    altitude, radius of the earth, zero fill factor and bnr file name 
    and copying the bandpass, resolution from the mdb file 
    NB ratio flag set to 0 as default so no ratio file name 
    """
    outfile = 'pspec.input'
    zerofill = 2 # what is this?
    ratioflag = 0 # what is this?
    
    with open(outfile, 'w') as fout:
      fout.write('# prepspec.input\n')
      fout.write('# input file for prepspec.f90\n')
      fout.write('# prepspec loops through this file and creates an ascii \n')
      fout.write('# file with a block for each spectrum that will be fit\n')
      fout.write('# Latitude of Observation [+N, 90 - -90]\n')
      fout.write(str(lat)+'\n')
      fout.write('# Longitude of Observation[+E, 0 - 360]\n')
      fout.write(str(long)+'\n')
      fout.write('# Altitude of Observation [masl]\n')
      fout.write(str(alt)+'\n')
      fout.write('#\n')
      fout.write('# filter bands and regions for calculating SNR\n')
      fout.write('6\n')
      fout.write('f1  4038.727 4038.871\n')
      fout.write('f2  3381.155 3381.536\n')
      fout.write('f3  2924.866 2925.100\n')
#      fout.write('f4  2526.228 2526.618\n')
      fout.write('f5  1985.260 1985.510\n')
      fout.write('f6  1139.075 1139.168\n')
      fout.write('f8  907.854  907.977\n')
      fout.write('#\n')
      fout.write('# number of data blocks in the output ascii file\n')
      fout.write('# = [# binary formatted spectra] X [# fit regions] (often but not necessarily)\n')
      fout.write(str(1)+'\n')
      fout.write('# each block contains:\n')
      fout.write('# bnr spectra file name\n')
      fout.write('# radius of earth, zero fill factor, ratioflag\n')  
      fout.write('# ratio file name (bnr format) if ratioflag eq 1, skip if 0\n')
      fout.write(bnrfile+'\n')
      line = '%1.2f    %1.0f    %0.1d\n' % (rearth, zerofill, ratioflag)
      fout.write(line)
    os.system('/net/glacier/data/1/sconway/sfit4/pspec')
    with open('t15asc.4', 'r') as f:
      lines = f.readlines()
    fov = float(lines[2][51:56])
    ctl = sfit4.sfit4_ctl()
    ctl.read_ctl_file('sfit4.ctl')
    nband = len(ctl.value['band'].rstrip().lstrip().split())
    for i in range(1,nband):
       ctl.replace_in_file('sfit4.ctl','band.'+str(i)+'.omega',str(fov))

def Cleanup():
    """ 
    Moves files to a new folder to avoid files being overwritten.
    """
    odir = base2 + gas + '/' +  bfile[0:12] + '/'
    if (os.path.isdir(base2 + gas) is False):
        os.system('mkdir ' + base2 + gas)
    if (os.path.isdir(odir) is False):
        os.system('mkdir ' + odir)
    os.system('cp * ' + odir)
    os.system('rm ' + odir + '*.pyc')
    os.system('rm ' + odir + '*.hbin')
    os.system('rm ' + odir + '*.sh') 
    os.system('rm ' + odir + '*.sfit4.ctl')    
  
#if the script is executed as if
if __name__ == '__main__':

    # Analysis
    bfile = sys.argv[1]
    filter = bfile[6:8] 
    zptfile = base + 'zpt/43_layer/' + bfile[0:6]+'.zpt'  
    mdbfile = base + 'mdb/49_layer/' + filter + '/'+ bfile[0:-4]+'.mdb'
    bnrfile = base + 'bnr/' + filter + '/' + bfile

    PrepareRef(zptfile,vmrfile)

    PrepareSpec(bnrfile)

    os.system('/net/glacier/data/1/sconway/sfit4/sfit4')

    Cleanup()
