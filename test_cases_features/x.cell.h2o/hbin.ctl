# hbin.input for Testing HITRAN 2016
#
# implemented in version b3.99:
# the linelist directory structure is the key to the gas names and the molecule id numbers
# those id's and names must be the same in the reference.prf file
# eg a files containing hitran lines is read from one subdir in linelist then the molid will be changed
# to the 2digit integer 0NN of the subdir name and assumed to be for gas 'abcdef' from subdir 0NN_abcdef
#
# Save an ascii (HITRAN format) list file (True / False)
file.out.ascii = T
#
# Path to the directory tree where the gas subdirectories are
#
file.in.linelist = /home/mathias/linelist-core/
#
# Then the next lines are paths to each gas file that will be searched for lines in the
# desired wavenember region.  The id numbers are in sfit order - not HITRAN, BUT KEEP
# these directory names the files are in HITRAN format.  File names can be anything.
#
# Tip: put a blank after the '/' to skip that gas or just remove it (and decrease hitran.nr)
#
# 049_O2CIA/o2cia_20060420.101 - this file is a special case, molid & isotope
#  ids are dealt with differently
#
# Number of path/files to look for
hitran.nr = 1
#
#
hitran.files =
001_H2O/001_H2O.hit16.20181107
#
#
# Galatry parameters
# molecule id numbers in these files have to match the sfit molecule id
aux = gal sdv lm
# gal sdv lm
aux.gal.nr = 0
aux.gal.files =
#
# Speed Dependent Voigt parameter files --- not yet!
aux.sdv.nr = 0
aux.sdv.files =
#
#
# CO2 Line mixing parameters for Boone implementation
aux.lm.nr = 0
aux.lm.files =
#
#
