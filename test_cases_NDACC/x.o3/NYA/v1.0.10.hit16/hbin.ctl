# hbin.input for Testing ATM 20200512
#
# Update SFIT4 v1.0
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
#file.in.linelist = /bira-iasb/projects/FTIR/tools/programs/sfit4/linelists/linelist-core/
#file.in.linelist = /Users/jamesw/FDP/sfit/400/linelist-core/
#file.in.linelist = /data/ebaumer/Code/linelist-core/
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
hitran.nr = 99
#
#
hitran.files =
001_H2O/01_hit09.par
002_CO2/02_hit08_f53.par
003_O3/03_hit08.par
004_N2O/04_hit08.par
005_CO/05_hit08.par
006_CH4/06_hit08.par
007_O2/07_hit08.par
008_NO/08_hit08.par
009_SO2/09_hit08.par
010_NO2/10_hit08.par
011_NH3/11_hit08.par
012_HNO3/12_hit08.par
013_OH/13_hit08.par
014_HF/14_hit08.par
015_HCL/15_hit08.par
016_HBR/16_hit08.par
017_HI/17_hit08.par
018_CLO/18_hit08.par
019_OCS/19_hit08.par
020_H2CO/20_hit08.par
021_HOCL/21_hit08.par
022_HO2/33_hit08.par
023_H2O2/25_hit08.par
024_HONO/
025_HO2NO2/
026_N2O5/2007.sudo.n2o5
027_CLONO2/2007.sudo.clono2
028_HCN/23_hit08.par
029_CH3F/
030_CH3CL/24_hit08_f53.par
031_CF4/2007.sudo.cf4
032_CCL2F2/2007.sudo.ccl2f2
033_CCL3F/2007.sudo.ccl3f
034_CH3CCL3/
035_CCL4/2007.sudo.ccl4
036_COF2/29_hit08.par
037_COCLF/2007.sudo.coclf
038_C2H6/c2h6_2720_3100.101.27
039_C2H4/38_hit08.par
040_C2H2/26_hit08.par
041_N2/22_hit08.par
042_CHF2CL/2007.sudo.chf2cl
043_COCL2/2007.sudo.cocl2
044_CH3BR/40_hit08.par
045_CH3I/
046_HCOOH/32_hit08.par
047_H2S/31_hit08.par
048_CHCL2F/
049_O2CIA/o2cia_20060420.101
050_SF6/2007.sudo.sf6
051_NF3/
052_N2CIA/
053_OTHER/
054_OTHER/
055_PH3
056_OTHER/
057_OTHER/
058_OCLO/
059_F134A/
060_C3H8/2005ppn.sudo.c3h8
061_F142B/2007.sudo.f142b
062_CFC113/2007.sudo.f113
063_F141B/
064_CH3OH/39_hit08.par
065_OTHER
066_OTHER
067_PAN/2007.sudo.pan
068_CH3CHO/2007.sudo.ch3cho
069_CH3CN/41_hit08.par
070_OTHER/
071_CH3COOH/
072_C5H8/
073_MVK/
074_MACR/
075_C3H6/
076_C4H8/
077_OTHER/
078_OTHER/
079_OTHER/
080_OTHER/
081_OTHER/
082_OTHER/
083_OTHER/
084_OTHER/
085_OTHER/
086_OTHER/
087_OTHER/
088_OTHER/
089_OTHER/
090_OTHER/
091_OTHER/
092_OTHER/
093_OTHER/
094_OTHER/
095_OTHER/
096_OTHER/
097_OTHER/
098_OTHER/
099_OTHER/
#
#
# These are all set off (*nr=0) until you verify compatibility of lines.
#
# Galatry parameters
# molecule id numbers in these files have to match the sfit molecule id
aux = gal sdv lm
aux.gal.nr = 2
aux.gal.files =
014_HF/14_hit16_Galatry.txt
015_HCL/15_hit16_Galatry.txt
#
# Speed Dependent Voigt parameter files
aux.sdv.nr = 1
aux.sdv.files =
005_CO/05_hit16_SDV.txt
#
#
# CO2 Line mixing parameters (e.g. CCl4 retrieval)
aux.lm.nr = 1
aux.lm.files =
002_CO2/002_CO2.hit16_LM1ST.par
#
#
