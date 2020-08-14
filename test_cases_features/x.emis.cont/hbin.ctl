# hbin.input for Testing HITRAN 2016 (d/l 20181107)
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
#file.in.linelist = /home/mathias/linelist-core/
#file.in.linelist = /bira-iasb/projects/FTIR/tools/programs/sfit4/linelists/linelist-core/
#file.in.linelist = /Users/jamesw/FDP/sfit/400/linelist-core/
file.in.linelist = ~happy_place/sfit/linelist-core/
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
001_H2O/001_H2O.hit16.20181107
002_CO2/002_CO2.hit16.20181107
003_O3/003_O3.hit16.20181107
004_N2O/004_N2O.hit16.20181107
005_CO/005_CO.hit16.20181107
006_CH4/006_CH4.hit16.20181107
007_O2/007_O2.hit16.20181107
008_NO/008_NO.hit16.20181107
009_SO2/009_SO2.hit16.20181107
010_NO2/010_NO2.hit16.20181107
011_NH3/011_NH3.hit16.20181107
012_HNO3/012_HNO3.hit16.20181107
013_OH/013_OH.hit16.20181107
014_HF/014_HF.hit16.20181107
015_HCL/015_HCl.hit16.20181107
016_HBR/016_HBr.hit16.20181107
017_HI/017_HI.hit16.20181107
018_CLO/018_ClO.hit16.20181107
019_OCS/019_OCS.hit16.20181107
020_H2CO/020_H2CO.hit16.20181107
021_HOCL/021_HOCl.hit16.20181107
022_HO2/033_HO2.hit16.20181107
023_H2O2/025_H2O2.hit16.20181107
024_HONO/
025_HO2NO2/
026_N2O5/2007.sudo.n2o5
027_CLONO2/2007.sudo.clono2
028_HCN/023_HCN.hit16.20181107
029_CH3F/
030_CH3CL/024_CH3Cl.hit16.20181107
031_CF4/2007.sudo.cf4
032_CCL2F2/2007.sudo.ccl2f2
033_CCL3F/2007.sudo.ccl3f
034_CH3CCL3/
035_CCL4/2007.sudo.ccl4
036_COF2/029_COF2.hit16.20181107
037_COCLF/2007.sudo.coclf
038_C2H6/c2h6_2720_3100.101
039_C2H4/038_C2H4.hit16.20181107
040_C2H2/026_C2H2.hit16.20181107
041_N2/022_N2.hit16.20181107
042_CHF2CL/2007.sudo.chf2cl
043_COCL2/049_COCl2.hit16.20181107
044_CH3BR/040_CH3Br.hit16.20181107
045_CH3I/
046_HCOOH/032_HCOOH.hit16.20181107
047_H2S/031_H2S.hit16.20181107
048_CHCL2F/
049_O2CIA/o2cia_20060420.101
050_SF6/2007.sudo.sf6
051_NF3/nf3.101.511.txt
052_N2CIA/n2cia.20060420.101
053_OTHER/
054_OTHER/
055_PH3/055_PH3.atm.20200512
056_OTHER/
057_OTHER/
058_OCLO/
059_F134A/
060_C3H8/2005ppn.sudo.c3h8
061_F142B/2007.sudo.f142b
062_CFC113/2007.sudo.f113
063_F141B/
064_CH3OH/039_CH3OH.hit16.20181107
065_OTHER/
066_OTHER/
067_PAN/2007.sudo.pan
068_CH3CHO/2007.sudo.ch3cho
069_CH3CN/041_CH3CN.hit16.20181107
070_OTHER/
071_CH3COOH/ch3cooh_1100_pll.txt
072_C5H8/c5h8_885_pll.txt
073_MVK/mvk_910_pll.txt
074_MACR/macr_890_pll.txt
075_C3H6/c3h6_850_pll.txt
076_C4H8/c4h8_750_pll.txt
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
aux.gal.nr = 0
aux.gal.files =
014_HF/14_hit16_Galatry.txt
015_HCL/15_hit16_Galatry.txt
#
# Speed Dependent Voigt parameter files
aux.sdv.nr = 0
aux.sdv.files =
005_CO/05_hit16_SDV.txt
#
#
# CO2 Line mixing parameters (e.g. CCl4 retrieval)
aux.lm.nr = 0
aux.lm.files =
002_CO2/002_CO2.hit16_LM1ST.par
#
