#-----------------------------------------------------------------
# Name:
#      sb.ctl
#
# Purpose:
#      This is the ctl file for the error analysis module located 
#      in Layer1Mods.py
#
# Notes:
#       1) Comments are denoted by '#'
#       2) Sb are specified in either native or fractional units. See below
#        
#
# Units:
#    -- Sb for Temperature, SZA, and FOV can either be specified in native 
#       or fractional units. This is controlled by the units flag.
#    -- Following list shows units for various Sb's:
#
#     Parameter Name          Description                           Units
#     ---------------         --------------------------            ----------------------------------
#     temperature             Temperature                           Native [Kelvin] or fractional                      
#     solshft                 Solar line shift                      Native [cm^-1]
#     solstrnth               Solar line strength                   Fractional       
#     phase                   Phase                                 Native [Radians]
#     wshift                  Wavelength shift                      Fractional 
#     dwshift                 Differential Wavelength shift         ******Currently NOT implemented***** 
#     sza                     Solar zenith angle                    Native [degrees] or fractional ??
#     lineInt                 Line intensity                        Fractional
#     lineTAir                Line temperature broadening           Fractional
#     linePAir                Line pressure broadening              Fractional
#     slope                   Background slope                      Native [cm^-1]
#     curvature               Background curvature                  Native [cm^-2]
#     apod_fcn                Empirical apodization Function        Fractional
#     phase_fcn               Empirical phase function              Fractional
#     omega                   Field of view                         Native [milliradians] or fractional
#     max_opd                 Optical path difference               Fractional
#     zshift                  Zero level                            Native [0-1]
#     profile.gas             VMR of retrieval gas                  Fractional
#
# Notes:
#  1) phase and phase_fcn are different ways to describe the same parameter. 
#     It is not recommended to calculate an error on both simultaneously. 
#  2) dwshift is currently not functional. Do NOT use.
#     
#-----------------------------------------------------------------

                        #-------#
                        # Flags #
                        #-------#
#-------------
# Output flags
#-------------
VMRoutFlg                  = T      # T = output error covariance matrices in VMR
MolsoutFlg                 = T      # T = output error covariance matrices in molecules cm^-2
out.total                  = T      # T = write out total random error covariance matrix
out.srandom                = T      # T = write out random error covariance matrix
out.ssystematic            = T      # T = write out systematic error covariance matrix 

#------------
# Input Flags
#------------
SeInputFlg = F                # This flag determines where the Se matrix is read in.
                              # If = T, the Se matrix is read in from sfit output file: file.out.seinv_vector
                              # This method takes into account de-weighting of the SNR set in the sfit4.ctl file
                              # If = F, the Se is taken from the summary file. These are the actual SNR values
                              # taken from the t15asc files.

#-----------------------------------------------
# Units flag indicate whether the Sb is given in
# native units or scaled
#     F = Native Units, T = Fractional
#-----------------------------------------------
sb.temperature.random.scaled                = F         # If = T (fractional) --> scaled by a priori
sb.temperature.systematic.scaled            = F         # If = T (fractional) --> scaled by a priori
sb.sza.random.scaled                        = T
sb.sza.systematic.scaled                    = T
sb.omega.random.scaled                      = T
sb.omega.systematic.scaled                  = T     # FOV Change!!


                    #-------------------#
                    # Output file names #
                    #-------------------#
file.out.total                        = Stotal.output
file.out.total.vmr                    = Stotal.vmr.output
file.out.srandom                      = Srandom.output
file.out.srandom.vmr                  = Srandom.vmr.output
file.out.ssystematic                  = Ssystematic.output
file.out.ssystematic.vmr              = Ssystematic.vmr.output
file.out.error.summary                = Errorsummary.output
file.out.avk                          = avk.output


                    #-------------------#
                    #     Sb values     #
                    #-------------------#
#-----------------------------------------------
# Sb for temperature & profile.gas:
# -- Specify diagonals of Sb matrix, i.e.
#    one value for each layer (Descending, first
#    value is top layer)  
#------------------------------------------------
sb.temperature.random       =   
9  9  9  9  9  9  9  7  7  7  7  6  6  5  5  2  2  2  2  2  2  2  2  2
2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

sb.temperature.systematic   =
9  9  9  9  9  9  9  7  7  7  7  6  6  5  5  2  2  2  2  2  2  2  2  2
2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

sb.profile.H2O.random       = 
9  9  9  9  9  9  9  7  7  7  7  6  6  5  5  2  2  2  2  2  2  2  2  2
2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

sb.profile.H2O.systematic   = 
9  9  9  9  9  9  9  7  7  7  7  6  6  5  5  2  2  2  2  2  2  2  2  2
2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2

#----------------------------------------------------------------------------
# Micro-window dependent Sb:
# -- Number of entries corresponds to the number of bands
#    *** The order of the entries for sb must correspond to the order of
#        bands specified in "band = " in the sfit4.ctl file ***
#----------------------------------------------------------------------------
sb.omega.random         = 0.001 
sb.omega.systematic     = 0.001 

sb.sza.random           = 0.05  
sb.sza.systematic       = 0.01  

sb.phase.random         = 0.001 
sb.phase.systematic     = 0.001 

sb.wshift.random        = 0.001 
sb.wshift.systematic    = 0.001 

sb.slope.random         = 0.001 
sb.slope.systematic     = 0.001 

sb.curvature.random     = 0.001 
sb.curvature.systematic = 0.001 

sb.max_opd.random       = 0.00001 
sb.max_opd.systematic   = 0.00001 

#-----------------------------------------------------------
# Sb for zshift is micro window dependent. However, Sb's are 
# only to be specified in microwindows where zshift is not
# retrieved. For example, if you have two microwindows and 
# you retrieve the first (1), then you would specify:
#  sb.band.2.zshift.random and sb.band.2.zshift.systematic
# The number corresponds to the band number in the sfit4.ctl
# file.
#-----------------------------------------------------------
sb.band.1.zshift.random      = 0.01
sb.band.1.zshift.systematic  = 0.01

#----------------
# Single value Sb
#----------------
sb.solshft.random       = 0.005
sb.solshft.systematic   = 0.005
sb.solstrnth.random     = 0.001
sb.solstrnth.systematic = 0.001

#------------------------------------------------------------
# If the apodization and phase function are used in the forward model,
# (fw.apod_fcn.type and fw.phase_fcn.type = 2 or 3 ONLY) then the 
# Sb's are given for each term of the polynomial or fourier series.
# Example: 2nd order polynomial would require 2 Sbs.
# If the apodization or phase function is not included in the forward
# model then the default value for calculating Kb is a 3rd order
# polynomial. Therefore 3 values of Sb must be given. 
#------------------------------------------------------------
sb.apod_fcn.random      = 0.05 0.05 0.05 0.05
sb.apod_fcn.systematic  = 0.05 0.05 0.05 0.05
sb.phase_fcn.random     = 0.05 0.05 0.05 0.05  
sb.phase_fcn.systematic = 0.05 0.05 0.05 0.05

#----------------------------------------------------------
# Sb's for lineInt, lineTair, and linePair are specific to
# an individual gas. They should be specified as:
#      lineInt_<GAS>    example: lineInt_H2O
#      linePAir_<GAS>   example: linePAir_H2O
#      lineTAir_<GAS>   example: lineTAir_H2O
#
# Sb's (lineInt, linePAir, lineTAir) should be specified 
# for gases given in kb.line.gas
#----------------------------------------------------------
sb.lineInt_HNO3.random                   = 0.0
sb.lineTAir_HNO3.random                  = 0.0
sb.linePAir_HNO3.random                  = 0.0
sb.lineInt_HNO3.systematic               = 0.1     
sb.lineTAir_HNO3.systematic              = 0.1 
sb.linePAir_HNO3.systematic              = 0.1 

#--------------------------------------------------------------
# DO NOT USE!!!!!
# Kbs are calculated for dwshift for all interfering species.
# (i.e. all gases, except the first gas listed (the primary gas)
#  They should be specified as:
#      dwshift_<GAS>    example: dwshift_H2O
#-------------------------------------------------------------- 
#sb.dwshift_H2O.random               = 0.1 
#sb.dwshift_H2O.systematic           = 0.1


#---------------------------------------------------------
# These flags indicate which errors are included in the 
# total random and systematic error budget. If the flag is
# set as F or if it is missing than it is NOT included in
# total error
#---------------------------------------------------------
sb.total.lineInt_HNO3                 = T
sb.total.lineTAir_HNO3                = T
sb.total.linePAir_HNO3                = T
sb.total.temperature                 = T
sb.total.sza                         = T
sb.total.measurement                 = T
sb.total.smoothing                   = F
sb.total.interfering_species         = T
sb.total.retrieval_parameters        = T







