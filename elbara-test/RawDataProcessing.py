#!c:/Python25/python.exe
#-*-coding: utf-8-*-
"""
Program to process ELBARAII raw data to calibrated brightness temperatures
including RFI analysis in the time and frequency domain, as well as checks
for the plausibility of the derived brightness temperatures and for problems
with the thermal stabilization of the instrument during the measurements.

usage: RawDataProcessing <pf> [p2s]
        pf:  name of parameter file (compulsory)
        p2s: print detailed output to standard out (optional)
              1 = yes (default)
              0 = no

author: Ingo Voelksch, 21.10.2011
        Andreas Wiesmann, 20160524 changed to independent calibration per channel
"""

from libdataproc import *
import os

print('*** Processing ELBARAII raw data to calibrated brightness temperatures ***')
print('*** author: Ingo Voelksch, Andreas Wiesmann 24.05.2016 ***\n')


# Read and test command line arguments and set 'ParameterFile' and 'p2s' accordingly
# ----------------------------------------------------------------------------------
# test for correct number of command line arguments and abort if wrong number if found
if ( len(sys.argv)<2 or len(sys.argv)>3 ):
    print('Wrong number of command line arguments!')
    print('usage: RawDataProcessing <pf> [p2s]')
    print """
    pf:  name of parameter file (compulsory)
    p2s: print detailed output to standard out (optional)
           1 = yes (default)
           0 = no
    """
    sys.exit(-1)
    
# set 'p2s' either according to to command line argument or to default, if forbidden
# values for 'p2s' are found
p2s = 1    # default
if ( len(sys.argv)==3 ):
    try:
        p2s = int(sys.argv[2])
        if (p2s!=0) and (p2s!=1):
            print('Warning: Forbidden value \'' + sys.argv[2] + '\' for [p2s] found. [p2s] is set to default.')
            p2s = 1
    except ValueError:
        print('Warning: Forbidden character \'' + sys.argv[2] + '\' for [p2s] found. [p2s] is set to default.')

# try to read parameter file and abort if not existing
ParameterFile = sys.argv[1]
try:
    pf = open(ParameterFile, 'r')
    pf.close()
except:
    print('Error: Could not open parameter File \'' + sys.argv[1] + '\'. Abort.')
    sys.exit(-1)


# set predefined variables
# ------------------------
MeasureMode = 'scene'  # measurement mode
L_FC        = 0.2      # transmission loss of FC [dB]


# Read the instrument settings and evaluation parameters from the parameter file
# and print parameter values to screen 
# ------------------------------------------------------------------------------
print('***************************************************************\n')
print('Read parameter file \'' + ParameterFile + '\'.')
# open parameter file and extract dictionary with keywords and values
EP = ReadParameterFile(ParameterFile, MeasureMode)

# extract the instrument settings from the dictionary
T_set = EP['t_set']              # set-point temperature [K]

# extract evaluation parameters from the dictionary
MoSNoRFI    = EP['mosnorfi']     # mean/std expected for Gaussian distribution
MoSDiffMax  = EP['mosdiffmax']   # maximum allowed deviation of mean/std from MoSNoRFI
SKDiffMax   = EP['skdiffmax']    # maximum allowed deviation of skewness from 0
KUDiffMax   = EP['kudiffmax']    # maximum allowed deviation of kurtosis from 3
FDMax       = EP['fdmax']        # maximum allowed channel difference [K]
TempDiffMax = EP['tempdiffmax']  # maximum allowed difference between T_0 and T_set [K]
TB_HPMin    = EP['tb_hpmin']     # minimum acceptable brightness temperature at HP [K]
TB_HPMax    = EP['tb_hpmax']     # maximum acceptable brightness temperature at HP [K]
TB_VPMin    = EP['tb_vpmin']     # maximum acceptable brightness temperature at VP [K]
TB_VPMax    = EP['tb_vpmax']     # minimum acceptable brightness temperature at VP [K]
# If 'TB_pMax' is set to '1', the air temperature 'T_ext' measured simultaenously with
# the brightness temperatures is used as upper limit for the plausibility check.
# Otherwise, the value given in 'TB_pMax' is used.)

# extract other parameters from the dictionary
Directory   = EP['dir2evaluate'] # directory with files to evaluate
ACSTempFile = EP['acstempfile']  # file with time series of ACS temperatures
OutFileName = EP['outfilename']  # file name for output of results

L_FC = EP['trans_loss']
if L_FC == '':
  L_FC        = 0.2    # transmission loss of FC [dB]
  print('Warning: Transmission loss of Feed Cable not provided, default of L_FC: %.2f dB assumed'%(L_FC))



# Print instrument setting and evaluation parameters to screen
# ------------------------------------------------------------
print('\nEvaluating ELBARAII data files in directory ' + Directory + '/')
print('with the following instrument settings:')
print('     set-point temperature:       ' + str(T_set) + ' K')
print('and evaluation parameters:')
print('     ACS temperature file:        ' + ACSTempFile)
print('     MoS expected for Gaussian:   ' + str(MoSNoRFI))
print('     maximum allowed MoSDiff:     ' + str(MoSDiffMax))
print('     maximum allowed SkewDiff:    ' + str(SKDiffMax))
print('     maximum allowed KurtDiff:    ' + str(KUDiffMax))
print('     maximum allowed FD:          ' + str(FDMax) + ' K')
print('     maximum allowed SysTempDiff: ' + str(TempDiffMax) + ' K')
if (TB_HPMax != 1):
    print('     allowed range of TB_HP:      ' + str(TB_HPMin) + ' - ' + str(TB_HPMax) + ' K')
else:
    print('     allowed range of TB_HP:      ' + str(TB_HPMin) + ' K - T_air')
if (TB_VPMax != 1):
    print('     allowed range of TB_VP:      ' + str(TB_VPMin) + ' - ' + str(TB_VPMax) + ' K')
else:
    print('     allowed range of TB_VP:      ' + str(TB_VPMin) + ' K - T_air')
print('     detailed output to stdout:   ' + str(p2s))
print('     name basis for output files: ' + OutFileName)
print('\n***************************************************************')


# *************************************************************************
# * From here on, the actual evaluation program starts. Change only, if   *
# * you know what you are doing or if you really want to, of course ;- )  *
# *************************************************************************

# Load time series of ACS temperatures derived from the calibration
# measurements with 'ACSCalibration.py'
# -----------------------------------------------------------------
print('')
print('Load ACS temperature time series from ' + ACSTempFile)

ACSTimeSeries   = load(ACSTempFile)                               # load ACS temperature time series
Flagged         = nonzero(ACSTimeSeries[:,6])                     # find row number of flagged T_ACS values
ACSTimeSeries   = delete(ACSTimeSeries, Flagged, axis=0)          # remove rows with flagged T_ACS values
T_ACS_Dates     = DatetimeObject2DatetimeNum(ACSTimeSeries[:,0])  # extract dates and convert to numerical format
T_ACS           = ACSTimeSeries[:,1].tolist()                     # extract T_ACS values and convert to list
T_ACSHP1        = ACSTimeSeries[:,2].tolist()                     # extract T_ACS values and convert to list
T_ACSHP2        = ACSTimeSeries[:,3].tolist()                     # extract T_ACS values and convert to list
T_ACSVP1        = ACSTimeSeries[:,4].tolist()                     # extract T_ACS values and convert to list
T_ACSVP2        = ACSTimeSeries[:,5].tolist()                     # extract T_ACS values and convert to list

print('     ' + str(len(T_ACS)) + ' ACS temperatures between ' + ACSTimeSeries[0,0].strftime("%d-%b-%Y %H:%M:%S") + ' and ' + ACSTimeSeries[len(T_ACS)-1,0].strftime("%d-%b-%Y %H:%M:%S") + ' found.')
print('')
print('***************************************************************')


# Initialize arrays for output
# ----------------------------
TempOut = ones( (0,11) )
VoltOut = ones( (0,41) )


# Create list of files to evaluate
# --------------------------------
FilesToEvaluate = os.listdir(Directory)
NrOfFiles       = len(FilesToEvaluate)


# Do loop over all files to evaluate
# ----------------------------------
for nr in range(0, NrOfFiles):
    print('')
    print('Am now processing ' + Directory + '/' + FilesToEvaluate[nr])


    # 1.) Read data from file and split into individual variables
    # -----------------------------------------------------------
    LoadsMeasured, n_Loads, n_Samples, DatesMeasured, T_plate, T_ext, ElevationAngle, nps, Voltages = ReadDataFile(Directory + "/" + FilesToEvaluate[nr])
    

    # 2.) Perform some plausibility checks for the elevation angle
    # ------------------------------------------------------------
    PlausibilityChecks(ElevationAngle, MeasureMode)


    # 3.) Split Voltages into voltages per noise sourece and calculate statistics
    # ---------------------------------------------------------------------------
    U_HL1_Mean,  U_HL1_STD,  U_HL1_MoS,  U_HL1_SK,  U_HL1_KU  = StatsOverRecord(Voltages[: ,     0: nps])
    U_HL2_Mean,  U_HL2_STD,  U_HL2_MoS,  U_HL2_SK,  U_HL2_KU  = StatsOverRecord(Voltages[: ,   nps: 2*nps])
    U_ACS1_Mean, U_ACS1_STD, U_ACS1_MoS, U_ACS1_SK, U_ACS1_KU = StatsOverRecord(Voltages[: , 2*nps: 3*nps])
    U_ACS2_Mean, U_ACS2_STD, U_ACS2_MoS, U_ACS2_SK, U_ACS2_KU = StatsOverRecord(Voltages[: , 3*nps: 4*nps])
    U_RS1_Mean,  U_RS1_STD,  U_RS1_MoS,  U_RS1_SK,  U_RS1_KU  = StatsOverRecord(Voltages[: , 4*nps: 5*nps])
    U_RS2_Mean,  U_RS2_STD,  U_RS2_MoS,  U_RS2_SK,  U_RS2_KU  = StatsOverRecord(Voltages[: , 5*nps: 6*nps])
    U_HP1_Mean,  U_HP1_STD,  U_HP1_MoS,  U_HP1_SK,  U_HP1_KU  = StatsOverRecord(Voltages[: , 6*nps: 7*nps])
    U_HP2_Mean,  U_HP2_STD,  U_HP2_MoS,  U_HP2_SK,  U_HP2_KU  = StatsOverRecord(Voltages[: , 7*nps: 8*nps])
    U_VP1_Mean,  U_VP1_STD,  U_VP1_MoS,  U_VP1_SK,  U_VP1_KU  = StatsOverRecord(Voltages[: , 8*nps: 9*nps])
    U_VP2_Mean,  U_VP2_STD,  U_VP2_MoS,  U_VP2_SK,  U_VP2_KU  = StatsOverRecord(Voltages[: , 9*nps:10*nps])


    # 4.) RFI analysis in the time domain
    # -----------------------------------
    # no RFI analysis in the TD is carried out, when less than 10 values are found
    # in one data sample
    if (nps < 10):
        KUFlag_HP  = ones( (n_Samples,1) ) * -99
        KUFlag_VP  = ones( (n_Samples,1) ) * -99
        SKFlag_HP  = ones( (n_Samples,1) ) * -99
        SKFlag_VP  = ones( (n_Samples,1) ) * -99
        MoSFlag_HP = ones( (n_Samples,1) ) * -99
        MoSFlag_VP = ones( (n_Samples,1) ) * -99
        print("  No RFI analysis in the time domain possible (nps < 10). All TD Flags are set to -99.")
    else:
        # 4.1) RFI analysis by means of the kurtosis criterion
        # H-pol.
        KUFlag_HP1 = CheckValuePlausibility(-KUDiffMax, KUDiffMax, U_HP1_KU)
        KUFlag_HP2 = CheckValuePlausibility(-KUDiffMax, KUDiffMax, U_HP2_KU)
        KUFlag_HP  = KUFlag_HP1 + KUFlag_HP2; KUFlag_HP[KUFlag_HP>0] = 1
        # V-pol.
        KUFlag_VP1 = CheckValuePlausibility(-KUDiffMax, KUDiffMax, U_VP1_KU)
        KUFlag_VP2 = CheckValuePlausibility(-KUDiffMax, KUDiffMax, U_VP2_KU)
        KUFlag_VP  = KUFlag_VP1 + KUFlag_VP2; KUFlag_VP[KUFlag_VP>0] = 1
    
        # 4.2) RFI analysis by means of the skewness criterion
        # H-pol.
        SKFlag_HP1 = CheckValuePlausibility(-SKDiffMax, SKDiffMax, U_HP1_SK)
        SKFlag_HP2 = CheckValuePlausibility(-SKDiffMax, SKDiffMax, U_HP2_SK)
        SKFlag_HP  = SKFlag_HP1 + SKFlag_HP2; SKFlag_HP[SKFlag_HP>0] = 1
        # V-pol.
        SKFlag_VP1 = CheckValuePlausibility(-SKDiffMax, SKDiffMax, U_VP1_SK)
        SKFlag_VP2 = CheckValuePlausibility(-SKDiffMax, SKDiffMax, U_VP2_SK)
        SKFlag_VP  = SKFlag_VP1 + SKFlag_VP2; SKFlag_VP[SKFlag_VP>0] = 1

        # 4.3) RFI analysisis by means of the mean over sigma criterion
        # H-pol.
        MoSFlag_HP1 = CheckValuePlausibility(MoSNoRFI-MoSDiffMax, MoSNoRFI+MoSDiffMax, U_HP1_MoS)
        MoSFlag_HP2 = CheckValuePlausibility(MoSNoRFI-MoSDiffMax, MoSNoRFI+MoSDiffMax, U_HP2_MoS)
        MoSFlag_HP  = MoSFlag_HP1 + MoSFlag_HP2; MoSFlag_HP[MoSFlag_HP>0] = 1
        # V-pol.
        MoSFlag_VP1 = CheckValuePlausibility(MoSNoRFI-MoSDiffMax, MoSNoRFI+MoSDiffMax, U_VP1_MoS)
        MoSFlag_VP2 = CheckValuePlausibility(MoSNoRFI-MoSDiffMax, MoSNoRFI+MoSDiffMax, U_VP2_MoS)
        MoSFlag_VP  = MoSFlag_VP1 + MoSFlag_VP2; MoSFlag_VP[MoSFlag_VP>0] = 1


    # 5.) Calculate calibrated brightness temperatures for both polarizations
    #     and frequency channels from the radiometer raw data (voltages and
    #     physical temperatures) and the ACS temperature time series
    # -----------------------------------------------------------------------
    # 5.1) Derive ACS temperature for date and time of measurement from ACS
    #      temperature time series
    # compute interpolated T_ACS values
    T_ACS_interp = interp(DatetimeObject2DatetimeNum(DatesMeasured), T_ACS_Dates, T_ACS, left=T_ACS[0], right=T_ACS[len(T_ACS)-1])
    T_ACSHP1_interp = interp(DatetimeObject2DatetimeNum(DatesMeasured), T_ACS_Dates, T_ACSHP1, left=T_ACSHP1[0], right=T_ACSHP1[len(T_ACSHP1)-1])
    T_ACSHP2_interp = interp(DatetimeObject2DatetimeNum(DatesMeasured), T_ACS_Dates, T_ACSHP2, left=T_ACSHP2[0], right=T_ACSHP2[len(T_ACSHP2)-1])
    T_ACSVP1_interp = interp(DatetimeObject2DatetimeNum(DatesMeasured), T_ACS_Dates, T_ACSVP1, left=T_ACSVP1[0], right=T_ACSVP1[len(T_ACSVP1)-1])
    T_ACSVP2_interp = interp(DatetimeObject2DatetimeNum(DatesMeasured), T_ACS_Dates, T_ACSVP2, left=T_ACSVP2[0], right=T_ACSVP2[len(T_ACSVP2)-1])
    # convert to (n_Samples x 1) array
    T_ACS_interp = transpose([T_ACS_interp])
    T_ACSHP1_interp = transpose([T_ACSHP1_interp])
    T_ACSHP2_interp = transpose([T_ACSHP2_interp])
    T_ACSVP1_interp = transpose([T_ACSVP1_interp])
    T_ACSVP2_interp = transpose([T_ACSVP2_interp])
    # print warning, if 'DatesMeasured' are outside time period for which ACS temperatures are available
    if (min(DatetimeObject2DatetimeNum(DatesMeasured)) < T_ACS_Dates[0]) or (max(DatetimeObject2DatetimeNum(DatesMeasured)) > T_ACS_Dates[len(T_ACS)-1]):
        print('    Warning: At least one sample is outside the time period for which ACS temperatures are available. The nearest available ACS temperature is used instead of an interpolated value.')

    # 5.2) Calculate brightness temperatures 'TB_p_in' received at radiometer
    #      input port
    TB_HP1_in = U2TB_in(T_plate, T_ACSHP1_interp, U_RS1_Mean, U_ACS1_Mean, U_HP1_Mean)
    TB_HP2_in = U2TB_in(T_plate, T_ACSHP2_interp, U_RS2_Mean, U_ACS2_Mean, U_HP2_Mean)
    TB_VP1_in = U2TB_in(T_plate, T_ACSVP1_interp, U_RS1_Mean, U_ACS1_Mean, U_VP1_Mean)
    TB_VP2_in = U2TB_in(T_plate, T_ACSVP2_interp, U_RS2_Mean, U_ACS2_Mean, U_VP2_Mean)
    T_HL1 = U2TB_in(T_plate, T_ACSVP1_interp, U_RS1_Mean, U_ACS1_Mean, U_HL1_Mean)
    T_HL2 = U2TB_in(T_plate, T_ACSVP2_interp, U_RS2_Mean, U_ACS2_Mean, U_HL2_Mean)


    # 5.3) Correct for cable losses to derive  brightness temperature 'TB_p'
    #      entering the antenna aperture
    TB_HP1 = TB_in2TB(TB_HP1_in, T_ext, L_FC)
    TB_HP2 = TB_in2TB(TB_HP2_in, T_ext, L_FC)
    TB_VP1 = TB_in2TB(TB_VP1_in, T_ext, L_FC)
    TB_VP2 = TB_in2TB(TB_VP2_in, T_ext, L_FC)


    # 6.) RFI analysis in the frequency domain
    # ----------------------------------------
    FD_HP, FDFlag_HP = CheckTempDiff(TB_HP1, TB_HP2, FDMax)
    FD_VP, FDFlag_VP = CheckTempDiff(TB_VP1, TB_VP2, FDMax)


    # 7.) Calculate mean brightness temperature of both frequency channels
    # --------------------------------------------------------------------
    TB_HP = (TB_HP1 + TB_HP2) / 2
    TB_VP = (TB_VP1 + TB_VP2) / 2


    # 8.) Perform further checks for the plausibility and potential corruption
    #     of the derived brightness temperatures and set corresponding flags
    # ------------------------------------------------------------------------
    # 8.1) Calculate deviation of internal system temperature from set-point
    #      temperature and set corresponding flags
    SysTempDiff, SysTempFlag = CheckTempDiff(T_set, T_plate, TempDiffMax)
    
    # 8.2) Compare derived brightness temperatures to minimum and maximum
    #      acceptable values 'TBmin' and 'TBmax' and set corresponding flags
    # H-polarization
    if (TB_HPMax == 1):
        TBFlag_HP = CheckValuePlausibility(TB_HPMin, T_ext, TB_HP)
    else:
        TBFlag_HP = CheckValuePlausibility(TB_HPMin, TB_HPMax, TB_HP)
    # V-polarization
    if (TB_VPMax == 1):
        TBFlag_VP = CheckValuePlausibility(TB_VPMin, T_ext, TB_VP)
    else:
        TBFlag_VP = CheckValuePlausibility(TB_VPMin, TB_VPMax, TB_VP)


    # 9.) Combine all flags into unambigous value 'AllFlag'
    # -----------------------------------------------------
    AllFlag_HP = CombineFlags(KUFlag_HP, SKFlag_HP, MoSFlag_HP, FDFlag_HP, TBFlag_HP, SysTempFlag)
    AllFlag_VP = CombineFlags(KUFlag_VP, SKFlag_VP, MoSFlag_VP, FDFlag_VP, TBFlag_VP, SysTempFlag)


    # 10.) Prepare output of results
    # ------------------------------
    # 10.1) Brightness temperatures, channel differences, flags, and physical
    #       temperatures
    # write results to 'NextRows' (to be appended to 'TempOut')
    NextRows = transpose([DatesMeasured])
    NextRows = append(NextRows, ElevationAngle, axis=1)
    NextRows = append(NextRows, TB_HP, axis=1)
    NextRows = append(NextRows, FD_HP, axis=1)
    NextRows = append(NextRows, AllFlag_HP, axis=1)
    NextRows = append(NextRows, TB_VP, axis=1)
    NextRows = append(NextRows, FD_VP, axis=1)
    NextRows = append(NextRows, AllFlag_VP, axis=1)
    NextRows = append(NextRows, T_plate, axis=1)    
    NextRows = append(NextRows, SysTempDiff, axis=1)
    NextRows = append(NextRows, T_ext, axis=1)    
    # append 'NextRows' to 'TempOut'
    TempOut  = append(TempOut, NextRows, axis=0)

    # 10.2) Voltage statistics for ACS, HS, HP, and VP
    # write results to 'NextRows' (to be appended to 'VoltOut')
    NextRows  = transpose([DatesMeasured])
    # append data for ACS, channel 1
    NextRows = append(NextRows, U_ACS1_Mean, axis=1)
    NextRows = append(NextRows, U_ACS1_STD, axis=1)
    NextRows = append(NextRows, U_ACS1_MoS, axis=1)
    NextRows = append(NextRows, U_ACS1_SK, axis=1)
    NextRows = append(NextRows, U_ACS1_KU, axis=1)
    # append data for ACS, channel 2
    NextRows = append(NextRows, U_ACS2_Mean, axis=1)
    NextRows = append(NextRows, U_ACS2_STD, axis=1)
    NextRows = append(NextRows, U_ACS2_MoS, axis=1)
    NextRows = append(NextRows, U_ACS2_SK, axis=1)
    NextRows = append(NextRows, U_ACS2_KU, axis=1)
    # append data for RS, channel 1
    NextRows = append(NextRows, U_RS1_Mean, axis=1)
    NextRows = append(NextRows, U_RS1_STD, axis=1)
    NextRows = append(NextRows, U_RS1_MoS, axis=1)
    NextRows = append(NextRows, U_RS1_SK, axis=1)
    NextRows = append(NextRows, U_RS1_KU, axis=1)
    # append data for RS, channel 2
    NextRows = append(NextRows, U_RS2_Mean, axis=1)
    NextRows = append(NextRows, U_RS2_STD, axis=1)
    NextRows = append(NextRows, U_RS2_MoS, axis=1)
    NextRows = append(NextRows, U_RS2_SK, axis=1)
    NextRows = append(NextRows, U_RS2_KU, axis=1)
    # append data for HP, channel 1
    NextRows = append(NextRows, U_HP1_Mean, axis=1)
    NextRows = append(NextRows, U_HP1_STD, axis=1)
    NextRows = append(NextRows, U_HP1_MoS, axis=1)
    NextRows = append(NextRows, U_HP1_SK, axis=1)
    NextRows = append(NextRows, U_HP1_KU, axis=1)
    # append data for HP, channel 2
    NextRows = append(NextRows, U_HP2_Mean, axis=1)
    NextRows = append(NextRows, U_HP2_STD, axis=1)
    NextRows = append(NextRows, U_HP2_MoS, axis=1)
    NextRows = append(NextRows, U_HP2_SK, axis=1)
    NextRows = append(NextRows, U_HP2_KU, axis=1)
    # append data for VP, channel 1
    NextRows = append(NextRows, U_VP1_Mean, axis=1)
    NextRows = append(NextRows, U_VP1_STD, axis=1)
    NextRows = append(NextRows, U_VP1_MoS, axis=1)
    NextRows = append(NextRows, U_VP1_SK, axis=1)
    NextRows = append(NextRows, U_VP1_KU, axis=1)
    # append data for VP, channel 2
    NextRows = append(NextRows, U_VP2_Mean, axis=1)
    NextRows = append(NextRows, U_VP2_STD, axis=1)
    NextRows = append(NextRows, U_VP2_MoS, axis=1)
    NextRows = append(NextRows, U_VP2_SK, axis=1)
    NextRows = append(NextRows, U_VP2_KU, axis=1)
    # append 'NextRows' to 'VoltOut'
    VoltOut = append(VoltOut, NextRows, axis=0)


    # 11.) Print results of the computations, RFI and plausibility checks to
    #      screen, if p2s is set. (No calculations are carried out here.)
    # ----------------------------------------------------------------------
    if (p2s == 1):
        # Some basic information about the data file and the corresponding measurements
        print('     date of measurement:           %02d.%02d.%4d %02d:%02d:%02d to %02d.%02d.%4d %02d:%02d:%02d'%(
	     DatesMeasured[0].day, DatesMeasured[0].month, DatesMeasured[0].year, DatesMeasured[0].hour, DatesMeasured[0].minute, DatesMeasured[0].second,
             DatesMeasured[n_Samples-1].day, DatesMeasured[n_Samples-1].month, DatesMeasured[n_Samples-1].year, DatesMeasured[n_Samples-1].hour, DatesMeasured[n_Samples-1].minute, DatesMeasured[n_Samples-1].second))
        print('     nr. of samples:                ' + str(n_Samples))
        print('     values per sample (nps):       ' + str(nps))
        print('     nr. of measured loads:         ' + str(n_Loads))
        print('     elevation angle(s):            %d to %d deg.'%( amin(ElevationAngle), amax(ElevationAngle) ))
        # Voltage statistics
        print('     voltage statistics: mean, standard deviation, mean/std, skewness, kurtosis (mean of all samples)')
        if (nps < 10):
            print('     (Standard deviation, skewness, and kurtosis were set to -99, since less than 10 values per noise source were found.')
        print('        ACS1:                       %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_ACS1_Mean), mean(U_ACS1_STD), mean(U_ACS1_MoS), mean(U_ACS1_SK), mean(U_ACS1_KU)))
        print('        ACS2:                       %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_ACS2_Mean), mean(U_ACS2_STD), mean(U_ACS2_MoS), mean(U_ACS2_SK), mean(U_ACS2_KU)))
        print('        RS1:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_RS1_Mean),  mean(U_RS1_STD),  mean(U_RS1_MoS),  mean(U_RS1_SK),  mean(U_RS1_KU)))
        print('        RS2:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_RS2_Mean),  mean(U_RS2_STD),  mean(U_RS2_MoS),  mean(U_RS2_SK),  mean(U_RS2_KU)))
        print('        HL1:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HL1_Mean),  mean(U_HL1_STD),  mean(U_HL1_MoS),  mean(U_HL1_SK),  mean(U_HL1_KU)))
        print('        HL2:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HL2_Mean),  mean(U_HL2_STD),  mean(U_HL2_MoS),  mean(U_HL2_SK),  mean(U_HL2_KU)))
        print('        HP1:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HP1_Mean),  mean(U_HP1_STD),  mean(U_HP1_MoS),  mean(U_HP1_SK),  mean(U_HP1_KU)))
        print('        HP2:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HP2_Mean),  mean(U_HP2_STD),  mean(U_HP2_MoS),  mean(U_HP2_SK),  mean(U_HP2_KU)))
        print('        VP1:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_VP1_Mean),  mean(U_VP1_STD),  mean(U_VP1_MoS),  mean(U_VP1_SK),  mean(U_VP1_KU)))
        print('        VP2:                        %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_VP2_Mean),  mean(U_VP2_STD),  mean(U_VP2_MoS),  mean(U_VP2_SK),  mean(U_VP2_KU)))
        # Interpolated ACS temperatures
        print('     interpolated ACS temperatures: %.2f K - %.2f K'%(amin(T_ACS_interp), amax(T_ACS_interp)))
        print('     interpolated HL temperatures: ')
        print('        T_HL, channel 1:            %.2f K - %.2f K'%(amin(T_HL1), amax(T_HL1)))
        print('        T_HL, channel 2:            %.2f K - %.2f K'%(amin(T_HL2), amax(T_HL2)))
        print('     brightness temperatures at radiometer input port')
        print('        HP, channel 1:              %.2f K - %.2f K'%(amin(TB_HP1_in), amax(TB_HP1_in)))
        print('        HP, channel 2:              %.2f K - %.2f K'%(amin(TB_HP2_in), amax(TB_HP2_in)))
        print('        VP, channel 1:              %.2f K - %.2f K'%(amin(TB_VP1_in), amax(TB_VP1_in)))
        print('        VP, channel 2:              %.2f K - %.2f K'%(amin(TB_VP2_in), amax(TB_VP2_in)))
        # Brightness temperatures at antenna aperture
        print('     brightness temperatures at antenna aperture')
        print('        HP, channel 1:              %.2f K - %.2f K'%(amin(TB_HP1), amax(TB_HP1)))
        print('        HP, channel 2:              %.2f K - %.2f K'%(amin(TB_HP2), amax(TB_HP2)))
        print('        VP, channel 1:              %.2f K - %.2f K'%(amin(TB_VP1), amax(TB_VP1)))
        print('        VP, channel 2:              %.2f K - %.2f K'%(amin(TB_VP2), amax(TB_VP2)))
        # RFI analysis in the the time domain
        print('     RFI analysis in the time domain (kurtosis, skewness, mean/std)')
        if (nps > 10):
            print('        nr. of KUFlags for HP:      ' + str(size(flatnonzero(KUFlag_HP))))
            print('        nr. of KUFlags for VP:      ' + str(size(flatnonzero(KUFlag_VP))))
            print('        nr. of SKFlags for HP:      ' + str(size(flatnonzero(SKFlag_HP))))
            print('        nr. of SKFlags for VP:      ' + str(size(flatnonzero(SKFlag_VP))))
            print('        nr. of MoSFlags for HP:     ' + str(size(flatnonzero(MoSFlag_HP))))
            print('        nr. of MoSFlags for VP:     ' + str(size(flatnonzero(MoSFlag_VP))))
        else:
            print('        No RFI analysis in the time domain possible (nps < 10). All TD flags are set to -99.')
        # RFI Analysis in the frequency domain
        print('     RFI analysis in the frequency domain')
        print('        FD at HP:                   %.2f K - %.2f K'%(amin(FD_HP), amax(FD_HP)))
        print('        nr. of FDFlags at HP:       '  + str(size(flatnonzero(FDFlag_HP))))
        print('        FD at VP:                   %.2f K - %.2f K'%(amin(FD_VP), amax(FD_VP)))
        print('        nr. of FDFlags at VP:       ' + str(size(flatnonzero(FDFlag_VP))))
        # System temperature check
        print('     system temperature check')
        print('        T_set - T0:                 %.2f K - %.2f K'%(amin(SysTempDiff), amax(SysTempDiff)))
        print('        nr. of SysTempFlags:        ' + str(size(flatnonzero(SysTempFlag))))
        # Brightness temperature plausibility check
        print('     brightness temperature plausibility check');
        print('        nr. of TBFlags at HP:       ' + str(size(flatnonzero(TBFlag_HP))))
        print('        nr. of TBFlags at VP:       ' + str(size(flatnonzero(TBFlag_VP))))


# 12.) Save 'TempOut' and 'VoltOut' to file
# -----------------------------------------
print('\n***************************************************************\n')

# 12.1) save 'TempOut' and 'VoltOut' to numpy binary files
print('Saving results to ' + OutFileName + 'Temp.npy and ' + OutFileName + 'Volt.npy.')
save(OutFileName + 'Temp', TempOut)
save(OutFileName + 'Volt', VoltOut)

# 12.2) save 'TempOut' and 'VoltOut' to text file
print('Saving results to ' + OutFileName + 'Temp.txt and ' + OutFileName + 'Volt.txt.')
WriteResults2Text(OutFileName, TempOut, VoltOut)

# 13.) Print plot time series of results ('TempOut')
# --------------------------------------------------
PlotTimeSeries(TempOut)
