#!/usr/bin/python
#-*-coding: utf-8-*-
"""
Program to process ELBARAII raw data of sky calibration measurements to
calibrated ACS temperatures needed as input for 'RawDataProcessing'.
Included are also RFI analysis in the time and frequency domain, as well
as checks for the plausibility of the derived ACS temperatures and for
problems with the thermal stabilization of the instrument during the
measurements.

usage: ACSCalibration <pf> [p2s]
        pf:  name of parameter file (compulsory)
        p2s: print detailed output to standard out (optional)
              1 = yes (default)
              0 = no

author: Ingo Voelksch, 31.10.2011
        Andreas Wiesmann, 20160524 changed to independent calibration per channel
"""

from libdataproc import *
import os

print('*** Processing ELBARAII raw data of sky measurements to ACS temperatures ***')
print('*** author: Ingo Voelksch, Andreas Wiesmann 24.05.2016 ***\n')


# Read and test command line arguments and set 'ParameterFile' and 'p2s' accordingly
# ----------------------------------------------------------------------------------
# test for correct number of command line arguments and abort if wrong number if found
if ( len(sys.argv)<2 or len(sys.argv)>3 ):
    print('Wrong number of command line arguments!')
    print('usage: ACSCalibration <pf> [p2s]')
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
MeasureMode = 'sky'  # measurement mode


# Read the instrument settings and evaluation parameters from the parameter file
# and print parameter values to screen 
# ------------------------------------------------------------------------------
print('***************************************************************\n')
print('Read parameter file ' + ParameterFile + '.')
# open parameter file and extract dictionary with keywords and values
EP = ReadParameterFile(ParameterFile, MeasureMode)

# extract the instrument settings from the dictionary
Z     = EP['altitude']           # altitude of ELBARAII [km]
T_set = EP['t_set']              # set-point temperature [K]

# extract evaluation parameters from the dictionary
MoSNoRFI    = EP['mosnorfi']     # mean/std expected for Gaussian distribution
MoSDiffMax  = EP['mosdiffmax']   # maximum allowed deviation of mean/std from MoSNoRFI
SKDiffMax   = EP['skdiffmax']    # maximum allowed deviation of skewness from 0
KUDiffMax   = EP['kudiffmax']    # maximum allowed deviation of kurtosis from 3
FDmax       = EP['fdmax']        # maximum allowed channel difference [K]
TempDiffMax = EP['tempdiffmax']  # maximum allowed difference between T_0 and T_set [K]
T_ACSmin    = EP['t_acsmin']     # minimum acceptable ACS temperature [K]
T_ACSmax    = EP['t_acsmax']     # maximum acceptable ACS temperature [K]

# extract other parameters from the dictionary
Directory   = EP['dir2evaluate'] # directory with files to evaluate
OutFileName = EP['outfilename']  # file name for output of results

L_FC = EP['trans_loss']
if L_FC == '':
  L_FC        = 0.2    # transmission loss of FC [dB]
  print('Warning: Transmission loss of Feed Cable not provided, default of L_FC: %.2f dB assumed'%(L_FC))

# Print instrument setting and evaluation parameters to screen
# ------------------------------------------------------------
print('\nEvaluating ELBARAII data files in directory ' + Directory + '/')
print('with the following instrument settings:')
print('     Feed Cable loss:             ' + str(L_FC) + ' dB')
print('     altitude:                    ' + str(Z) + ' km')
print('     set-point temperature:       ' + str(T_set) + ' K')
print('and evaluation parameters:')
print('     MoS expected for Gaussian:   ' + str(MoSNoRFI))
print('     maximum allowed MoSDiff:     ' + str(MoSDiffMax))
print('     maximum allowed SkewDiff:    ' + str(SKDiffMax))
print('     maximum allowed KurtDiff:    ' + str(KUDiffMax))
print('     maximum allowed FD:          ' + str(FDmax) + ' K')
print('     maximum allowed SysTempDiff: ' + str(TempDiffMax) + ' K')
print('     allowed range of T_ACS:      ' + str(T_ACSmin) + '-' + str(T_ACSmax) + ' K')
print('     detailed output to stdout    ' + str(p2s))
print('     name basis for output files: ' + OutFileName)
print('\n***************************************************************')


# *************************************************************************
# * From here on, the actual evaluation program starts. Change only, if   *
# * you know what you are doing or if you really want to, of course ;- )  *
# *************************************************************************

# Create list of files to evaluate
# --------------------------------
FilesToEvaluate = sorted(os.listdir(Directory))
NrOfFiles       = len(FilesToEvaluate)


# Initialize arrays for output
# ----------------------------
TACSOut = ones( (0,7) )
TempOut = ones( (0,11) )
VoltOut = ones( (0,41) )


# Do loop over all files to evaluate
# ----------------------------------
for nr in range(0, NrOfFiles):
    print('')
    print("Now processing " + Directory + "/" + FilesToEvaluate[nr])


    # 1.) Read data from file and split into individual variables
    # -----------------------------------------------------------
    LoadsMeasured, n_Loads, n_Samples, DatesMeasured, T_plate, T_ext, ElevationAngle, nps, Voltages = ReadDataFile(Directory + "/" + FilesToEvaluate[nr])


    # 2.) Perform some plausibility checks for the elevation angle
    # ------------------------------------------------------------
    PlausibilityChecks(ElevationAngle, MeasureMode)


    # 3.) Split Voltages into Voltages per noise source and calculate statistics
    # --------------------------------------------------------------------------
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
        KUFlag  = ones( (n_Samples,1) ) * -99
        SKFlag  = ones( (n_Samples,1) ) * -99
        MoSFlag = ones( (n_Samples,1) ) * -99
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

        # 4.4) If one polarization is flagged as corrupted by one of the three
        #      checks, the corresponding TD flag is principally set, since the
        #      mean voltage of both polarizations is used in the further
        #      calculations (see next calculation step)
        KUFlag  = KUFlag_HP + KUFlag_VP; KUFlag[KUFlag>0] = 1
        SKFlag  = SKFlag_HP + SKFlag_VP; SKFlag[SKFlag>0] = 1
        MoSFlag = MoSFlag_HP + MoSFlag_VP; MoSFlag[MoSFlag>0] = 1


    # 5.) Calculation of ACS temperature from voltages measured
    # ---------------------------------------------------------
    # 5.1) Calculate mean of H- and V-pol per noise source and store in 'U_skyin':
    U_skyin1 = (U_HP1_Mean + U_VP1_Mean) / 2
    U_skyin2 = (U_HP2_Mean + U_VP2_Mean) / 2

    # 5.2) Calculate theoretical sky brightness 'T_sky' entering the antenna aper-
    #      ture with the approach of   from the radiometer
    #      altitute 'Z', the external air temperature 'T_ext', and the zenith
    #      angle, which is calculated as (180-ElevationAngle)
    T_sky = TSkyPellarin(Z, T_ext, (180-ElevationAngle))

    # 5.3) Correct for cable losses to derive the brightness 'T_skyin' which is
    #      actually received at the radiometer input port
    T_skyin = TB2TB_in(T_sky, T_ext, L_FC)
    

    # 5.4) Calculate the ACS temperatures for both channels from the sky measurements
    T_ACS1 = U2TB_in(T_plate, T_skyin, U_RS1_Mean, U_skyin1, U_ACS1_Mean)
    T_ACS2 = U2TB_in(T_plate, T_skyin, U_RS2_Mean, U_skyin2, U_ACS2_Mean)
    T_ACSHP1 = U2TB_in(T_plate, T_skyin, U_RS1_Mean, U_HP1_Mean, U_ACS1_Mean)
    T_ACSHP2 = U2TB_in(T_plate, T_skyin, U_RS2_Mean, U_HP2_Mean, U_ACS2_Mean)
    T_ACSVP1 = U2TB_in(T_plate, T_skyin, U_RS1_Mean, U_VP1_Mean, U_ACS1_Mean)
    T_ACSVP2 = U2TB_in(T_plate, T_skyin, U_RS2_Mean, U_VP2_Mean, U_ACS2_Mean)

    T_HL1 = U2TB_in(T_plate, T_skyin, U_RS1_Mean, U_HP1_Mean, U_HL1_Mean)
    T_HL2 = U2TB_in(T_plate, T_skyin, U_RS2_Mean, U_HP2_Mean, U_HL2_Mean)

    # 6.) RFI analysis in the frequency domain
    # ----------------------------------------
    FD, FDFlag = CheckTempDiff(T_ACS1, T_ACS2, FDmax)


    # 7.) Calculate mean ACS temperature of both channels
    # ---------------------------------------------------
    T_ACS = (T_ACS1 + T_ACS2) / 2


    # 8.) Perform further checks for the plausibility of the derived ACS tempera-
    #     tures and potential corruption of the data and set corresponding flags
    # ---------------------------------------------------------------------------
    # 8.1) Calculate deviation of internal system temperature from set-point
    #      temperature and set corresponding flags
    SysTempDiff, SysTempFlag = CheckTempDiff(T_set, T_plate, TempDiffMax)

    # 8.2) Compare derived ACS temperature to minimum and maximum acceptable
    #      values 'T_ACSmin' and 'T_ACSmax' and set corresponding flags
    ACSTempFlag = CheckValuePlausibility(T_ACSmin, T_ACSmax, T_ACS)


    # 9.) Combine all Flags into one unambigious value
    # ------------------------------------------------
    AllFlag = CombineFlags(KUFlag, SKFlag, MoSFlag, FDFlag, ACSTempFlag, SysTempFlag)


    # 10.) Prepare output of results
    # -------------------------------
    # 10.1) Create time series T_ACS(time) needed later on as input for raw data
    #       processing: TACSOut = [DatesMeasured T_ACS Flag]
    #       The mean ACS temperature of all unflagged samples is written to
    #       'TempOut'. If no unflagged samples are found, the mean ACS temp.
    #       of all samples is calculated and 'Flag' is set to 1.
    # find indices, where AllFlag is either 0 or -999000 (i.e., unflagged samples)
    i = flatnonzero( any([AllFlag == 0, AllFlag == -999000], axis=0) )  
    # calculate mean ACS temperature
    if (size(i) == 0):  # i.e., no unflagged samples were found
        # calculate mean ACS temperature of all samples
        T_ACS_mean = mean(T_ACS, axis=0)
        T_ACSHP1_mean = mean(T_ACSHP1, axis=0)
        T_ACSHP2_mean = mean(T_ACSHP2, axis=0)
        T_ACSVP1_mean = mean(T_ACSVP1, axis=0)
        T_ACSVP2_mean = mean(T_ACSVP2, axis=0)
        # write results to 'NextRow' (to be appended to 'TACSOut')
        NextRow = [DatesMeasured[0], T_ACS_mean[0], T_ACSHP1_mean[0], T_ACSHP2_mean[0], T_ACSVP1_mean[0], T_ACSVP2_mean[0], 1]
    else:  #  i.e., at least one unflagged sample was found
        # calculate mean ACS temperature of all unflagged samples
        T_ACS_mean = mean(T_ACS[i,:], axis=0)
        T_ACSHP1_mean = mean(T_ACSHP1[i,:], axis=0)
        T_ACSHP2_mean = mean(T_ACSHP2[i,:], axis=0)
        T_ACSVP1_mean = mean(T_ACSVP1[i,:], axis=0)
        T_ACSVP2_mean = mean(T_ACSVP2[i,:], axis=0)
        # write results to 'NextRow' (to be appended to 'TACSOut')
        NextRow = [DatesMeasured[0], T_ACS_mean[0], T_ACSHP1_mean[0], T_ACSHP2_mean[0], T_ACSVP1_mean[0], T_ACSVP2_mean[0], 0]
    # append 'NextRow' to 'TACSOut'
    TACSOut = append(TACSOut, [NextRow], axis=0)

    # 10.2) Create detailed output 'TempOut' and 'VoltOut' for each sample.
    # 10.2.a) ACS temperatures, channel differences, flags, and physical
    #         temperatures
    # write results to 'NextRows' (to be appended to 'TempOut')
    NextRows = transpose([DatesMeasured])
    NextRows = append(NextRows, ElevationAngle, axis=1)
    NextRows = append(NextRows, T_ACS, axis=1)
    NextRows = append(NextRows, FD, axis=1)
    NextRows = append(NextRows, AllFlag, axis=1)
    NextRows = append(NextRows, T_ACS, axis=1)
    NextRows = append(NextRows, FD, axis=1)
    NextRows = append(NextRows, AllFlag, axis=1)
    NextRows = append(NextRows, T_plate, axis=1)
    NextRows = append(NextRows, SysTempDiff, axis=1)
    NextRows = append(NextRows, T_ext, axis=1)
    # append 'NextRows' to 'TempOut'
    TempOut  = append(TempOut, NextRows, axis=0)

    # 10.2.b) voltage statistics for ACS, HS, HP, and VP
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


    # 11.) Print results of the evaluations and the different plausibility
    #      checks to screen, if p2s is set
    # --------------------------------------------------------------------
    if p2s == 1:
        # some basic information about the data file and the corresponding measurements
        print("     date of measurement:         %02d.%02d.%4d %02d:%02d:%02d to %02d.%02d.%4d %02d:%02d:%02d"%(
	     DatesMeasured[0].day, DatesMeasured[0].month, DatesMeasured[0].year, DatesMeasured[0].hour, DatesMeasured[0].minute, DatesMeasured[0].second,
             DatesMeasured[n_Samples-1].day, DatesMeasured[n_Samples-1].month, DatesMeasured[n_Samples-1].year, DatesMeasured[n_Samples-1].hour, DatesMeasured[n_Samples-1].minute, DatesMeasured[n_Samples-1].second))
        print("     nr. of samples:              " + str(n_Samples))
        print("     values per sample (nps):     " + str(nps))
        print("     nr. of measured loads:       " + str(n_Loads))
        print("     elevation angle(s):          %d to %d deg."%( amin(ElevationAngle), amax(ElevationAngle) ))
        # voltage statistics
        print('     voltage statistics: mean, standard deviation, mean/std, skewness, kurtosis (mean of all samples)')
        if (nps < 10):
            print('     (Standard deviation, skewness, and kurtosis were set to -99, since less than 10 values per noise source were found.')
        print('        ACS1:                     %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_ACS1_Mean), mean(U_ACS1_STD), mean(U_ACS1_MoS), mean(U_ACS1_SK), mean(U_ACS1_KU)))
        print('        ACS2:                     %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_ACS2_Mean), mean(U_ACS2_STD), mean(U_ACS2_MoS), mean(U_ACS2_SK), mean(U_ACS2_KU)))
        print('        RS1:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_RS1_Mean),  mean(U_RS1_STD),  mean(U_RS1_MoS),  mean(U_RS1_SK),  mean(U_RS1_KU)))
        print('        RS2:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_RS2_Mean),  mean(U_RS2_STD),  mean(U_RS2_MoS),  mean(U_RS2_SK),  mean(U_RS2_KU)))
        print('        HL1:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HL1_Mean),  mean(U_HL1_STD),  mean(U_HL1_MoS),  mean(U_HL1_SK),  mean(U_HL1_KU)))
        print('        HL2:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HL2_Mean),  mean(U_HL2_STD),  mean(U_HL2_MoS),  mean(U_HL2_SK),  mean(U_HL2_KU)))
        print('        HP1:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HP1_Mean),  mean(U_HP1_STD),  mean(U_HP1_MoS),  mean(U_HP1_SK),  mean(U_HP1_KU)))
        print('        HP2:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_HP2_Mean),  mean(U_HP2_STD),  mean(U_HP2_MoS),  mean(U_HP2_SK),  mean(U_HP2_KU)))
        print('        VP1:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_VP1_Mean),  mean(U_VP1_STD),  mean(U_VP1_MoS),  mean(U_VP1_SK),  mean(U_VP1_KU)))
        print('        VP2:                      %.3f V, %.4f V, %.4f, %.4f, %.2f'%(mean(U_VP2_Mean),  mean(U_VP2_STD),  mean(U_VP2_MoS),  mean(U_VP2_SK),  mean(U_VP2_KU)))
        # theoretical sky brigthness at antenna aperture and input port, ACS temperatures
        print('     theoretical sky brightness at antenna aperture and input port, ACS temperatures')
        print('        T_sky:                    %.2f K - %.2f K'%(amin(T_sky), amax(T_sky)))
        print('        T_skyin:                  %.2f K - %.2f K'%(amin(T_skyin), amax(T_skyin)))
        print('        T_ACS, channel 1:         %.2f K - %.2fK'%(amin(T_ACS1), amax(T_ACS1)))
        print('        T_ACS, channel 2:         %.2f K - %.2fK'%(amin(T_ACS2), amax(T_ACS2)))
        print('        T_HL, channel 1:          %.2f K - %.2fK'%(amin(T_HL1), amax(T_HL1)))
        print('        T_HL, channel 2:          %.2f K - %.2fK'%(amin(T_HL2), amax(T_HL2)))
        # RFI analysis in the the time domain
        print('     RFI analysis in the time domain (kurtosis, skewness, mean/std)')
        if (nps > 10):
            print('        nr. of KUFlags:           ' + str(size(flatnonzero(KUFlag))))
            print('        nr. of SKFlags:           ' + str(size(flatnonzero(SKFlag))))
            print('        nr. of MoSFlags:          ' + str(size(flatnonzero(MoSFlag))))
        else:
            print('        No RFI analysis in the time domain possible (nps < 10). All TD flags are set to -99.')
        # RFI Analysis in the frequency domain
        print('     RFI analysis in the frequency domain')
        print('        FD:                       %.2f K - %.2f K'%(amin(FD), amax(FD)))
        print('        nr. of FDFlags:           ' + str(size(flatnonzero(FDFlag))) )
        # System temperature check
        print('     system temperature check')
        print('        T_set - T0:               %.2f K - %.2f K'%(amin(SysTempDiff), amax(SysTempDiff)))
        print('        nr. of SysTempFlags:      ' + str(size(flatnonzero(SysTempFlag)))  )
        # ACS temperature plausibility check
        print('     ACS temperature plausibility check')
        print('        nr. of ACSTempFlags:      ' + str(size(flatnonzero(ACSTempFlag))) )
        # values written to time series to be used for raw data processing
        print('     ACS temperature time series')
        print('        nr. of unflagged samples: ' + str(size(i)))
     
        print('        results written to ACS temperature time series: %02d.%02d.%4d %02d:%02d:%02d, T_ACS = %.4f K, Flag = %d'%
        (NextRow[0].day, NextRow[0].month, NextRow[0].year, NextRow[0].hour, NextRow[0].minute, NextRow[0].second, NextRow[1], NextRow[2]))




# 12.) Save 'TACSOut', 'TempOut' and 'VoltOut' to file
# ----------------------------------------------------
print('')
print('***************************************************************')
print('')

# 12.1) save 'TACSOut', 'TempOut' and 'VoltOut' to numpy binary file
print('Saving results to ' + OutFileName + 'TACS.npy, ' + OutFileName + 'Temp.npy, and ' + OutFileName + 'Volt.npy.')
save(OutFileName + 'TACS', TACSOut)
save(OutFileName + 'Temp', TempOut)
save(OutFileName + 'Volt', VoltOut)

# 12.2) save 'TACSOut', 'TempOut' and 'VoltOut' to text file
print('Saving results to ' + OutFileName + 'TACS.txt, ' + OutFileName + 'Temp.txt, and ' + OutFileName + 'Volt.txt.')
WriteResults2Text(OutFileName, TempOut, VoltOut, TACSOut)


# 13.) Print plot time series of results ('TACSOut' and 'TempOut')
# ----------------------------------------------------------------
PlotTimeSeries(TempOut, TACSOut)
