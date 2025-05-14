#!c:/Python25/python.exe
#-*-coding: utf-8-*-

from numpy import *
from scipy.stats import skew
from scipy.stats import kurtosis
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
import librad
import sys


def CombineFlags(Flag1, Flag2, Flag3, Flag4, Flag5, Flag6):
    """
    Combines the given flags 'Flag1' to 'Flag6' into one unambigous value.

    usage: Flag = CombineFlags(Flag1, Flag2, Flag3, Flag4, Flag5, Flag6)
    inputs:
        Flag1  array of Flag1
        Flag2  array of Flag2
        ...
    output:
        Flag  array of corresponding combined flags

    author: Ingo Voelksch, 6.9.2011
    """

    # compute Flag
    if (Flag1[0] == -99):
        Flag = -(999000 + (100*Flag4) + (10*Flag5) + (1*Flag6))
    else:
        Flag = (100000*Flag1) + (10000*Flag2) + (1000*Flag3) + (100*Flag4) + (10*Flag5) + (1*Flag6)

    # return Flag to invoking function
    return Flag


# **************************************************************************** 


def SplitFlags(Flags):
    """
    Splits the six-digit flag 'Flag' into its individual components 'Flag1' to 'Flag6'.

    usage: Flag1, Flag2, Flag3, Flag4, Flag5, Flag6 = SplitFlags(Flag)
    input:
        Flag   array of six-digit flag
    outputs:
        Flag1  array of corresponding Flag1
        Flag2  array of corresponding Flag2
        ...

    author: Ingo Voelksch, 31.10.2011
    """

    # This is a quite crude approach to splitting 'Flags' into its individual components
    # 'Flag1' to 'Flag6', since for some reason, I could not make 'rint' work on the quick.
    # Anyway, it works like this also and so I leave it at that for the moment.

    # add 0.1 to Flags to avoid problems when converting floats to int below
    Flags = Flags+0.1

    # split 'Flags' into individual values 'Flag1' to 'Flag6'
    # extract 'Flag1'
    F  = (Flags/100000)
    Flag1 = F*1
    for i in range(0, len(Flag1)):   # loose all decimals by converting to int
       Flag1[i] = int(Flag1[i])

    # extract 'Flag2'
    F = (F-Flag1)*10
    Flag2 = F*1
    for i in range(0, len(Flag2)):
       Flag2[i]=int(Flag2[i])

    # extract 'Flag3'
    F = (F-Flag2)*10
    Flag3 = F*1
    for i in range(0, len(Flag3)):
       Flag3[i]=int(Flag3[i])

    # extract 'Flag4'
    F = (F-Flag3)*10
    Flag4 = F*1
    for i in range(0, len(Flag4)):
       Flag4[i]=int(Flag4[i])

    # extract 'Flag5'
    F = (F-Flag4)*10
    Flag5 = F*1
    for i in range(0, len(Flag5)):
       Flag5[i]=int(Flag5[i])

    # extract 'Flag6'
    F = (F-Flag5)*10
    Flag6 = F*1
    for i in range(0, len(Flag6)):
       Flag6[i]=int(Flag6[i])

    # return Flag1 to Flag6 to invoking function
    return Flag1, Flag2, Flag3, Flag4, Flag5, Flag6


# **************************************************************************** 


def CheckValuePlausibility(minValue, maxValue, Value):
    """
    Determines, whether a given value lies between minValue and maxValue and
    sets the corresponding flags. This function is intended for:
    i.)   checking the plausibility of the derived ACS temperature
    ii.)  checking the plausibility of the derived brightness temperature
    iii.) checking the mean over standard deviation, skewness, and kurtosis
          criterion for  RFI analysis in the time domain

    usage: Flag = CheckValuePlausibility(minValue, maxValue, Value)
    inputs:
        minValue  array or scalar of minimum allowed value
        maxValue  array or scalar of maximum allowed value
        Value     array of value to be checked
    output:
        Flag  array of corresponding flag
              Flag = 1: Value < minValue or Value > maxValue
              Flag = 0: minValue <= Value <= maxValue

    author: Ingo Voelksch, 5.9.2011
    """

    # initialize array Flag and set to 0 for minValue <= Value <= maxValue
    Flag = ones( (size(Value,0),1) )
    if len(Flag[0]) == 1:
      Flag[all([Value>=minValue, Value<=maxValue], axis=0)] = 0
    else:
      Flag[all([Value>=minValue, Value<=maxValue], axis=0),:] = 0

    # return Flag to invoking function
    return Flag


# **************************************************************************** 


def U2TB_in(T_HS, T_CS, U_HS, U_CS, U):
    """
    Calculates the brightness temperature at the radiometer input port from
    the measured radiometer raw data (physical temperatures and voltages
    measured for the cold and hot calibration target, voltages measured for
    the scene).
    
    usage: TB_in = U2TB_in(T_HS, T_CS, U_HS, U_CS, U)
    inputs:
        T_HS  array of physical temperature [K] of the hot calibraton target
        T_CS  array of physical temperature [K] of the cold calibration target
        U_HS  array of voltages [V] measured for the hot calibration target
        U_CS  array of voltages [V] measured for the cold calibration target
        U     array of voltages [V] measured for the scene
    outputs:
        TB_in  array of scene brightness temperatures [K] at radiometer input
               port

    author:  Ingo Voelksch, 2.9.2011
    """

    TB_in = ( (T_HS - T_CS) / (U_HS - U_CS) ) * (U - U_CS) + T_CS

    # return TB_in to invoking function
    return TB_in


# **************************************************************************** 


def StatsOverRecord(DataRecord):
    """
    Caculates the mean value, standard deviation, mean over standard deviation,
    skewness, and kurtosis of each row in the array 'DataRecord'.
    
    usage: MeanValue, StandardDeviation, MeanOverSTD, Skewness, Kurtosis = 
           StatsOverRecord(DataRecord)
    inputs:
        DataRecord:  array of voltages [V] measured for one noise source
    outputs:
        MeanValue          array of mean value [V]
        StandardDeviation  array of standard deviation [V]
        MeanOverSTD        array of (MeanValue/StandardDeviation)
        Skewness           array of skewness
        Kurtosis           array of kurtosis (using Fisher's definition)

    author: Ingo Voelksch, 31.8.2011
    """

    # initialize arrays and set to -99
    MeanValue         = ones( (size(DataRecord,0),1) ) * -99
    StandardDeviation = ones( (size(DataRecord,0),1) ) * -99
    MeanOverSTD       = ones( (size(DataRecord,0),1) ) * -99
    Skewness          = ones( (size(DataRecord,0),1) ) * -99
    Kurtosis          = ones( (size(DataRecord,0),1) ) * -99

    # calculate statistics only if at least 3 values are found, otherwise
    # values of -99 are kept (except for mean value)
    MeanValue[:,0] = mean(DataRecord, axis=1)
    if ( size(DataRecord, 1) >= 3):
        StandardDeviation[:,0] = std(DataRecord, axis=1, ddof=0)  # ddof=0 => std is calculated with number of elements as divisor
        MeanOverSTD[:,0]       = MeanValue[:,0] / StandardDeviation[:,0]
        Skewness[:,0]          = skew(DataRecord, axis=1)
        Kurtosis[:,0]          = kurtosis(DataRecord, axis=1, fisher=1)  # fisher=1 => 3.0 is substracted from kurtosis to give 0.0 for normal distribution

    # return values to invoking function
    return MeanValue, StandardDeviation, MeanOverSTD, Skewness, Kurtosis


# ****************************************************************************


def ReadDataFile(File2Read):
    """
    Opens an ELBARAII measurement data file and extracts the values needed
    for further processing.
    
    usage: LoadsMeasured, n_Loads, n_Samples, DatesMeasured, T_plate, T_ext,
           ElevationAngle, nps, Voltages = ReadDataFile(File2Read)
    inputs:
        File2Read  name of file to read
    outputs:
        LoadsMeasured   list of string with identifiers of noise sources measured
                        (HL, CL, AL, HP, VP)
        n_Loads         number of noise sources measured
        n_Samples       number of samples (rows) in the data file
        DatesMeasured   list of measurement dates and times as datetime.datetime
                        objects (year, month, day, hour, minute, second)
        T_plate         array of internal system temperatures [K]
        T_ext           array of external air temperatures [K]
        ElevationAngle  array of elevation angles [deg.]
        nps             number of individual values per noise source in data file
        Voltages        array of voltages [V] (all noise sources and samples)
    
    author: Ingo Voelksch, 30.8.2011
    """

    # open data file
    # --------------
    DataFile  = open(File2Read, 'rb')
    
    # read file header and derive parameters needed for further processing
    # --------------------------------------------------------------------
    adc_samp_freq, adc_num_samp, adc_presum, num_cycles, cycle_delay, start_elev_ang, elev_step, elev_num_step, measure_modes, hdr_len, rec_len = librad.Read_FHeader(DataFile)
    # Wenn nötig, könnte ich anhand dieser Sachen dann noch jede Menge Info zur
    # gerade prozessiertem Datei ausgeben. Ich habe jetzt mal alles, was ich nicht
    # gerade brauche, weggeschmissen.
    LoadsMeasured = measure_modes
    n_Loads       = len(LoadsMeasured)
    nps           = adc_num_samp / adc_presum
    n_Samples     = num_cycles*elev_num_step

    # read data from file and split into individual variables
    # -------------------------------------------------------
    # initialize arrays
    DatesMeasured  = []
    ElevationAngle = zeros( (n_Samples, 1) )
    T_plate        = zeros( (n_Samples, 1) )
    T_ext          = zeros( (n_Samples, 1) )
    Voltages       = zeros( (n_Samples, 2*nps*len(LoadsMeasured)) )

    # read data from file
    i = 0
    DataFile.seek(librad.HDR_SIZE, 0)    # move to start of the data in data file
    while 1:  # read data line by line
        try:
           # read measurement date, physical temperatures, etc.
           year, month, day, sod, elev_ang, t_plate, t_sink, t_ext, t_cable, dc_offset0, dc_offset1, nps = librad.Read_LHeader(DataFile)
           HH, MM, SS = librad.Sod2hms(sod)    # convert seconds of day to hours, minutes, seconds
        except:  # quit loop, when end of file is reached
            break
        # read voltages for all noise sources
        volts = fromfile(file=DataFile, dtype=float32, count=2*nps*len(LoadsMeasured))
        # assign values to corresponding arrays
        DatesMeasured.append( datetime(year, month, day, int(HH), int(MM), int(SS)) )
        ElevationAngle[i, :]    = elev_ang
        T_plate[i, :]           = t_plate + 273.15  # convert to [K]
        T_ext[i, :]             = t_ext   + 273.15  # convert to [K]
        Voltages[i,:]           = volts
        i = i + 1
        
    # close data file and return values to invoking function
    DataFile.close()
    return LoadsMeasured, n_Loads, n_Samples, DatesMeasured, T_plate, T_ext, ElevationAngle, nps, Voltages


# ****************************************************************************


def ReadParameterFile(FileName, MeasureMode):
    """
    Read file with evaluation parameters necessary for ACSCalibration.py and
    RawDataProcessing.py and return dictionary with corresponding values.

    usage: EvalPar = ReadParameterFile(FileName, MeasureMode)
    inputs:
        FileName     name of file with evaluation parameters
        MeasureMode  measurement mode
    outputs:
        EvalPar  dictionary with evaluation parameter values

    author: Ingo Voelksch, 21..2011
            Andreas Wiesmann, 20160524 transmission loss added
    """

    # define dictionary
    EvalPar = {
        'trans_loss'   : '',
        't_set'        : '',
        'dir2evaluate' : '',
        'outfilename'  : '',
        'mosnorfi'     : '',
        'mosdiffmax'   : '',
        'skdiffmax'    : '',
        'kudiffmax'    : '',
        'fdmax'        : '',
        'tempdiffmax'  : '',
        }
    if MeasureMode == 'sky':  # add specific parameters for ACS calibration
        EvalPar.update({'altitude':'', 't_acsmin':'', 't_acsmax':''})
    elif MeasureMode == 'scene':  # add specific parameters for raw data processing
        EvalPar.update({'acstempfile':'', 'tb_hpmin':'', 'tb_hpmax':'', 'tb_vpmin':'', 'tb_vpmax':''})
    else:
        print('Unknown measurement mode. Abort.')
        sys.exit(-1)
    
    pf     = open(FileName, 'r')  # open parameter file for reading
    ParSet = 0  # set number of assigned parameter values to 0

    # read parameter file and extract parameter values
    NextLine = pf.readline()
    LineNr   = 1  # set line counter to 1
    while (NextLine):
        # skip this line, if 'NextLine' is empty or starts with '#'
        if (NextLine[0] == '#') or (len(NextLine.strip())) == 0:
            NextLine = pf.readline()
            LineNr   = LineNr + 1
            continue

        # split 'NextLine' into its individual parts, which are separated by a colon
        kw_val = NextLine.split(':')

        # ignore this line, if no or more than 2 colons are found
        if (len(kw_val) == 1) or (len(kw_val) > 3):
            print('   Don\'t understand line ' + str(LineNr) + '. Am just ignoring this line.')
            NextLine = pf.readline()
            LineNr   = LineNr + 1
            continue

        # extract keyword 'kw' and corresponding parameter value 'val' from 'kw_val'
        kw  = kw_val[0].strip(); kw = kw.lower()
        val = kw_val[1].strip()
        # If three entries are in 'kw_val', this is probably due to a directory
        # name in Windows format => Concatenate second and third entry, separated
        # by a colon.
        if (len(kw_val) == 3):
            val = val + ':' + kw_val[2].strip()

        # ignore this line, if no entry besides an end of line is found in 'val'
        if len(val) == 0:
            print('   Something seems to be missing in line ' + str(LineNr) + '. Am just ignoring this line.')
            NextLine = pf.readline()
            LineNr   = LineNr + 1
            continue

        # assign parameter value to corresponding keyword in dictionary
        if kw in EvalPar:
            if (kw == 'dir2evaluate') or (kw == 'outfilename') or (kw == 'acstempfile'):
                EvalPar[kw] = val
            else:  # convert to float for numerical values
                EvalPar[kw] = float(val)
            ParSet = ParSet + 1
        # if keyword is unknown, ignore this line
        else:
            print('   Don\'t understand keyword in line ' + str(LineNr) + '. Am just ignoring this line.')

        # read next line and increase line counter
        NextLine = pf.readline()
        LineNr   = LineNr + 1

    # close parameter file
    pf.close()

    # test if all parameters are set
    if MeasureMode == 'sky':
        if ParSet < 12:
            print('Not enough evaluation parameters for ACSCalibration found! Abort program.')
            sys.exit(-1)
    if MeasureMode == 'scene':
        if ParSet < 14:
            print('Not enough evaluation parameters for RawDataProcessing found! Abort program.')
            sys.exit(-1)

    # return dictionary 'EvalPar' to invoking function
    return(EvalPar)


# ****************************************************************************


def TSkyPellarin(Z, T_2m, theta):
    """
    Calculates the theoretical sky brightness after Pellarin et al. (2003).

    usage: T_sky_Pellarin = TSkyPellarin(Z, T_2m, theta)
    inputs:
        Z      altitude in [km]
        T_2m   array of air temperature 2 m above the ground in [K]
        theta  array of zenith angle in [deg.]
    outputs:
        T_sky_Pellarin  array of theoretical sky brightness in [K]

    author: Ingo Voelksch, 29.8.2011
    """

    # convert theta from degree to radian
    theta = deg2rad(theta)
  
    # compute sky brightness
    tau_ATM        = exp( -3.9262 - 0.2211*Z - 0.00369*T_2m )    # atmosphere optical thickness
    T_ATMeq        = exp( 4.9274 + 0.002195*T_2m )               # atmosphere equivalent temperature
    TB_COSD        = 2.7 * exp( -tau_ATM/cos(theta) )            # downwelling cosmic emission
    TB_ATMD        = T_ATMeq * (1 - exp( -tau_ATM/cos(theta)) )  # downward atmospheric brightness
    T_sky_Pellarin = TB_ATMD + TB_COSD                           # entire sky radiation received

    # return T_sky_Pellarin to invoking function
    return(T_sky_Pellarin)


# ****************************************************************************


def PlausibilityChecks(ElevationAngle, mode):
    """
    Performs some plausibility checks for the data within one measurement file
    and aborts the invoking function, if inplausibilities are encountered.

    usage: PlausibilityChecks(ElevationAngle, mode)
    inputs:
        ElevationAngle         array of elevation angle [deg.]
        mode                   measurement mode ("sky" = sky calibration,
                                                 "scene" = normal measurement)
    outputs:
        none
    
    author: Ingo Voelksch, 31.8.2011
    """

    # checks for sky measurements
    if (mode == "sky"):
        # check, if the elevation angle changes within one measurement file
        if (amin(ElevationAngle) != amax(ElevationAngle)):
            print('    Warning: Different elevation angles found within the data file!')
        # check, if elevation angle is sensible for sky measurements
        if ( (amin(ElevationAngle) < 95) or (amax(ElevationAngle) > 265) ):
            print('    Warning: Elevation angles < 95 deg found! You might see more than just sky with your antenna.')

    # checks for scene measurements
    elif (mode == "scene"):
        # check if the elevation angle is meaningful
        if (amin(ElevationAngle) < 0) or (amax(ElevationAngle) > 360):
            print("Meaningless elevation angle outside 0 deg <= alpha <= 360 deg found. Abort.")
            sys.exit(-1)

    else:
        print("Unrecognized measurement mode. Abort program.")
        sys.exit(-1)


# ****************************************************************************


def TB2TB_in(TB, T_FC, L_FC):
    """
    Calculates the noise power added to the brightness temperature entering
    the antenna aperture due to the lossy feed cables, and then calculates
    the actual brightness temperature received at the radiometer input port.

    usage: TB_in = TB2TB_in(TB, T_FC, L_FC)
    inputs:
        TB    array of brightness temperature [K] entering the antenna
        T_FC  array of feed-cable temperature [K]
        L_FC  transmission loss of the feed cable [dB]
    outputs:
        TB_in  array of brightness temperature [K] at radiometer input port

    author: Ingo Voelksch, 30.8.2011
    """

    t_FC    = power(10, (-L_FC/10))     # transmission factor of feed cable
    DeltaTB = (1 - t_FC) * (T_FC - TB)  # noise added due to lossy feed cable
    TB_in   = TB + DeltaTB              # noise temperature at input port
    # return TB_in to invoking function
    return(TB_in)


# ****************************************************************************


def TB_in2TB(TB_in, T_FC, L_FC):
    """
    Calculates the actual brightness temperature entering the antenna
    aperture from the noise temperature received at the radiometer input
    port (i.e., the cable loss is substracted).

    usage: TB = TB_in2TB(TB_in, T_FC, L_FC)
    inputs:
        TB_in  array of noise temperature at radiometer input port [K]
        T_FC   array of temperature of the feed cable [K]
        L_FC   transmission loss of the feed cable [dB]
    outputs:
        TB  array of brightness temperature entering antenna aperture [K]

    author: Ingo Voelksch, 30.8.2011
    """

    t_FC = power(10, (-L_FC/10))           # transmission factor of feed cable
    TB   = (TB_in - (1-t_FC)*T_FC) / t_FC  # noise temperature entering antenna
    # return TB to invoking function
    return(TB)


# ****************************************************************************


def CheckTempDiff(T1, T2, DeltaTMax):
    """
    Calculates the temperature difference between T1 and T2 and the corres-
    ponding flags. This script is intended for:
    i.)  checking the deviation of T_set from T_0 (T1 = T_set, T2 = T_0), and
    ii.) the evaluation of the channel difference (T1 = T_B1, T2 = T_B2).

    usage: DeltaT, Flag = CheckTempDiff(T1, T2, DeltaTMax)
    inputs:
        T1         array or scalar of temperature 1 [K]
        T2         array of temperature 2 [K]
        DeltaTMax  maximum allowed temperature difference [K]
    outputs:
        DeltaT  array of temperature difference [K]
        Flag    array of corresponding flag:
                Flag = 1: abs(DeltaT) >  DeltaTMax
                Flag = 0: abs(DeltaT) <= DeltaTMax

    author: Ingo Voelksch, 6.9.2011
    """

    # compute temperature difference
    DeltaT = T1 - T2

    # initialize array 'Flag' and set to 1 for abs(DeltaT)>DeltaTMax
    Flag                               = zeros( (DeltaT.shape[0],1) )
    Flag[absolute(DeltaT) > DeltaTMax] = 1

    # return DeltaT and Flag to invoking function
    return DeltaT, Flag


# ****************************************************************************


def WriteResults2Text(FileName, TempData, VoltData, TACSTimeSeries = 0):
    """
    Writes the results of 'ACSCalibration.py' and 'RawDataProcessing.py' to
    text files.

    usage: WriteResults2Text(FileName, TempData, VoltData, TACSTimeSeries)
    inputs:
        FileName        name of file to write
        TempData        array with brightnes temperatures etc. ('TempOut')
        VoltData        array with voltage statistics data ('VoltOut')
        TACSTimeSeries  array with time series of ACS temperatures to be used
                        as input for 'RawDataProcessing.py'
                        This argument is optional

    author: Ingo Voelksch, 16.9.2011 
    """

    # write TempData
    fid1 = open(FileName + 'Temp.txt', 'w')
    fid1.write('Date; Elevation Angle; TB_HP; FD_HP; Flag_HP; TB_VP; FD_VP; Flag_VP; T_plate; SysTempDiff; T_ext\n')
    for line in range(0, size(TempData, 0)):
        fid1.write(TempData[line,0].strftime('%d-%m-%Y %H:%M:%S') + '; %.1f; %.4f; %.4f; %06d; %.4f; %.4f; %06d; %.3f; %.3f; %.3f'%
          (TempData[line,1], TempData[line,2], TempData[line,3], TempData[line,4], TempData[line,5], TempData[line,6], TempData[line,7], TempData[line,8], TempData[line,9], TempData[line,10]))
        fid1.write('\n')
    fid1.close()
    
    # write  VoltData
    fid2 = open(FileName + 'Volt.txt', 'w')
    fid2.write(
    'Date; U_ACS1_Mean; U_ACS1_STD; U_ACS1_MoS; U_ACS1_SK; U_ACS1_KU; U_ACS2_Mean; U_ACS2_STD; U_ACS2_MoS; U_ACS2_SK; U_ACS2_KU; U_RS1_Mean; U_RS1_STD; U_RS1_MoS; U_RS1_SK; U_RS1_KU; U_RS2_Mean; U_RS2_STD; U_RS2_MoS; U_RS2_SK; U_RS2_KU; U_HP1_Mean; U_HP1_STD; U_HP1_MoS; U_HP1_SK; U_HP1_KU; U_HP2_Mean; U_HP2_STD; U_HP2_MoS; U_HP2_SK; U_HP2_KU; U_VP1_Mean; U_VP1_STD; U_VP1_MoS; U_VP1_SK; U_VP1_KU; U_VP2_Mean; U_VP2_STD; U_VP2_MoS; U_VP2_SK; U_VP2_KU\n'
    )
    for line in range(0, size(VoltData, 0)):
        fid2.write(VoltData[line,0].strftime('%d-%m-%Y %H:%M:%S')
          + '; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f; %.4f'%
          (VoltData[line,1],  VoltData[line,2],  VoltData[line,3],  VoltData[line,4],  VoltData[line,5],  VoltData[line,6],  VoltData[line,7],  VoltData[line,8],
           VoltData[line,9],  VoltData[line,10], VoltData[line,11], VoltData[line,12], VoltData[line,13], VoltData[line,14], VoltData[line,15], VoltData[line,16],
           VoltData[line,17], VoltData[line,18], VoltData[line,19], VoltData[line,20], VoltData[line,21], VoltData[line,22], VoltData[line,23], VoltData[line,24],
           VoltData[line,25], VoltData[line,26], VoltData[line,27], VoltData[line,28], VoltData[line,29], VoltData[line,30], VoltData[line,31], VoltData[line,32],
           VoltData[line,33], VoltData[line,34], VoltData[line,35], VoltData[line,36], VoltData[line,37], VoltData[line,38], VoltData[line,39], VoltData[line,40]))
        fid2.write('\n')
    fid2.close()

    # write ACS temperature time series, if corresponding data is given
    if ( type(TACSTimeSeries) == ndarray ):
        fid3 = open(FileName + 'TACS.txt', 'w')
        fid3.write('Date; T_ACS; T_ACS_HP1; T_ACS_HP2; T_ACS_VP1; T_ACS_VP2; Flag\n')
        for line in range(0, size(TACSTimeSeries, 0)):
            fid3.write(TACSTimeSeries[line,0].strftime('%d-%m-%Y %H:%M:%S') + '; %.4f; %.4f; %.4f; %.4f; %.4f; %d'%(TACSTimeSeries[line,1], TACSTimeSeries[line,2], TACSTimeSeries[line,3], TACSTimeSeries[line,4], TACSTimeSeries[line,5], TACSTimeSeries[line,6]))
            fid3.write('\n')
        fid3.close()


# ****************************************************************************


def DatetimeObject2DatetimeArray(Datetime_dto):
    """
    Converts a list of datetime objects into a corresponding array.

    usage: Datetime_arr = DatetimeObject2DatetimeArray(Datetime_dto)
    input:
        Datetime_dto  list with datetime objects
    output:
        Datetime_arr  array of corresponding dates and times: [year, month,
                      day, hour, minute, second]

    author:  Ingo Voelksch, 8.9.2009
    """

    # initialize array for output
    Datetime_arr = ones( (len(Datetime_dto),6) ) * -99
    # fill columns and rows of 'Datetime_arr)
    for i in range(0, len(Datetime_dto)):
        Datetime_arr[i,0] = Datetime_dto[i].year
        Datetime_arr[i,1] = Datetime_dto[i].month
        Datetime_arr[i,2] = Datetime_dto[i].day
        Datetime_arr[i,3] = Datetime_dto[i].hour
        Datetime_arr[i,4] = Datetime_dto[i].minute
        Datetime_arr[i,5] = Datetime_dto[i].second
    # return 'Datetime_arr' to invoking function
    return Datetime_arr

    
# ****************************************************************************


def DatetimeObject2DatetimeNum(Datetime_dto):
    """
    Converts a list of datetime objects into a list with numerical datetimes.

    usage: Datetime_num = DatetimeObject2DatetimeNum(Datetime_dto)
    input:
        Datetime_dto  list with datetime objects
    output:
        Datetime_num  list with numerical datetime

    author: Ingo Voelksch, 8.9.2011
    """

    # initialize 'Datetime_num'
    Datetime_num = []
    # convert datetimes
    for i in range(0, len(Datetime_dto)):
        date_num = Datetime_dto[i].toordinal()  # date (year, month, day) as number
        time_num = (float(Datetime_dto[i].hour)/24) + (float(Datetime_dto[i].minute)/1440) + (float(Datetime_dto[i].second)/86400)  # time (hour, minute, second) as number
        Datetime_num.append( (float(date_num) + time_num) )
    # return 'Datetime_num' to invoking function
    return Datetime_num

    
# ****************************************************************************


def DatetimeNum2DatetimeObject(Datetime_num):
    """
    Converts a list with numerical datetimes into a list of datetime objects.

    usage: Datetime_dto = DatetimeNum2DatetimeObject(Datetime_num)
    input:
        Datetime_num  list with numerical datetime
    output:
        Datetime_dto  list with datetime objects

    author: Ingo Voelksch, 8.9.2011
    """

    # initialize date => 'Datetime_dto'
    Datetime_dto = []
    # convert datetimes
    for i in range(0, len(Datetime_num)):
        # extract numerical date and convert to datetime object
        date_num = int(Datetime_num[i])
        date_num = datetime.fromordinal(date_num)
        # extract numerical time and convert to hour, minute, second
        time_num = Datetime_num[i] - int(Datetime_num[i])
        hour     = int(time_num*24)
        minute   = int((time_num - (float(hour)/24)) * 1440)
        second   = int((time_num - (float(hour)/24) - (float(minute)/1440)) * 86400)
        # replace hour, minute, and second in 'datenum' with previously calculated values
        date_num = date_num.replace(date_num.year, date_num.month, date_num.day, hour, minute, second)
        Datetime_dto.append(date_num)
    # return 'Datetime_dto' to invoking function
    return Datetime_dto


# ****************************************************************************


def PlotTimeSeries(Temp, TACS=-1):
    """
    Plots the results (ACS / brightness temperatures and corresponding flags)
    of 'ACSCalibration.py' and 'RawDataProcessing.py'

    usage: PlotTimeSeries(Temp, TACS)
    inputs:
        Temp  array 'TempOut' from 'RawDataProcessing.py' or 'ACSCalibration.py'
        TACS  arry 'TACSOut' from 'ACSCalibration.py' (optional argument)

    author: Ingo Voelksch, 31.10.2011
            Andreas Wiesmann, 20160524 plot beautification
    """

    # adjust matploblibrc file to change some default properties for all figures
    mpl.rc('legend', fancybox=True, numpoints=1, fontsize=10, borderpad=0.25)   # use fancy legend, show 1 symbol only, fontsize of legend entries, border whitespace
    mpl.rc('axes', grid=True, labelsize=12, titlesize=13)                       # display grid, fontsize of axis titles, fontsize of subplot title
    mpl.rc('xtick', labelsize=11)                                               # fontsize of x-tick labels
    mpl.rc('ytick', labelsize=11)                                               # fontsize of y-tick labels
    mpl.rc('figure.subplot', hspace=0)                                          # set vertical distance between subplots to 0


    # Plots for results of ACSCalibration
    if ( type(TACS) == ndarray ):
        # 1.) plot T_ACS time series for RawDataProcessing into Fig1
        # ----------------------------------------------------------
        # 1.1) find indices of unflagged / flagged measurements
        Unflagged = flatnonzero(TACS[:,6]==0)
        Flagged   = flatnonzero(TACS[:,6])
    
        # 1.2) plot T_ACS data into Fig1
        Fig1 = plt.figure(1, figsize=(8,8))
        plt.plot(TACS[Unflagged,0], TACS[Unflagged,2], 'k.', label=r'$T_{\mathrm{\mathrm{ACS\_HP1}}}$ unflagged')
        plt.plot(TACS[Flagged,0],   TACS[Flagged,2],   'r.', label=r'$T_{\mathrm{\mathrm{ACS\_HP1}}}$ flagged')
        plt.plot(TACS[Unflagged,0], TACS[Unflagged,3], 'g.', label=r'$T_{\mathrm{\mathrm{ACS\_HP2}}}$ unflagged')
        plt.plot(TACS[Flagged,0],   TACS[Flagged,3],   'b.', label=r'$T_{\mathrm{\mathrm{ACS\_HP2}}}$ flagged')
        plt.plot(TACS[Unflagged,0], TACS[Unflagged,4], 'c.', label=r'$T_{\mathrm{\mathrm{ACS\_VP1}}}$ unflagged')
        plt.plot(TACS[Flagged,0],   TACS[Flagged,4],   'm.', label=r'$T_{\mathrm{\mathrm{ACS\_VP1}}}$ flagged')
        plt.plot(TACS[Unflagged,0], TACS[Unflagged,5], 'b.', label=r'$T_{\mathrm{\mathrm{ACS\_VP2}}}$ unflagged')
        plt.plot(TACS[Flagged,0],   TACS[Flagged,5],   'y.', label=r'$T_{\mathrm{\mathrm{ACS\_VP2}}}$ flagged')

        # 1.3) format plots
        plt.ylabel(r'$T_{\mathrm{ACS}}$ [K]')   # title of y-axis
        plt.xlabel('date of sky measurement')   # title of x-axis
        Fig1.autofmt_xdate()                    # format x-labels, so that dates are readable
        plt.legend()                            # show legend using labels given in plot command
        plt.title(r'Time series $T_{\mathrm{\mathrm{ACS}}}(time)$ as input for $\mathtt{RawDataProcessing}$')   # title of plot
        plt.subplots_adjust(top=0.88, bottom = 0.15)   # adjust top and bottom position of plots in figure

        # 2.) plot T_ACS of all data samples and corresponding flags into Fig2
        # --------------------------------------------------------------------
        # 2.1) find indices of unflagged / flagged measurements
        # first, find all samples (if any) where no RFI analysis in the TD was carried out (i.e., flags start with -999) and
        # remove the TD flags
        NoTDFlags = flatnonzero(Temp[:,4]<0)
        Temp[NoTDFlags,4] = (Temp[NoTDFlags,4]*-1) - 999000
        # then, test which flags are 0
        Unflagged = flatnonzero(Temp[:,4]==0)
        Flagged   = flatnonzero(Temp[:,4])

        # 2.2) create Fig2 for plotting T_ACS and flags
        Fig2 = plt.figure(2, figsize=(8,10))
    
        # 2.3) plot T_ACS into subplot 'PlotTACS'
        PlotTACS = plt.subplot2grid((4,1),(0,0), rowspan=3)   # subplot for T_ACS, spanning 3 rows
        plt.plot(Temp[Unflagged,0], Temp[Unflagged,2], 'k.', label=r'$T_{\mathrm{\mathrm{ACS\_1}}}$ unflagged')
        plt.plot(Temp[Flagged,0],   Temp[Flagged,2],   'r.', label=r'$T_{\mathrm{\mathrm{ACS\_1}}}$ flagged')
        plt.plot(Temp[Unflagged,0], Temp[Unflagged,5], 'g.', label=r'$T_{\mathrm{\mathrm{ACS\_2}}}$ unflagged')
        plt.plot(Temp[Flagged,0],   Temp[Flagged,5],   'b.', label=r'$T_{\mathrm{\mathrm{ACS\_2}}}$ flagged')

        # 2.4) format 'PlotTACS'
        plt.ylabel(r'$T_{\mathrm{ACS}}$ [K]')
        plt.legend()
        plt.title(r'Calibrated $T_{\mathrm{\mathrm{ACS}}}$ of all samples and corresponding flags for distortion')

        # 2.5) plot flags into 'PlotFlags'
        PlotFlags = plt.subplot2grid((4,1),(3,0), rowspan=1, sharex = PlotTACS)   # subplot for flags, spanning 1 row, shared x-axis with 'PlotTACS'
        # first, derive all individual flags with 'SplitFlags'
        F1, F2, F3, F4, F5, F6 = SplitFlags(Temp[:,4])
        # then, plot individual flags into 'PlotFlags'
        plt.plot(Temp[:,0], F6*5,   linestyle='none', marker='x', markeredgecolor='m',      markerfacecolor='w', label=r'$T_{\mathrm{inst}}$ flag')   # F6
        plt.plot(Temp[:,0], F5*4,   linestyle='none', marker='x', markeredgecolor='orange', markerfacecolor='w', label=r'$T_{\mathrm{ACS}}$ flag')    # F5
        plt.plot(Temp[:,0], F4*3,   linestyle='none', marker='o', markeredgecolor='k',      markerfacecolor='w', label=r'$FD$ flag')                  # F4
        plt.plot(Temp[:,0], F3*2,   linestyle='none', marker='o', markeredgecolor='r',      markerfacecolor='w', label=r'$w_{\mathbf{U}}$ flag')      # F3
        plt.plot(Temp[:,0], F2*1.5, linestyle='none', marker='o', markeredgecolor='g',      markerfacecolor='w', label=r'$s_{\mathbf{U}}$ flag')      # F2
        plt.plot(Temp[:,0], F1,     linestyle='none', marker='o', markeredgecolor='b',      markerfacecolor='w', label=r'$k_{\mathbf{U}}$ flag')      # F1

        # 2.5) format 'PlotFlags'
        PlotFlags.set_yticks((1, 1.5, 2, 3, 4, 5))   # set position of y-tick labels
        plt.ylim(0.5, 5.5)   # set y-axis limits
        labels = PlotFlags.set_yticklabels((r'$k_{\mathbf{U}}$', r'$s_{\mathbf{U}}$', r'$w_{\mathbf{U}}$', r'$FD$', r'$T_{\mathrm{ACS}}$', r'$T_{\mathrm{inst}}$'))   # define y-labels
        plt.xlabel('date of sky measurement')
        plt.ylabel('distortion check')
        plt.legend()
        Fig2.autofmt_xdate()   # format x-labels, so that dates are readable
        plt.subplots_adjust(bottom = 0.13)   # adjust top and bottom position of plots in figure

        plt.show()
        plt.close('all')

    # plots for results of RawDataProcessing
    else:
        # 1.) plot TB's of all data samples and corresponding flags into Fig1
        # ------------------------------------------------------------------------
        # 1.1) create Fig1 for plotting
        Fig1 = plt.figure(1, figsize=(15.76, 10))
        Fig1.suptitle('Calibrated brightness temperatures of all samples and corresponding flags for distortion', fontsize=13)
    
        # 1.2) plot TB_H's into upper subplot of left column of Fig1
        # find indices of unflagged / flagged measurements:
        #    first, find all samples (if any) where no RFI analysis in the TD was carried out (i.e., flags start with -999) and remove the TD flags
        NoTDFlags = flatnonzero(Temp[:,4]<0)
        Temp[NoTDFlags,4] = (Temp[NoTDFlags,4]*-1) - 999000
        #    then, test which flags are 0
        Unflagged = flatnonzero(Temp[:,4]==0)
        Flagged   = flatnonzero(Temp[:,4])

        # plot TB_H's into subplot 'PlotTBH'
        PlotTBH = plt.subplot2grid((4,2),(0,0), rowspan=3)   # subplot for TB_H's in column 1, spanning 3 rows
        plt.plot(Temp[Unflagged,0], Temp[Unflagged,2], 'k.', label=r'$T_{\mathrm{B}}^{\mathrm{ H}}$ unflagged')
        plt.plot(Temp[Flagged,0],   Temp[Flagged,2],   'r.', label=r'$T_{\mathrm{B}}^{\mathrm{ H}}$ flagged')

        # format 'PlotTBH'
        plt.ylabel(r'$T_{\mathrm{B}}^\mathrm{ H}$ [K]')
        plt.legend()
        plt.title(r'$T_{\mathrm{B}}^{\mathrm{ H}}$ at horizontal polarization')

        # 1.3) plot flags at H-pol into lower subplot 'PlotFlagsH' of left column of Fig1
        PlotFlagsH = plt.subplot2grid((4,2),(3,0), rowspan=1, sharex = PlotTBH)   # subplot for flags at H-pol in column 1, spanning 1 row, shared x-axis with 'PlotTBH'
        # first derive all individual flags with 'SplitFlags'
        F1, F2, F3, F4, F5, F6 = SplitFlags(Temp[:,4])
        # then, plot individual flags into 'PlotFlagsH'
        plt.plot(Temp[:,0], F6*5,   linestyle='none', marker='x', markeredgecolor='m',      markerfacecolor='w', label=r'$T_{\mathrm{inst}}$ flag')              # Flag6
        plt.plot(Temp[:,0], F5*4,   linestyle='none', marker='x', markeredgecolor='orange', markerfacecolor='w', label=r'$T_{\mathrm{B}}^{\mathrm{ H}}$ flag')   # Flag5
        plt.plot(Temp[:,0], F4*3,   linestyle='none', marker='o', markeredgecolor='k',      markerfacecolor='w', label=r'$FD$ flag')                             # Flag4
        plt.plot(Temp[:,0], F3*2,   linestyle='none', marker='o', markeredgecolor='r',      markerfacecolor='w', label=r'$w_{\mathbf{U}}$ flag')                 # Flag3
        plt.plot(Temp[:,0], F2*1.5, linestyle='none', marker='o', markeredgecolor='g',      markerfacecolor='w', label=r'$s_{\mathbf{U}}$ flag')                 # Flag2
        plt.plot(Temp[:,0], F1,     linestyle='none', marker='o', markeredgecolor='b',      markerfacecolor='w', label=r'$k_{\mathbf{U}}$ flag')                 # Flag1
        # then format 'PlotFlagsH'
        PlotFlagsH.set_yticks((1, 1.5, 2, 3, 4, 5))   # set position of y-tick labels
        plt.ylim(0.5, 5.5)   # set y-axis limits
        labels = PlotFlagsH.set_yticklabels((r'$k_{\mathbf{U}}$', r'$s_{\mathbf{U}}$', r'$w_{\mathbf{U}}$', r'$FD$', r'$T_{\mathrm{B}}^{\mathrm{ H}}$', r'$T_{\mathrm{inst}}$'))   # define y-labels
        plt.xlabel('date of sky measurement')
        plt.ylabel('distortion check')
        plt.legend()
        Fig1.autofmt_xdate()   # format x-labels, so that dates are readable

        # 1.4) plot TB_V's into upper subplot of right column of Fig1
        # find indices of unflagged / flagged measurements:
        #    first, find all samples (if any) where no RFI analysis in the TD was carried out (i.e., flags start with -999) and remove the TD flags
        NoTDFlags = flatnonzero(Temp[:,7]<0)
        Temp[NoTDFlags,7] = (Temp[NoTDFlags,7]*-1) - 999000
        #    then, test which flags are 0
        Unflagged = flatnonzero(Temp[:,7]==0)
        Flagged   = flatnonzero(Temp[:,7])

        # plot TB_V's into subplot 'PlotTBV'
        PlotTBV = plt.subplot2grid((4,2),(0,1), rowspan=3)   # subplot for TB_Vs in column 2, spanning 3 rows
        plt.plot(Temp[Unflagged,0], Temp[Unflagged,5], 'k.', label=r'$T_{\mathrm{B}}^{\mathrm{ V}}$ unflagged')
        plt.plot(Temp[Flagged,0],   Temp[Flagged,5],   'r.', label=r'$T_{\mathrm{B}}^{\mathrm{ V}}$ flagged')

        # format 'PlotTBV'
        plt.ylabel(r'$T_{\mathrm{B}}^\mathrm{ V}$ [K]')
        plt.legend()
        plt.title(r'$T_{\mathrm{B}}^{\mathrm{ V}}$ at vertical polarization')

        # 1.5) plot flags at V-pol into lower subplot 'PlotFlagsV' of right column of Fig1
        PlotFlagsV = plt.subplot2grid((4,2),(3,1), rowspan=1, sharex = PlotTBV)   # subplot for flags at V-pol in column 2, spanning 1 row, shared x-axis with 'PlotTBV'
        # first, derive all individual flagswith 'SplitFlags'
        F1, F2, F3, F4, F5, F6 = SplitFlags(Temp[:,7])
        #then, plot individual flags into 'PlotFlagsV'
        plt.plot(Temp[:,0], F6*5,   linestyle='none', marker='x', markeredgecolor='m',      markerfacecolor='w', label=r'$T_{\mathrm{inst}}$ flag')              # Flag6
        plt.plot(Temp[:,0], F5*4,   linestyle='none', marker='x', markeredgecolor='orange', markerfacecolor='w', label=r'$T_{\mathrm{B}}^{\mathrm{ V}}$ flag')   # Flag5
        plt.plot(Temp[:,0], F4*3,   linestyle='none', marker='o', markeredgecolor='k',      markerfacecolor='w', label=r'$FD$ flag')                             # Flag4
        plt.plot(Temp[:,0], F3*2,   linestyle='none', marker='o', markeredgecolor='r',      markerfacecolor='w', label=r'$w_{\mathbf{U}}$ flag')                 # Flag3
        plt.plot(Temp[:,0], F2*1.5, linestyle='none', marker='o', markeredgecolor='g',      markerfacecolor='w', label=r'$s_{\mathbf{U}}$ flag')                 # Flag2
        plt.plot(Temp[:,0], F1,     linestyle='none', marker='o', markeredgecolor='b',      markerfacecolor='w', label=r'$k_{\mathbf{U}}$ flag')                 # Flag1
        # then, format 'PlotFlagsV'
        PlotFlagsV.set_yticks((1, 1.5, 2, 3, 4, 5))   # set position of y-tick labels
        plt.ylim(0.5, 5.5)   # set y-axis limits
        labels = PlotFlagsV.set_yticklabels((r'$k_{\mathbf{U}}$', r'$s_{\mathbf{U}}$', r'$w_{\mathbf{U}}$', r'$FD$', r'$T_{\mathrm{B}}^{\mathrm{ V}}$', r'$T_{\mathrm{inst}}$'))   # define y-labels
        plt.xlabel('date of sky measurement')
        plt.ylabel('distortion check')
        plt.legend()
        Fig1.autofmt_xdate()   # format x-labels, so that dates are readable

        plt.subplots_adjust(bottom = 0.13)   # adjust top and bottom position of all subplots in Fig1
        plt.show()
        plt.close('all')
