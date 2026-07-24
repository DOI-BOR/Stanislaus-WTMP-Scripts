#version 2.0
#modified 03-28-2023 by Scott Burdick-Yahya

from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath
from com.rma.model import Project
import os,shutil,copy,sys
from java.util import Vector, Date

import datetime
from hec.heclib.util import HecTime  # Replace with actual import


def standardize_interval(tsm, interval, makePerAver=True):
    
    """
    Resamples a time-series math object to a target interval and optionally sets its type to PER-AVER.

    If the time series is already at the target interval, it is returned unchanged.
    Otherwise, it is transformed using time-averaging to the requested interval.

    Inputs:
      tsm         -- HEC TimeSeriesMath object to standardize
      interval    -- target interval string: '1hour', '1day', or '1week' (case-insensitive)
      makePerAver -- if True and a transform is needed, sets the type to PER-AVER before resampling (default True)

    Output:
      Returns a TimeSeriesMath object at the target interval.
      Calls sys.exit(-1) if the interval string is not recognized.
    """
    
    # Extract the underlying TimeSeriesContainer to inspect the current interval
    tsc = tsm.getData()
    
    # Map the interval string to its integer equivalent in minutes
    if interval.lower()=='1hour':
        intint = 60
    elif interval.lower()=='1day':
        intint=1440
    elif interval.lower()=='1week':
        intint=10800
    else:
        # Unsupported interval: log and abort
        print('interval not supported:',interval)
        sys.exit(-1)

    # Only transform if the current interval differs from the target
    if tsc.interval != intint:
        if makePerAver:
            #tsc.type = 'PER-AVER'  # make sure it's per-aver ... we are
            # Set the data type to PER-AVER before resampling
            tsm.setType('PER-AVER')
        # Resample to the target interval using time-averaging
        return tsm.transformTimeSeries(interval, "", "AVE")
    else:
        # Already at the target interval; return as-is
        return tsm


def data_from_dss(dss_file,dss_rec,starttime_str, endtime_str,):
    """
    Reads a DSS time-series record and returns its values array.

    Inputs:
      dss_file      -- full path to the DSS file to read from
      dss_rec       -- DSS pathname string of the record to read
      starttime_str -- start date/time string in HEC format (e.g., '01Jan2014 0000'), or None for full record
      endtime_str   -- end date/time string in HEC format, or None for full record

    Output:
      Returns the values array from the TimeSeriesContainer for the requested record and time range.
    """
    
    # Open the DSS file, read the record for the given time range, and extract values
    dssFm = HecDss.open(dss_file)        
    tsc = dssFm.read(dss_rec, starttime_str, endtime_str, False).getData()
    dssFm.close()
    return tsc.values
 

def hectime_to_datetime(tsc):
    """
    Converts HEC integer time values from a TimeSeriesContainer to a list of Python datetime objects.

    Inputs:
      tsc -- HEC TimeSeriesContainer object containing numberValues and getHecTime() method

    Output:
      Returns a list of Python datetime objects corresponding to each time step in the container.
    """
    
    # Initialize the output list of datetime objects
    dtt = []
    for j in range(tsc.numberValues):
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        # Convert the HEC time at index j to a Java Date object (UTC offset 0)
        java_date = tsc.getHecTime(j).getJavaDate(0)  
        
        # Convert Java Date to Python datetime via Unix timestamp (milliseconds to seconds)
        timestamp = (java_date.getTime() / 1000)
        dtt.append(datetime.datetime.fromtimestamp(timestamp))

    return dtt

def fixInputLocationFpart(currentAlternative, tspath):
    """
    Updates the F-part of a DSS pathname to match the current alternative's input F-part prefix,
    preserving the last colon-separated token of the original F-part.

    Inputs:
      currentAlternative -- WAT scripting alternative object providing getInputFPart()
      tspath             -- DSS pathname string whose F-part is to be updated

    Output:
      Returns the updated DSS pathname string with the corrected F-part.
    """
    
    # Build the new F-part prefix from all but the last colon-separated token of the alternative's input F-part
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    
    # Split the pathname into its slash-delimited parts
    tspath = tspath.split('/')
    
    # Extract the existing F-part (index 6) and its last colon-separated token
    fpart = tspath[6]
    fpart_split = fpart.split(':')
    
    # Combine the new prefix with the original trailing token
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    
    # Replace the F-part in the pathname parts and reassemble
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    return tspath

def appendAPart(current_path, ApartAppend):
    """
    Appends a suffix to the A-part of a DSS pathname, separated by an underscore.
    If the A-part is empty, the suffix becomes the entire A-part.

    Inputs:
      current_path  -- DSS pathname string whose A-part is to be modified
                       NOTE: contains a bug - uses 'tspath' instead of 'current_path' (NameError at runtime)
      ApartAppend   -- string to append to the existing A-part

    Output:
      Returns the updated DSS pathname string with the modified A-part.
    """
    
    # Split the pathname into slash-delimited parts
    # BUG (pre-existing): 'tspath' is not defined; should be 'current_path'
    tspath = tspath.split('/')
    
    # Extract the existing A-part (index 1)
    Apart = tspath[1]
    
    # Build the new A-part: use the suffix alone if A-part is empty, otherwise append with underscore
    if len(Apart) == 0:
        new_Apart = ApartAppend
    else:
        new_Apart = Apart + '_' + ApartAppend
    
    # Replace the A-part and reassemble the pathname
    tspath[1] = new_Apart
    tspath = '/'.join(tspath)
    return tspath

def getDataLocationDSSInfo(location, currentAlternative, computeOptions):
    """
    Resolves the DSS pathname and DSS file path for a given WAT data location,
    handling both model-linked and externally linked locations.

    Inputs:
      location           -- WAT DataLocation object describing the data source
      currentAlternative -- WAT scripting alternative object for resolving linked time series
      computeOptions     -- WAT compute options object providing the compute DSS filename

    Output:
      Returns a tuple (tspath, dsspath) where:
        tspath   -- the resolved DSS pathname string
        dsspath  -- the full path to the DSS file containing the record
    """
    
    # Branch on whether the location is linked to output from a previous model
    if location.isLinkedToPreviousModel():
        # Resolve the DSS pathname from the alternative's linked time series
        tspath = str(currentAlternative.loadTimeSeries(location))
        
        # Fix the F-part to match the current alternative's input context
        tspath = fixInputLocationFpart(currentAlternative, tspath)
        
        # Use the compute DSS file as the data source
        dsspath = computeOptions.getDssFilename()
    else:
        # Use the externally specified DSS path and file from the linked location
        tspath = location.getLinkedToLocation().getDssPath()
        rundir = Project.getCurrentProject().getProjectDirectory()
        dsspath = location.getLinkedToLocation().get_dssFile()
        
        # Resolve the DSS file path relative to the project directory
        dsspath = os.path.join(rundir, dsspath)
    return tspath, dsspath

def strip_templateID_and_rename_records(dssFilePath,currentAlt):
    """
        Strips the first 4 characters from the F-part of all DSS record pathnames in a file,
        effectively removing a template ID prefix. Creates a backup of the original DSS file first.

        If any record's F-part does not contain a hyphen, the function returns early without changes.

        Inputs:
          dssFilePath -- full path to the DSS file whose records are to be renamed
          currentAlt  -- WAT scripting alternative object used for logging renamed paths

        Output:
          No return value. Renames records in-place in the DSS file; a backup is saved at dssFilePath+'.bak'.
    """

    # make copy of dss file
    shutil.copyfile(dssFilePath,dssFilePath+'.bak')

    # Open the DSS file and retrieve all record pathnames for processing
    # rename all records, stripping first 4 chars from f-part
    dss = HecDss.open(dssFilePath)
    rec_names = dss.getPathnameList()
    
    # Prepare a Java Vector to hold the renamed pathnames (required by HecDss.renameRecords)
    new_rec_names = Vector()
    #currentAlt.addComputeMessage(type(rec_names).__name__)
    
    
    for i,r in enumerate(rec_names):                
        
        # Split the pathname into its slash-delimited parts
        parts = r.split('/')
        
        # If the F-part contains no hyphen, assume no template ID prefix; abort early
        if not '-' in parts[-2]:
            return
        
        # Strip the first 4 characters (template ID) from the F-part
        parts[-2] = parts[-2][4:]
        
        # Build the renamed pathname and add to the output vector
        new_rec_names.add('/'.join(parts))
        currentAlt.addComputeMessage('Fixing path: '+r+' --> '+new_rec_names[-1])
    
    # Perform the batch rename operation on the DSS file
    dss.renameRecords(rec_names, new_rec_names)
    
    dss.close()

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    """
    Reads multiple DSS time-series records and writes their element-wise sum to a new DSS record.

    Inputs:
      currentAlt   -- WAT scripting alternative object used for logging
      dssFile      -- full path to the DSS file containing the input records
      timewindow   -- WAT run time window object providing start/end time strings
      input_data   -- list of DSS pathname strings to read and sum
      output_path  -- DSS pathname string for the output summed record

    Output:
      Returns 0 on successful completion.
      Writes the summed time-series record to dssFile at output_path.
    """
    
    # Extract start and end time strings from the run time window
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the DSS file for reading input records
    dssFm = HecDss.open(dssFile)
    output_data = []
    
    # Read and accumulate each input DSS record element-wise
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        dsstype = ts.type
        
        # First record initializes the accumulator; subsequent records are added element-wise
        if len(output_data) == 0:
            output_data = values
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val
    
    # Build the output TimeSeriesContainer with the summed values    
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_path
    tsc.values = output_data
    tsc.startTime = times[0]
    tsc.units = units
    tsc.type = dsstype
    tsc.endTime = times[-1]
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    
    # Write the summed record to the DSS file and close
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod, prepend_first_value=False):
    """
    Resamples a DSS time-series record to a new time interval using time-averaging.
    Optionally prepends the first source value to the resampled output to prevent
    loss of the initial time step during upsampling (e.g., 1DAY -> 1HOUR).

    Inputs:
      inputDSSFile        -- full path to the source DSS file
      inputRec            -- DSS pathname string of the record to resample
      timewindow          -- WAT run time window object providing start/end time strings
      outputDSSFile       -- full path to the destination DSS file for the resampled record
      newPeriod           -- target interval string (e.g., '1HOUR', '1DAY')
      prepend_first_value -- if True, prepend the first source value one step before the
                             resampled start time to ensure compute completeness (default False)

    Output:
      No return value. Writes the resampled time-series record to outputDSSFile.
      NOTE: contains a pre-existing bug - 'dssFm_out' is used but not defined in the
      prepend_first_value branch; should be 'dssFmout'.
    """
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR'''

    # Open the source DSS file and read the record for the run time window
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    
    # Resample to the target interval using time-averaging
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")
    if not prepend_first_value:
        # simple tsm transform
        # Write the resampled record directly to the output DSS file
        dssFm.close()
        dssFmout = HecDss.open(outputDSSFile)
        dssFmout.write(tsm_new)
        dssFmout.close()
    else:
       # sometimes when upsampling, like from 1Day -> 1Hour, the tsm transform lops off the first hour value that 
       # simulations need. Here we append the previous value as the first value of the new record.  This comes in handy
       # to make sure the compute record os complete.

       # Retrieve the original source container to access its first value
       tsc_orig = tsm.getData()
       dssFm.close()
       
       # Retrieve the resampled container and compute the time step size
       tsc_new = tsm_new.getData()
       steptime = tsc_new.times[1]-tsc_new.times[0]
       
       # Build a new container starting one step before the resampled start, prepending the source's first value
       tsc = TimeSeriesContainer()
       tsc.startTime = tsc_new.times[0] - steptime
       tsc.interval = tsc_new.interval
       tsc.fullName = tsc_new.fullName
       tsc.values = [tsc_orig.values[0],] + tsc_new.values
       tsc.units = tsc_new.units
       tsc.type = 'PER-AVER'
       tsc.numberValues = len(tsc.values)
       
       # BUG (pre-existing): 'dssFm_out' is not defined here; should be 'dssFmout'
       dssFm_out.write(tsc)



def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    """
    Applies an air temperature lapse rate correction to a DSS temperature record
    and writes the result with an updated F-part to an output DSS file.

    If the source record's units are in Fahrenheit, the lapse value is converted
    from Celsius to Fahrenheit before being applied.

    Inputs:
      dss_file   -- full path to the source DSS file containing the temperature record
      dss_rec    -- DSS pathname string of the temperature record to adjust
      lapse_in_C -- lapse rate correction value in degrees Celsius (added to each time step)
      dss_outfile-- full path to the destination DSS file for the corrected record
      f_part     -- new F-part string to apply to the output DSS pathname

    Output:
      No return value. Writes the lapse-corrected temperature record to dss_outfile.
    """
    
    # Open the source DSS file and read the temperature time-series record
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    
    # Convert the lapse rate to Fahrenheit if the record's units are in Fahrenheit
    lapse = lapse_in_C
    if 'f' in tsm.getUnits().lower():
        lapse = lapse*9.0/5.0+32.0
    
    # Add the (possibly converted) lapse value to all time steps in the record
    tsm = tsm.add(lapse)
    tsc = tsm.getData()
    dss.close()

    # Update the F-part of the output pathname and write to the output DSS file
    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()

def min_ts_flow_cfs(dss_file,dss_rec,min_value_cfs,dss_outfile,f_part):
    """
    Enforces a minimum flow value on a DSS flow record and writes the result
    with an updated F-part to an output DSS file.

    Handles unit conversion from CMS to CFS when the source record is in metric units.

    Inputs:
      dss_file      -- full path to the source DSS file containing the flow record
      dss_rec       -- DSS pathname string of the flow record to process
      min_value_cfs -- minimum allowable flow value in CFS
      dss_outfile   -- full path to the destination DSS file for the corrected record
      f_part        -- new F-part string to apply to the output DSS pathname

    Output:
      No return value. Writes the minimum-constrained flow record to dss_outfile.
      Calls sys.exit(-1) if the record's units are neither 'cms' nor 'cfs'.
    """
    
    # Open the source DSS file and read the flow time-series record
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    tsc = tsm.getData()
    dss.close()

    # Convert the minimum threshold to match the record's units
    if tsc.units.lower() == 'cms':
        # Convert CFS minimum to CMS for comparison
        min_value = min_value_cfs*0.028316847      
    elif tsc.units.lower() == 'cfs':
        # No conversion needed
        min_value = min_value_cfs
    else:
        # Unsupported units: log and abort
        print('min_ts_flow_cfs: flow units not understood (%s).'%tsc.units)
        sys.exit(-1)        

    # Apply the minimum constraint to each value in the record
    for vi, v in enumerate(tsc.values):
        tsc.values[vi] = max(v, min_value)

    # Update the F-part of the output pathname and write to the output DSS file
    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()

def add_flows(currentAlt, timewindow, inflow_records, dss_file, output_dss_record_name, output_dss_file):
     
    """
    Reads multiple DSS flow records and writes their element-wise sum to a new DSS record.

    Trims records that extend beyond the run time window on either end.
    Converts CMS values to CFS automatically. Units are assumed to be homogeneous across records.

    Inputs:
      currentAlt             -- WAT scripting alternative object for logging
      timewindow             -- WAT run time window object providing start/end time strings
      inflow_records         -- list of DSS pathname strings to read and sum
      dss_file               -- full path to the DSS file containing the input records
      output_dss_record_name -- DSS pathname string for the output summed record
      output_dss_file        -- full path to the DSS file where the output record is written

    Output:
      No return value. Writes the summed flow record (in CFS) to output_dss_file.
      Calls sys.exit(-1) if any record cannot be read.
    """

    # Extract start and end time strings and convert to HecTime integer values for bounds-checking
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the source DSS file for reading
    dssFm = HecDss.open(dss_file)

    # Initialize accumulators for summed flow values and the shared time axis
    inflows = []
    times = []

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
            # Read the flow record for the run time window
            print(starttime_str, endtime_str)
            print(dss_file)
            ts = dssFm.read(pathname, starttime_str, endtime_str, False)
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            tstype = ts_data.type
            
            # Trim leading values that fall before the run time window start
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
                
            # Trim trailing values that fall after the run time window end    
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
                

        except HecMathException:
            # Fatal error: log and abort if the record cannot be read
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)
        
        # Convert CMS to CFS if the source record is in metric units
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # First record initializes the accumulator; subsequent records are summed element-wise
        if len(inflows) == 0:
            inflows = values
            times = hectimes 
        else:
            for vi, v in enumerate(values):
                inflows[vi] += v

    # Build the output TimeSeriesContainer with the summed flow values
    # Output record
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    tsc.units = 'cfs'
    tsc.type = tstype
    tsc.numberValues = len(inflows)
    
    # Open the output DSS file, write the result, and close both files
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()


def add_or_subtract_flows(currentAlt, timewindow, inflow_records, dss_file, operation,
                       output_dss_record_name, output_dss_file):
    """
    Reads multiple DSS flow records and combines them element-wise using per-record
    add or subtract operations, writing the result to a new DSS record.

    Trims records that extend beyond the run time window on either end.
    Converts CMS values to CFS automatically.

    Inputs:
      currentAlt             -- WAT scripting alternative object for logging
      timewindow             -- WAT run time window object providing start/end time strings
      inflow_records         -- list of DSS pathname strings to combine
      dss_file               -- full path to the DSS file containing the input records
      operation              -- list of booleans/None parallel to inflow_records:
                                  None or True  = add this record to the accumulator
                                  False         = subtract this record from the accumulator
                                  (first record always initializes the accumulator regardless)
      output_dss_record_name -- DSS pathname string for the output record
      output_dss_file        -- full path to the DSS file where the output record is written

    Output:
      No return value. Writes the combined flow record (in CFS) to output_dss_file.
      Calls sys.exit(-1) if any record cannot be read.
    """

    # Extract start and end time strings and convert to HecTime integer values for bounds-checking
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the source DSS file for reading
    dssFm = HecDss.open(dss_file)

    # Initialize accumulators for the combined flow values and the shared time axis
    inflows = []
    times = []

    # Read inflows
    print('Reading inflows')
    
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
            # Read the flow record for the run time window
            print(starttime_str, endtime_str)
            print(dss_file)
            ts = dssFm.read(pathname, starttime_str, endtime_str, False)
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            tstype = ts_data.type

            # Trim leading values that fall before the run time window start
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            
            # Trim trailing values that fall after the run time window end
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
                

        except HecMathException:
            # Fatal error: log and abort if the record cannot be read
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # Convert CMS to CFS if the source record is in metric units
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # First record always initializes the accumulator
        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            # Apply per-record add or subtract based on the operation flag
            if operation[j]:
                # operation[j] is True or None: add this record's values
                for vi, v in enumerate(values):
                    inflows[vi] += v
            else:
                # operation[j] is False: subtract this record's values
                for vi, v in enumerate(values):
                    inflows[vi] -= v
                    
    # Build the output TimeSeriesContainer with the combined flow values
    # Output record
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    tsc.units = 'cfs'
    tsc.type = tstype
    tsc.numberValues = len(inflows)
    
    # Open the output DSS file, write the result, and close both files
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()


def hec_str_time_to_dt(hec_str_time):
    """
    Converts a HEC-formatted date/time string to a Python datetime object.

    Handles the HEC convention of '2400' end-of-day timestamps by converting
    them to midnight of the following day.

    Inputs:
      hec_str_time -- HEC date/time string in the format 'ddMmmYYYY HHMM' (e.g., '01Jan2014 2400')

    Output:
      Returns a Python datetime object representing the equivalent date and time.
    """
    
    '''Convert HEC date time format to python datetime object'''

    # HEC datetime format string for strptime parsing
    dt_format = '%d%b%Y %H%M'
    
    # Flag to add one day when the HEC time uses the '2400' end-of-day convention
    add_day = False
    if hec_str_time.endswith('2400'):
        # Replace '2400' with '0000' so strptime can parse it, then add one day afterward
        my_hec_str_time = hec_str_time[:-4] + '0000'
        add_day = True
    else:
        # Standard time string; no adjustment needed
        my_hec_str_time = hec_str_time

    # Parse the (possibly adjusted) string to a datetime object
    dt = datetime.datetime.strptime(my_hec_str_time,dt_format)
    
    # Advance by one day to account for the '2400' end-of-day convention
    if add_day:
        dt = dt + datetime.timedelta(days=1)
    
    return dt


def create_constant_dss_rec(currentAlt, timewindow, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS', fpart='ZEROS'):
    
    """
    Creates and writes a DSS time-series record filled with a constant value
    for the run time window, padded by one day on each end.

    Supports flow, water temperature, and gate parameter types.

    Inputs:
      currentAlt      -- WAT scripting alternative object for logging error messages
      timewindow      -- WAT run time window object providing start/end time strings
      output_dss_file -- full path to the DSS file where the constant record is written
      constant        -- constant value to fill the record with (default 0.0)
      what            -- parameter type string: 'flow', 'temp-water', or 'gate' (default 'flow')
      dss_type        -- DSS data type string (e.g., 'PER-AVER', 'INST-VAL') (default 'PER-AVER')
      period          -- time interval string: '1HOUR' or '1DAY' (default '1HOUR')
      cpart           -- C-part (location) string for the DSS pathname (default 'ZEROS')
      fpart           -- F-part (version) string for the DSS pathname (default 'ZEROS')

    Output:
      Returns True on successful write.
      Returns False if 'what' or 'period' is not recognized (and logs an error message).
    """
    
    
    '''Create and write a dss record with a constant in it for the given time windows.
       what={'flow','temp-water'}
       period={'1HOUR','1DAY'}
    '''

    # Resolve units and parameter string based on the requested data type
    if what.lower()=='flow':
        units = 'cfs'
        parameter = 'flow'
    elif what.lower()=='temp-water':
        units = 'C'
        parameter = 'temp-water'
    elif what.lower()=='gate':
        units = 'n/a'
        parameter = 'gate'
    else:
        # Unsupported parameter type: log and return failure
        currentAlt.addComputeMessage('create_zero_dss_rec: what not known: %s'%what)
        return False

    # Validate the requested time period
    if period.lower()=='1hour':
        pass
    elif period.lower()=='1day':
        pass
    else:
        # Unsupported period: log and return failure
        currentAlt.addComputeMessage('create_zero_dss_rec: period not known: %s'%period)
        return False

    # HEC datetime format string for formatting padded start/end times
    dt_format = '%d%b%Y %H%M'
    
    # Extract the run time window start and end strings
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()

    # pad 1 day on records, in case these are used for lookbacks, or balance flow calcs, etc.
    # Pad the time window by one day on each end to support lookback and balance flow use cases
    starttime_dt = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=1)    
    endtime_dt = hec_str_time_to_dt(endtime_str) + datetime.timedelta(days=1)
    starttime_str_pad = starttime_dt.strftime(dt_format)
    endtime_str_pad = endtime_dt.strftime(dt_format)    
 
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    ########################
    # Zero-Flow Time Series
    ########################

    # Generate a regular interval time series filled with the constant value
    tsmath_zero_flow_day = tsmath.generateRegularIntervalTimeSeries(
        starttime_str_pad,
        endtime_str_pad,
        period, "0M", constant)
    
    # Set metadata on the generated time series: units, type, interval, location, parameter, and version
    tsmath_zero_flow_day.setUnits(units)
    tsmath_zero_flow_day.setType(dss_type)
    tsmath_zero_flow_day.setTimeInterval(period)
    tsmath_zero_flow_day.setLocation(cpart)
    tsmath_zero_flow_day.setParameterPart(parameter)
    tsmath_zero_flow_day.setVersion(fpart)

    # Open the output DSS file, write the constant record, and close
    dssFm = HecDss.open(output_dss_file)
    dssFm.write(tsmath_zero_flow_day)
    dssFm.close()

    return True


def calculate_relative_humidity(air_temp, dewpoint_temp):
    """
    Calculate Relative Humidity given the air temperature and dewpoint temperature - August-Roche-Magnus approximation

    :param air_temp: Air Temperature in degrees Celsius
    :param dewpoint_temp: Dew Point Temperature in degrees Celsius
    :return: Relative Humidity in percentage
    """
    
    # Compute the numerator and denominator of the Magnus-formula scaling factor
    numerator = (112.0 - 0.1 * dewpoint_temp + air_temp)
    denominator = (112.0 + 0.9 * air_temp)
    
    # Compute the exponent term: difference of saturation vapor pressure log-approximations
    exponent = ((17.62 * dewpoint_temp) / (243.12 + dewpoint_temp)) - ((17.62 * air_temp) / (243.12 + air_temp))
    
    # Combine terms to produce relative humidity, then clamp to the valid range [0.01, 100.0]
    relative_humidity = 100.0 * (numerator / denominator) * math.exp(exponent)
    return max(0.01, min(100.0, relative_humidity))


def relhum_from_at_dp(met_dss_file, at_path, dp_path):
    
    """
    Computes relative humidity from air temperature and dewpoint temperature DSS records
    and writes the derived record back to the meteorology DSS file.

    The output pathname is derived from the air temperature record's pathname:
      - Location (B-part) is truncated to 5 characters
      - Parameter (C-part) is set to 'RELHUM-FROM-AT-DP'
      - Version (F-part) is appended with '-DERIVED'

    Inputs:
      met_dss_file -- full path to the meteorology DSS file containing both input records
      at_path      -- DSS pathname string for the air temperature record (degrees C)
      dp_path      -- DSS pathname string for the dewpoint temperature record (degrees C)

    Output:
      No return value. Writes the derived relative humidity record (%) to met_dss_file.
    """
    
    # Open the meteorology DSS file and read the air temperature record
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    
    # Read the dewpoint temperature values using data_from_dss (no time window filtering)
    dp_data = DSS_Tools.data_from_dss(met_dss_file, dp_path, None, None)
    
    # Compute relative humidity at each time step from air temp and dewpoint temp
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_relative_humidity(tsc.values[i], dp_data[i])
    
    # Build the output pathname from the air temperature record's pathname parts
    parts = tsc.fullName.split('/')
    # Truncate B-part (location) to 5 characters
    parts[2] = parts[2][:5]
    # Set C-part (parameter) to the derived humidity label
    parts[3] = 'RELHUM-FROM-AT-DP'
    # Append '-DERIVED' to the F-part (version)
    parts[6] = parts[6] + '-DERIVED'
    new_pathname = '/'.join(parts)
    
    # Assign the new pathname and units to the output container
    tsc.fullName = new_pathname
    tsc.units = '%'
    
    # Write the derived relative humidity record to the DSS file
    print('writing: ', new_pathname)
    dss.write(tsc)
    dss.close()
