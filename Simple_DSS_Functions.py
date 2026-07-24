from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    """
    Reads multiple DSS time-series records from a single DSS file and writes
    their element-wise sum to a new record in the same file.

    Units and data type are taken from the last record read. All input records
    must share the same time axis and units for the result to be meaningful.

    Inputs:
      currentAlt   -- WAT scripting alternative object used for logging
      dssFile      -- full path to the DSS file containing both input and output records
      timewindow   -- WAT run time window object providing start/end time strings and HecTime values
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
    
    # Open the DSS file for both reading input records and writing
    dssFm = HecDss.open(dssFile)
    output_data = []
    
    # Read and accumulate each input DSS record element-wise
    for dsspath in input_data:
        print('reading', str(dsspath))
        
        # Read the record for the run time window and extract the underlying container
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        
        # First record initializes the accumulator; subsequent records are added element-wise
        if len(output_data) == 0:
            output_data = values
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val
      
    # Build the output TimeSeriesContainer with the summed values and metadata
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_path
    tsc.values = output_data
    tsc.startTime = times[0]
    tsc.units = units
    tsc.endTime = times[-1]
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    
    # Write the summed record to the DSS file, close, and log the result
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod, 
                    prepend_first_value=False,inst_val=False):
                        
    """
    Resamples a DSS time-series record to a new time interval using time-averaging.

    Optionally prepends the first value of the source record one step before the
    resampled output's start time to prevent loss of the initial time step during
    upsampling (e.g., 1DAY -> 1HOUR). Optionally forces the output data type to INST-VAL.

    Inputs:
      inputDSSFile        -- full path to the source DSS file containing the record to resample
      inputRec            -- DSS pathname string of the record to resample
      timewindow          -- WAT run time window object providing start/end time strings
      outputDSSFile       -- full path to the destination DSS file for the resampled record
      newPeriod           -- target interval string (e.g., '1HOUR', '1DAY')
      prepend_first_value -- if True, prepend the source record's first value to the resampled
                             output to ensure compute completeness (default False)
      inst_val            -- if True, sets the output record's data type to 'INST-VAL'
                             regardless of the source type; used for instantaneous records
                             such as reservoir stage (default False)

    Output:
      No return value. Writes the resampled time-series record to outputDSSFile.
      NOTE: contains a pre-existing tab/space indentation mix on the 'tsc.type = INST-VAL'
      line inside the prepend branch, which may cause an IndentationError in Python 3.
    """
    
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR'''
    
    # Open the source DSS file and read the record for the run time window
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    
    # Resample the source record to the target interval using time-averaging
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")

    # Open the output DSS file for writing the resampled record
    dssFmout = HecDss.open(outputDSSFile)
    if not prepend_first_value:
        # simple tsm transform
        # Write the resampled record directly to the output DSS file and close the source
        dssFm.close()
        dssFmout.write(tsm_new)
    else:
       # sometimes when upsampling, like from 1Day -> 1Hour, the tsm transform lops off the first hour value that 
       # simulations need. Here we append the previous value as the first value of the new record.  This comes in handy
       # to make sure the compute record os complete.

       # Retrieve the original source container to access its first value for prepending
       tsc_orig = tsm.getData()
       dssFm.close()
       
       # Retrieve the resampled container and compute the step size for reference
       tsc_new = tsm_new.getData()
       steptime = tsc_new.times[1]-tsc_new.times[0]
       
       # Build the output container, reusing the resampled record's time axis and metadata
       tsc = TimeSeriesContainer()
       tsc.times = tsc_new.times
       tsc.fullName = tsc_new.fullName
       
       # Extract resampled values for manual list construction
       new_values = tsc_new.values
       
       # Prepend the source record's first value; iterate remainder to convert array.array -> list
       values = [tsc_orig.values[0],]
       for i in range(1,len(new_values)):
           values.append(new_values[i]) # this is inane, but I'm having trouble getting this to work otherwise, converting form array.array - list
       
       # Assign the manually constructed values list to the output container
       tsc.values = values     
       tsc.units = tsc_new.units
       
       # Set data type: INST-VAL for instantaneous records (e.g., stage), or preserve source type
       if inst_val:
	       tsc.type = 'INST-VAL'
       else:
           tsc.type = tsc_new.type
       
       tsc.numberValues = len(tsc.values)
       
       # Write the prepended, resampled record to the output DSS file
       dssFmout.write(tsc)
    
    # Close the output DSS file after writing
    dssFmout.close()
