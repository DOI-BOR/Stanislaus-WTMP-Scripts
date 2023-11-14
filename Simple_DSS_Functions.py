from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    output_data = []
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        if len(output_data) == 0:
            output_data = values
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val
                
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
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod, 
                    prepend_first_value=False,inst_val=False):
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR'''
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")

    dssFmout = HecDss.open(outputDSSFile)
    if not prepend_first_value:
        # simple tsm transform
        dssFm.close()
        dssFmout.write(tsm_new)
    else:
       # sometimes when upsampling, like from 1Day -> 1Hour, the tsm transform lops off the first hour value that 
       # simulations need. Here we append the previous value as the first value of the new record.  This comes in handy
       # to make sure the compute record os complete.
       tsc_orig = tsm.getData()
       dssFm.close()
       tsc_new = tsm_new.getData()
       steptime = tsc_new.times[1]-tsc_new.times[0]
       
       tsc = TimeSeriesContainer()
       tsc.times = tsc_new.times
       tsc.fullName = tsc_new.fullName
       new_values = tsc_new.values
       values = [tsc_orig.values[0],]
       for i in range(1,len(new_values)):
           values.append(new_values[i]) # this is inane, but I'm having trouble getting this to work otherwise, converting form array.array - list
       tsc.values = values     
       tsc.units = tsc_new.units
       if inst_val:
	       tsc.type = 'INST-VAL'
       else:
           tsc.type = tsc_new.type
       tsc.numberValues = len(tsc.values)
       dssFmout.write(tsc)
    dssFmout.close()
