"""
Created on 4/28/2022
@note:
modified for jython to be used in WAT by SBY on 8/3/2023
"""

# import datetime as dt
# import numpy as np
# from pydsstools.heclib.dss import HecDss
# from pydsstools.core import TimeSeriesContainer
# from scipy import     
# import pandas as pd

import math
from hec.heclib.dss import HecDss
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.heclib.util import HecTime
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
import hec.hecmath.TimeSeriesMath as tsmath
from rma.util.RMAConst import MISSING_DOUBLE
import math
import sys
import datetime as dt
import os

from hec.io import DSSIdentifier
from hec.heclib.util import HecTime
from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

def linear_interpolation(x_values, y_values, x):
    
    """
    Performs linear interpolation to estimate a y value for a given x.

    Inputs:
      x_values -- list of x data points (must be sorted ascending, length >= 2)
      y_values -- list of corresponding y data points (same length as x_values)
      x        -- the x value at which to interpolate

    Output:
      Returns the interpolated y value (float).
      Raises ValueError if list lengths mismatch, fewer than 2 points, or x is out of range.
    """
    
    # Validate that input lists are the same length and have at least 2 points
    if len(x_values) != len(y_values) or len(x_values) < 2:
        raise ValueError("Input lists must have the same length and contain at least 2 data points.")

    # Iterate through consecutive pairs to find the bracketing interval for x
    for i in range(1, len(x_values)):
        if x <= x_values[i]:
            
            # Retrieve the lower and upper bounding (x, y) pairs
            x0, y0 = x_values[i - 1], y_values[i - 1]
            x1, y1 = x_values[i], y_values[i]

            # Perform linear interpolation
            y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)

            return y

    # If x is beyond the range of x_values, raise an error
    raise ValueError("Interpolation point is outside the range of provided data.")

def read_elev_storage_area_file(file_name, res_name):
    """
    Reads an elevation-storage-area lookup table from a CSV file.

    For most reservoirs, each row contains: elevation (ft), storage (acre-ft), area (acre).
    For Natoma reservoir, each row contains only: elevation (ft), area (acre) (no storage column).

    Inputs:
      file_name -- full path to the CSV file containing the elevation-storage-area table
      res_name  -- reservoir name string used to select the appropriate parsing format

    Output:
      Returns a dictionary with keys 'elev', 'stor', and 'area', each mapping to a list of floats.
      For Natoma, 'stor' will be an empty list.
    """
    
    # These are in [elev, stor, area] with units [ft, acre-ft, acre]
    elevstorarea = {} #avoid lists doing weird things like mixing up order..
    
    # Initialize separate lists for each column
    elev = []
    stor = []
    area = []
    
    # Print working directory for traceability during compute
    import os
    print('cwd: ' + os.getcwd())
    
    # Natoma has a two-column format (elev, area) - no storage column
    if res_name.lower() == 'natoma':
        with open(file_name, 'r') as fn:
            for line in fn:
                sline = line.strip().split(',')
                elev.append(float(sline[0]))
                area.append(float(sline[1]))
    else:
        # All other reservoirs use three-column format (elev, stor, area)
        with open(file_name, 'r') as fn:
            for line in fn:
                sline = line.strip().split(',')
                elev.append(float(sline[0]))
                stor.append(float(sline[1]))
                area.append(float(sline[2]))
                
    # Package results into a dictionary and return            
    elevstorarea['elev'] = elev
    elevstorarea['stor'] = stor
    elevstorarea['area'] = area
    return elevstorarea

def build_conic_storage_array(elev, area, firstStorageValue=0.0):
    """
    Computes cumulative storage at each elevation in the elevation-area curve
    using the conic (frustum) formula between consecutive measurement points.

    Adapted from storage.java in HEC ResSim (2022-06-17).

    Inputs:
      elev              -- list of elevation values (ft), sorted ascending
      area              -- list of surface area values (acre) corresponding to each elevation
      firstStorageValue -- storage value at the lowest elevation (default 0.0 acre-ft)

    Output:
      Returns a list of cumulative storage values (acre-ft), one per elevation point.
    """
    
    '''Find storage of slabs between measurement points on the elevation area curve,
    using a conic estimation.  Adapted from storage.java from HEC ResSim, 2022-06-17'''
    
    # Determine the number of elevation-area measurement pairs
    n_measures = len(elev)
    
    # Initialize storage array with the known base storage at the lowest elevation
    storage = []
    storage.append(firstStorageValue)
    
    # Apply conic frustum volume formula between each consecutive pair of elevation points
    for i in range(1, n_measures):
        # Height of the slab between consecutive elevation measurements
        h = elev[i] - elev[i-1]
        
        # Conic formula: V = h/3 * (A1 + A2 + sqrt(A1*A2)) + storage at base
        storage.append(h/3. * (area[i-1] + area[i] + math.sqrt(area[i-1] * area[i])) + storage[i-1])
        
    return storage


def conic_storage_interp(interpElev, elev, area, conicStorage, idx):
    """
    Interpolates storage at an arbitrary elevation using conic (frustum) interpolation
    between two bounding elevation-area measurement points.

    Adapted from storage.java in HEC ResSim (2022-06-17).

    Inputs:
      interpElev   -- the target elevation at which to interpolate storage (ft)
      elev         -- list of elevation values from the elevation-area table (ft)
      area         -- list of surface area values from the elevation-area table (acre)
      conicStorage -- list of precomputed conic storage values at each table elevation (acre-ft)
      idx          -- lower-bounding index into elev/area/conicStorage for interpElev

    Output:
      Returns the interpolated storage value at interpElev (acre-ft, float).
    """
    
    '''Find storage between measurement points on the elevation area curve,
    using interpolation between conic layers.  Adapted from storage.java from
    HEC ResSim, 2022-06-17'''
    
    # Fractional height within the bounding layer
    h = (interpElev - elev[idx]) / (elev[idx+1] - elev[idx])
    
    # Geometric mean of the two bounding areas (used in conic interpolation)
    geomMean = math.sqrt(area[idx] * area[idx+1])
    
    # Interpolate surface area at interpElev using conic area formula
    areaInterp = area[idx] + 2.*(geomMean - area[idx])*h + (area[idx] + area[idx+1] - 2.*geomMean)*h*h
    
    # Compute storage from the lower conic boundary up to interpElev using frustum formula
    storageInterp = (interpElev - elev[idx])/3. * (area[idx] + areaInterp + math.sqrt(area[idx] * areaInterp)) + conicStorage[idx]
    
    return storageInterp


def get_elev_layer_idx(elev, obs_elev, elev_stor_area):
    """
    Finds the lower-bounding index in the elevation-storage-area table for a given
    observed elevation (i.e., the index of the largest table elevation <= obs_elev).

    Inputs:
      elev          -- list of elevation values from the table (ft)
      obs_elev      -- the observed elevation to locate in the table (ft)
      elev_stor_area-- the full elevation-storage-area dictionary (used to verify bounding)

    Output:
      Returns the integer index of the lower-bounding elevation in the table,
      or -1 if a valid index cannot be determined.
    """
    
    # find lower bounding index of where elevation lands in elev-stor-area table
    # idx = np.argmin(np.abs(elev-obs_elev))
    # if elev_stor_area[idx, 0] > obs_elev:
    #     idx -= 1
    # return idx

    # Initialize index and minimum distance to UNDEFINED_DOUBLE sentinel
    idx = UNDEFINED_DOUBLE
    min_val = None
    
    # Iterate through all table elevations to find the closest one to obs_elev
    for i in range(len(elev)):
        valchk = abs(elev[i]-obs_elev) #TODO: is this multidimensional?
        
        # Handle NaN distance: mark as found but invalid
        if math.isnan(valchk):
            min_val = valchk
            idx = i
            
        # First valid distance initializes the minimum
        elif min_val == None:
            min_val = valchk
            idx = i
            
        # Update to a closer table elevation if found    
        elif valchk < min_val:
            min_val = valchk
            idx = i
    
    # If a closest index was found, step back one if the table value is above obs_elev
    if idx != UNDEFINED_DOUBLE:
        if elev_stor_area['elev'][idx] > obs_elev: 
            idx -= 1
    else:
        # No valid index could be found; return sentinel -1
        idx = -1
        
    return idx

def get_balance_period(balance_period):
    """
    Converts a balance period string to a duration in hours (float).

    Inputs:
      balance_period -- string describing the time step (e.g., '1HOUR', '1DAY', '30MIN')

    Output:
      Returns the equivalent period duration in hours as a float.
      Returns None implicitly if the string does not match any recognized unit.
    """
    
    # Parse hourly period: extract numeric portion and return directly as hours
    if 'hour' in balance_period.lower():
        return float(balance_period.lower().replace('hour', ''))
        
    # Parse daily period: convert days to hours by multiplying by 24
    elif 'day' in balance_period.lower():
        return float(balance_period.lower().replace('day', '')) * 24
        
    # Parse minute period: convert minutes to hours by dividing by 60    
    elif 'min' in balance_period.lower():
        return float(balance_period.lower().replace('min', '')) / 60

def check_dss_intervals(records, balance_period, currentAlt):
    """
    Validates that each DSS record pathname contains the expected time interval string.
    Exits the compute with an error message if any record does not match.

    Inputs:
      records        -- list of DSS pathname strings to validate
      balance_period -- expected time interval string (e.g., '1HOUR') to find in each pathname
      currentAlt     -- WAT alternative object used for posting error messages to the compute log

    Output:
      No return value. Calls sys.exit(-1) on the first mismatched record.
    """
    
    # Check each record pathname for the presence of the expected time interval token
    for r in records:
        if balance_period.lower() not in r.lower():
            # Log the mismatch and abort the compute
            currentAlt.addComputeMessage('DSS record {0} not matching time interval {1}'.format(r, balance_period))
            sys.exit(-1)


def read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str):
    """
    Reads a DSS time-series record, optionally from an alternate DSS file embedded
    in the pathname using '::' as a separator.

    If the pathname contains '::', the portion before '::' is treated as the DSS file
    path and a new connection is opened for that file. Otherwise, the provided dssFm
    connection is used.

    Inputs:
      dssFm         -- default open HecDss file manager object
      pathname      -- DSS pathname string, optionally prefixed with 'alt_dss_file::' 
      starttime_str -- start date/time string in HEC format (e.g., '01Jan2014 0000')
      endtime_str   -- end date/time string in HEC format

    Output:
      Returns the TimeSeriesContainer data object (.getData()) from the DSS read.
    """
    
    '''pathname may contain the dss filepath additionally before the dss ts path, separated by '::'
       If so, use that dss file.'''
       
    # Check whether the pathname embeds an alternate DSS file reference   
    if '::' in pathname:
        print('Splitting and reading:',pathname)
        
        # Split into the alternate file path and the clean DSS pathname
        alt_dss_file,pathname_clean = pathname.split('::')
        
        # Open the alternate DSS file and read the record
        dssFmRec = HecDss.open(alt_dss_file)
        tsc = dssFmRec.read(pathname_clean, starttime_str, endtime_str, False).getData()
        dssFmRec.close()
    else:
        # Use the default DSS file manager connection to read the record
        tsc = dssFm.read(pathname, starttime_str, endtime_str, False).getData()
    
    return tsc


def create_balance_flows(currentAlt, timewindow, res_name, inflow_records, outflow_records, stage_record, evap_record,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         storage_dss_record_name='', evap_dss_record_name='',
                         balance_period_str="1HOUR", use_conic=False, write_evap=False, write_storage=False,
                         alt_period=None,alt_period_string=None, lookback_padding=1440):

    """
    Computes water balance flows for a reservoir and writes results to a DSS file.

    For each time step, the balance flow is computed as:
      balance_flow = delta_storage_from_stage / dt - (inflow - outflow - evap_loss)

    The result represents the unexplained residual flow needed to close the water balance.

    Inputs:
      currentAlt             -- WAT scripting alternative; used for logging compute messages
      timewindow             -- WAT run time window object providing start/end time strings
      res_name               -- reservoir name string (used for CSV output naming)
      inflow_records         -- list of DSS pathnames for inflow time series (cfs, period-avg)
      outflow_records        -- list of DSS pathnames for outflow time series (cfs, period-avg)
      stage_record           -- DSS pathname for reservoir stage/elevation (ft, instantaneous)
      evap_record            -- DSS pathname for evaporation depth (ft, period-accumulated)
      elev_stor_area         -- dict with keys 'elev', 'stor', 'area' (lists of floats)
      dss_file               -- path to the source DSS file for reading inflows/stage/evap
      output_dss_record_name -- DSS pathname for the output balance flow record
      output_dss_file        -- path to the DSS file where output records are written
      shared_dir             -- directory for CSV diagnostic output
      storage_dss_record_name-- DSS pathname for optional storage output record (default '')
      evap_dss_record_name   -- DSS pathname for optional evaporation output record (default '')
      balance_period_str     -- time step string for balance computation (default '1HOUR')
      use_conic              -- if True, use conic interpolation for storage; else linear (default False)
      write_evap             -- if True, write derived evaporation flow to DSS (default False)
      write_storage          -- if True, write derived storage to DSS (default False)
      alt_period             -- optional alternate output period in minutes for resampling (default None)
      alt_period_string      -- string label for the alternate period (e.g., '1Day') (default None)
      lookback_padding       -- minutes of lookback padding (currently unused; default 1440)

    Output:
      Returns True on successful completion.
      Writes balance flow (and optionally evap/storage) DSS records to output_dss_file.
      Writes a CSV diagnostic file to shared_dir/<res_name>_balance_flow.csv.
    """
    
    # Validate that all input DSS records use the expected time interval
    check_dss_intervals(inflow_records, balance_period_str, currentAlt)
    check_dss_intervals(outflow_records, balance_period_str, currentAlt)
    check_dss_intervals([stage_record, evap_record], balance_period_str, currentAlt)
    
    # Convert the balance period string to a numeric duration in hours
    balance_period = get_balance_period(balance_period_str) # convert to (float) hours
    print('balance_period ' + str(balance_period))
    
    # Unit conversion factors between cfs (period-average) and acre-ft
    cfs_2_acreft = balance_period * 3600. / 43559.9
    acreft_2_cfs = 1. / cfs_2_acreft

    # Extract the run time window start and end as HEC-formatted strings
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #01Jan2014 0000

    # add lookback padding to enable ResSim to have balance flows on 1st timestep
    # starttime_hectime_obj = HecTime(starttime_str).add(lookback_padding)
    # starttime_str = starttime_hectime_obj.date()

    # Convert start and end time strings to integer HecTime values for bounds-checking
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    
    # Log the compute time window for traceability
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    
    # Open the source DSS file for reading inflow, stage, and evaporation records
    dssFm = HecDss.open(dss_file)

    # Initialize accumulators for inflows, outflows, and their shared time axis
    inflows = []
    outflows = []
    times = []

    # -----------------------------------------------------------------------
    # Read inflows
    # -----------------------------------------------------------------------
    
    
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        print('\nreading: ' + str(pathname))
        try:
            # Read the inflow time series from the DSS file for the run time window
            print(starttime_str, endtime_str)
            print(dss_file)
            tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units
            # print('num values {0}'.format(len(values)))
            # print('start {0}'.format(ts_data.getStartTime()))
            # print('end {0}'.format(ts_data.getEndTime()))
            
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
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)
        
        # Convert SI units (cms) to US customary (cfs) if needed
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting inflow to cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # First record initializes the accumulator; subsequent records are summed element-wise
        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            for vi, v in enumerate(values):
                inflows[vi] += v

    
    
    # -----------------------------------------------------------------------
    # Read outflows
    # -----------------------------------------------------------------------
    
    
    print('Reading outflow records')
    for j, outflow_record in enumerate(outflow_records):  # for each of the dss paths in inflow_records    
        pathname = outflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
        try:
            # Extract values, time axis, and units from the outflow time-series container
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units
            
            # Trim leading values that fall before the run time window start
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
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
            
            # Fatal error: log and abort compute if the outflow record cannot be read
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # Convert SI units (cms) to US customary (cfs) if needed
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting outflow cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # First record initializes the accumulator; subsequent records are summed element-wise
        if len(outflows) == 0:
            outflows = values
        else:
            for vi, v in enumerate(values):
                outflows[vi] += v


    # -----------------------------------------------------------------------
    # Compute net inflow minus outflow for each time step
    # -----------------------------------------------------------------------

    # Inflow minus outflow record
    inflow_outflow = []
    for i in range(len(inflows[1:])):
        inflow_outflow.append(inflows[i+1] - outflows[i+1])
   # this is in cfs (period avg vals)

    # -----------------------------------------------------------------------
    # Read stage (reservoir elevation)
    # -----------------------------------------------------------------------
    
    print('Reading stage')
    tsc = read_ts_rec_w_optional_fname(dssFm, stage_record, starttime_str, endtime_str)
    try:
        # Extract stage values and time axis from the container
        stage = tsc.values
        hectimes = tsc.times
        
        # Trim leading values that fall before the run time window start
        if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
            print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
            st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
            stage = stage[st_offset:]
            hectimes = hectimes[st_offset:]
            
        # Trim trailing values that fall after the run time window end
        if hectimes[-1] > endtime_hectime:
            print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
            st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
            stage = stage[:(len(hectimes) - st_offset)]
            hectimes = hectimes[:(len(hectimes) - st_offset)]
        print('Number Stage Values: {0}'.format(len(stage)))

        # Convert SI elevation units (meters) to US customary (feet) if needed
        if tsc.units.lower() == 'm':
            currentAlt.addComputeMessage('Converting stage m to ft')
            print('Converting stage cms to cfs')
            convvals = []
            for elev in stage:
                convvals.append(elev * 3.280839895)
            stage = convvals
        
    except HecMathException:
        
        # Fatal error: log and abort compute if the stage record cannot be read
        currentAlt.addComputeMessage('ERROR reading' + str(stage_record))
        sys.exit(-1)

    # -----------------------------------------------------------------------
    # Read evaporation depth
    # -----------------------------------------------------------------------
    
    print('Reading evap')
    tsc = read_ts_rec_w_optional_fname(dssFm, evap_record, starttime_str, endtime_str)
    try:
        # Extract evaporation values and time axis from the container
        evap = tsc.values
        hectimes = tsc.times
        
        # Trim leading values that fall before the run time window start
        if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
            print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
            st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
            evap = evap[st_offset:]
            hectimes = hectimes[st_offset:]
            
        # Trim trailing values that fall after the run time window end
        if hectimes[-1] > endtime_hectime:
            print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
            st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
            evap = evap[:(len(hectimes) - st_offset)]
            hectimes = hectimes[:(len(hectimes) - st_offset)]
        
        print('Number Evap Values: {0}'.format(len(evap)))
    
    except HecMathException:
        # Fatal error: log and abort compute if the evaporation record cannot be read
        currentAlt.addComputeMessage('ERROR reading' + str(evap_record))
        sys.exit(-1)

    # -----------------------------------------------------------------------
    # Build conic storage array for interpolation later
    # -----------------------------------------------------------------------

    # Build conic storage array for interpolation later
    # Precomputes cumulative storage at each elevation table point using the conic formula
    conic_storage = build_conic_storage_array(elev_stor_area['elev'], elev_stor_area['area'])

    # -----------------------------------------------------------------------
    # Balance flow calculation loop
    # -----------------------------------------------------------------------

    # Calculations
    
    n = len(stage) - 1
    flow_resid = []
    flow_evap = []
    
    # List to store storage values at the start of each time step (acre-ft)
    storage_record = []

    # Iterate over each consecutive pair of stage values to compute the balance flow
    for k in range(n):
        
        # Retrieve stage at start and end of current time step
        stage_start = stage[k]
        stage_end = stage[k+1]

        # Compute storage at start and end of the time step using chosen interpolation method
        if use_conic:
            # Conic interpolation: find the bounding table index then interpolate
            idx1 = get_elev_layer_idx(elev_stor_area['elev'], stage_start, elev_stor_area)
            storage_start = conic_storage_interp(stage_start, elev_stor_area['elev'], elev_stor_area['area'], conic_storage, idx1)
            idx2 = get_elev_layer_idx(elev_stor_area['elev'], stage_end, elev_stor_area)
            storage_end = conic_storage_interp(stage_end, elev_stor_area['elev'], elev_stor_area['area'], conic_storage, idx2)
        else:
            # Linear interpolation from the elevation-storage table
            storage_start = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], stage_start)
            storage_end = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], stage_end)

        # Change in storage implied by the stage change over this time step (acre-ft)
        delta_stor_from_stage = storage_end - storage_start  # in acre-ft
        
        # Convert storage change from acre-ft to an equivalent average flow rate (cfs)
        delta_stor_flow = delta_stor_from_stage * acreft_2_cfs # in cfs
        
        # Net inflow minus outflow for this time step (cfs, period average)
        inflow_minus_outflow = inflow_outflow[k]  # in cfs
        
        # Compute average surface area over the time step for evaporation loss calculation
        area_avg = 0.5 * (linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_start) +
                          linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_end))
        
        # Evaporation flow loss: evap depth (ft) * area (acre) converted to cfs
        evap_flow_loss = (evap[k] * area_avg) * acreft_2_cfs  # in cfs
        
        # Residual balance flow = stage-implied storage change rate minus (inflow - outflow - evap)
        resid = delta_stor_flow - (inflow_minus_outflow - evap_flow_loss)
        
        # Accumulate results for this time step
        flow_resid.append(resid)
        flow_evap.append(evap_flow_loss)
        storage_record.append(storage_start)

    # -----------------------------------------------------------------------
    # Write diagnostic CSV output
    # -----------------------------------------------------------------------
    
    if True:
        print('Writing to CSV')
        # dump to CSV if DSS is mis-behaving or if flows are -999. etc.
        with open(os.path.join(shared_dir, "{0}_balance_flow.csv".format(res_name)), 'w') as opf:
            opf.write('date, balance_flow [cfs]\n')
            for i in range(len(flow_resid)):
                # Write each time step's HecTime integer and balance flow to the CSV
                new_line = ','.join([str(times[i]), str(flow_resid[i]), '\n'])
                opf.write(new_line)
    
    # Open the output DSS file for writing balance flow (and optionally evap/storage) records    
    dssFm_out = HecDss.open(output_dss_file)
    
    # Output record

	# sometimes ResSim does not include the start record in period average simulations, so if one flow or elevation data
	# record is missing, the calc can sometimes go way off.  Constrain to realistic values, set invalid to zero.
    if flow_resid[0] > 300000.0 or flow_resid[0] < -300000.0:
        flow_resid[0] = 0.0

    # Compute the time step size in HecTime integer units (difference between first two timestamps)
    steptime = times[1]-times[0]
    
    # Build the output TimeSeriesContainer for the balance flow record
    tsc = TimeSeriesContainer()
    
    # Set the start time one step before the first computed balance flow
    # to align with ResSim's period-average convention
    tsc.startTime = times[0] - steptime
    tsc.interval = int(balance_period)*60
    tsc.fullName = output_dss_record_name

    # copy back 1st balance flow record 2 steps, instead of writing from 1st valid balance calc.
    # otherwise, time-averaging the balanece flows later leaves off the 1st time step needed for a ResSim run
    # best we can do I guess to make ResSim computes work
    tsc.values = [flow_resid[0],flow_resid[0]] + flow_resid
    tsc.units = 'CFS'
    tsc.type = 'PER-AVER'
    tsc.numberValues = len(tsc.values)
    
    # Write the balance flow time-series container to the output DSS file
    dssFm_out.write(tsc)

    # If an alternate output period is specified and differs from the balance period, resample
    if alt_period is not None:
        if alt_period_string.lower() != balance_period_str.lower():
            # Read the just-written record back and time-average to the alternate period
            tsm = dssFm_out.read(output_dss_record_name)
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            dssFm_out.write(tsm_new_interval)

    # Optionally write the evaporation flow loss record to DSS
    if write_evap:
        tsc = TimeSeriesContainer()
        tsc.times = times
        tsc.fullName = evap_dss_record_name
        tsc.values = flow_evap
        tsc.startTime = times[1]
        tsc.units = 'CFS'
        tsc.type = 'PER-AVER'
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)

    if write_storage:
        tsc = TimeSeriesContainer()
        tsc.times = times
        tsc.fullName = storage_dss_record_name
        tsc.values = storage_record
        tsc.startTime = times[1]
        tsc.units = "AC-FT"
        tsc.type = 'INST-VAL'  # is this right?
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)

    # Close both the source and output DSS file connections
    dssFm.close()
    dssFm_out.close()
    return True
