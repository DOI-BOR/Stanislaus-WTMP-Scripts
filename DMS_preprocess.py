
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
import os, sys

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

import DSS_Tools
reload(DSS_Tools)
import Simple_DSS_Functions as sdf
reload(sdf)

#  Units that require correction when reading from the DMS hydrology/met DSS files
units_need_fixing = ['tenths','m/s','deg','kph'] #'radians',]

def fix_DMS_types_units(dss_file):
    """
    Corrects DSS record data types and unit conversions for records sourced from the DMS.

    For each record in the DSS file:
      - Flow and 1-day records (excluding elevation) are set to PER-AVER data type.
      - Records with units in the 'units_need_fixing' list are converted or duplicated
        with corrected units, as follows:
          * 'tenths'  -> fractional (0-1) copy written with units 'FRAC' and '-FRAC' appended to parameter name
          * 'radians' -> degree copy written with units 'deg' and '-DEG' appended to parameter name
          * 'deg'     -> radian copy written with units 'radians' and '-RADIANS' appended to parameter name
          * 'kph'     -> converted in-place to m/s; a W2-link copy (further divided by 3.6) is also written
          * 'm/s'     -> a W2-link copy (divided by 3.6) is written with '-W2link' appended to parameter name

    Inputs:
      dss_file -- full path to the DSS file to be corrected

    Output:
      No return value. Modifies records in-place in the DSS file.
    """
    
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    
    # Open the DSS file for reading and writing
    dss = HecDss.open(dss_file)
    
    # Retrieve all DSS pathnames present in the file
    recs = dss.getPathnameList()
    
    # Iterate over each record pathname to check and correct type/units
    for r in recs:
        tsm = dss.read(r)
        rlow = r.lower()
        
        # Set flow and 1-day records (excluding elevation) to PER-AVER data type
        if "/flow" in rlow or "/1day/" in rlow:
            if not "/elev" in rlow:
                tsm.setType('PER-AVER')
                dss.write(tsm)
        
        # Check if the record's units require conversion or duplication
        if tsm.getUnits() in units_need_fixing:
            if tsm.getUnits() == 'tenths':
                # save off a copy of cloud record in 0-1 for ResSim
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-FRAC'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'FRAC'
                
                # Divide each value by 10 to convert tenths to a 0-1 fraction
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 10.0                
                dss.write(tsc)
                
            # 'radians': write a degree copy with '-DEG' appended to the parameter name    
            if tsm.getUnits() == 'radians':
                # save off a copy in deg
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-DEG'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'deg'
                
                # Convert radians to degrees: degrees = radians / (2*pi) * 360
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0  
                    
                dss.write(tsc)
                
            # 'deg': write a radian copy with '-RADIANS' appended to the parameter name    
            if tsm.getUnits() == 'deg':
                # save off a copy in redians
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-RADIANS'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'radians'
                
                # Convert degrees to radians: radians = degrees / 360 * (2*pi)
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)   
                    
                dss.write(tsc)
                
            # 'kph': convert in-place to m/s, then also write a W2-link copy (divided by 3.6 again)    
            if tsm.getUnits() == 'kph':
                # convert to m/s 
                tsc = tsm.getData()
                tsc.units = 'm/s'
                
                # Convert kph to m/s by dividing by 3.6
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 3.6
                dss.write(tsc)

                # also, add w2link
                rec_parts = tsc.fullName.split('/')
                # Only append '-W2link' if it hasn't already been added
                if not "w2link" in rec_parts[3].lower():
                    rec_parts[3] += '-W2link'
                    tsc.fullName = '/'.join(rec_parts)
                    
                    # Divide by 3.6 again to produce the W2-link scaling (already in m/s, divide further)
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    dss.write(tsc)
            
            # 'm/s': write a W2-link copy (divided by 3.6) with '-W2link' appended to the parameter name            
            if tsm.getUnits() == 'm/s':
                # make a copy divied by kph conversion as a hack to get W2 linking the wind speed correctly 
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                
                # Only append '-W2link' if it hasn't already been added
                if not "w2link" in rec_parts[3].lower():
                    rec_parts[3] += '-W2link'
                    tsc.fullName = '/'.join(rec_parts)
                    
                    # Divide by 3.6 to produce the W2-link scaled copy
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    dss.write(tsc)
    
    # Close the DSS file after all corrections are written                
    dss.close()

def DMS_fix_units_types(hydro_dss,met_dss_file):
    """
    Applies fix_DMS_types_units to both the hydrology and meteorology DSS files.

    Inputs:
      hydro_dss    -- full path to the DMS hydrology DSS file
      met_dss_file -- full path to the DMS meteorology DSS file

    Output:
      No return value. Both DSS files are corrected in-place.
    """
    
    # Fix data types and units in the hydrology DSS file
    fix_DMS_types_units(hydro_dss)
    
    # Fix data types and units in the meteorology DSS file
    fix_DMS_types_units(met_dss_file)

def compute_new_melones_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    """
    Computes combined inflow to New Melones reservoir by summing four tributary gauge records,
    and applies a minimum flow constraint to the New Melones generation release.

    Inputs:
      currentAlternative -- WAT scripting alternative object for logging and context
      rtw                -- WAT run time window object providing start/end time strings
      hydro_dss          -- full path to the DMS hydrology DSS file (source records)
      output_dss_file    -- full path to the pre-process DSS file where results are written

    Output:
      No return value. Writes the following DSS records to output_dss_file:
        - Combined inflow (sum of four tributaries) at 1-day resolution
        - New Melones generation release with a 4 cfs minimum applied
    """
    
    # Sum Stanislaus tribs to get total New Melones inflow
    inflow_records = ['/MR Stan.-New Melones/11293200 MF Stan. R BL Sanbar-Flow/Flow//1Day/240.62.125.1.1/',
                      '/MR Stan.-New Melones/11295250 Colliverville PP NR Murphys-Flow/Flow//1Day/240.6.125.1.1/',
                      '/MR Stan.-New Melones/11295300 NF Stan. R BL Beaver Creek-Flow/Flow//1Day/240.65.125.1.1/',
                      '/MR Stan.-New Melones/11295505 Stan. PP NR Hathaway Pines-Flow/Flow//1Day/240.67.125.1.1/']
    
    # Sum the four tributary inflow records and write the combined inflow to the output DSS file
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Stan.-New Melones/Combined Inflow/Flow//1Day/ResSim_PreProcess/', output_dss_file)

	# add min 4 cfs to generation flow from New Melones (to prevent zero flow)
    DSS_Tools.min_ts_flow_cfs(hydro_dss,'/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/240.1.125.1.1/', 4.0, output_dss_file, 'ResSim_PreProcess')

def compute_tulloch_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    
    """
    Computes and writes several derived flow records for Tulloch reservoir and Goodwin Dam,
    including combined outflows, hourly resampled releases, and W2/ResSim pre-process records.

    Inputs:
      currentAlternative -- WAT scripting alternative object for logging and context
      rtw                -- WAT run time window object providing start/end time strings
      hydro_dss          -- full path to the DMS hydrology DSS file (source records)
      output_dss_file    -- full path to the pre-process DSS file where results are written

    Output:
      No return value. Writes the following DSS records to output_dss_file:
        - New Melones combined outflow (gen + outlet release) at 1-hour for W2 pre-process
        - Tulloch generation release with a 4 cfs minimum, resampled to 1-hour
        - Tulloch control regulating flow resampled to 1-hour
        - Tulloch spillway release resampled to 1-hour
        - Tulloch control-regulating-plus-generation combined flow at 1-hour for W2
        - Tulloch combined outflow (gen + ctrl-reg + spillway) at 1-hour for ResSim pre-process
    """

	# Sum stanisalus outflows for independent W2 Tulloch simulatinos
    outflow_records = ['/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/ResSim_PreProcess/',  # not in pre-process file ...
                       '/MR Stan.-New Melones/NML-Outlet Release/Flow//1Hour/240.1.125.3.1/']
    
    # Sum New Melones generation and outlet releases for W2 Tulloch inflow
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, hydro_dss,
              '/MR Stan.-New Melones/Combined Outflow/Flow//1Hour/W2_PreProcess/', output_dss_file)	
	
    # add min 4 cfs to generation flow from New Melones (to prevent zero flow)
    DSS_Tools.min_ts_flow_cfs(hydro_dss,'/MR Stan.-Tulloch/TUL-Generation Release/Flow//1Day/241.1.125.2.1/', 4.0, output_dss_file, 'ResSim_PreProcess')

    # Resample Tulloch generation release, control regulating flow, and spillway release from daily to hourly
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Tulloch/TUL-Generation Release/Flow//1Day/ResSim_PreProcess/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1Day/241.1.125.4.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Tulloch/TUL-Spillway Release/Flow//1Day/241.1.125.3.1/',rtw,output_dss_file,'1HOUR')

	# combine Tulloch cnt-reg and gen flows, for W2
    outflow_records = ['/MR Stan.-Tulloch/TUL-Generation Release/Flow//1HOUR/ResSim_PreProcess/',
                       '/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1HOUR/241.1.125.4.1/']
    
    # Sum Tulloch control-regulating and generation flows for W2 pre-process
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, output_dss_file,
              '/MR Stan.-Tulloch/TUL-Ctrl-Reg-plus-Gen/Flow//1Hour/W2_PreProcess/', output_dss_file)	

    # Sum Tulloch outflows for Goodwin balance	    
    outflow_records = ['/MR Stan.-Tulloch/TUL-Generation Release/Flow//1HOUR/ResSim_PreProcess/',
                       '/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1HOUR/241.1.125.4.1/',
                       '/MR Stan.-Tulloch/TUL-Spillway Release/Flow//1HOUR/241.1.125.3.1/']
    
    # Sum all three Tulloch outflow components for the Goodwin Dam balance flow calculation
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, output_dss_file,
              '/MR Stan.-Tulloch/Combined Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

def compute_stanislaus_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    """
    Computes the Goodwin Dam balance flow by resampling daily Goodwin release components
    to hourly and combining them with Tulloch combined outflow using add/subtract logic.

    Inputs:
      currentAlternative -- WAT scripting alternative object for logging and context
      rtw                -- WAT run time window object providing start/end time strings
      hydro_dss          -- full path to the DMS hydrology DSS file (source records)
      output_dss_file    -- full path to the pre-process DSS file where results are written

    Output:
      No return value. Writes the following DSS records to output_dss_file:
        - Goodwin release, joint/south canal diversions, spillway, and outlet flows resampled to 1-hour
        - Goodwin Dam balance flow (release minus canal diversions plus Tulloch combined outflow)
    """
    
    # Goodwin Dam balance Flow
    # Resample all five Goodwin Dam release and diversion components from daily to hourly resolution
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Release to River/Flow//1Day/242.1.125.3.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Joint Canal Diversion/Flow//1Day/242.1.125.6.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-South Canal Diversion/Flow//1Day/242.1.125.5.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Spillway Flow/Flow//1Day/242.1.125.2.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Outlet Flow/Flow//1Day/242.1.125.1.1/',rtw,output_dss_file,'1HOUR')
    
    # Define the Goodwin flow records to combine: release is added, canal diversions are subtracted,
    # and Tulloch combined outflow is added (None = add, True = subtract, False = add)
    flow_records = ['/MR Stan.-Goodwin/GDW-Release to River/Flow//1HOUR/242.1.125.3.1/',
	                '/MR Stan.-Goodwin/GDW-Joint Canal Diversion/Flow//1HOUR/242.1.125.6.1/',
	                '/MR Stan.-Goodwin/GDW-South Canal Diversion/Flow//1HOUR/242.1.125.5.1/',
                    '/MR Stan.-Tulloch/Combined Outflow/Flow//1Hour/ResSim_PreProcess/',]  
    out_rec = '/MR Stan.-Goodwin/Goodwin Dam Balance Flow/Flow//1HOUR/ResSim_PreProcess/'
    
    # Compute Goodwin balance flow: release - joint diversion - south diversion + Tulloch outflow
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_records, output_dss_file, 
                                    [None,True,True,False], out_rec, output_dss_file)

	# Ripon Balance - shift goodwin release 1 day back before balance to account for some travel time
	#ripon_tsm = ripon_tsm_shift.shiftInTime(int timeShiftMinutes)

def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):
    """
    Placeholder function for computing additional DSS records used for plotting purposes.

    Inputs:
      currentAlternative -- WAT scripting alternative object for logging and context
      rtw                -- WAT run time window object providing start/end time strings
      hydro_dss          -- full path to the DMS hydrology DSS file (source records)
      output_dss_file    -- full path to the pre-process DSS file where results would be written

    Output:
      No return value. Currently a no-op (pass).
    """
    
    pass
    

def preprocess_W2_Stanislaus(currentAlternative, computeOptions):
    """
    Orchestrates the pre-processing workflow for the CE-QUAL-W2 (W2) Stanislaus model.

    Steps:
      1) Resolves file paths for the shared DSS files and project directory
      2) Fixes DMS data types and units in both hydrology and meteorology DSS files
      3) Creates constant DSS records (tiny flow, tens temperature, and zero placeholders)
      4) Computes New Melones inflow and release records
      5) Computes Tulloch outflow and derived records
      6) Computes any additional plotting records (currently a no-op)

    Inputs:
      currentAlternative -- WAT scripting alternative object providing time step,
                            compute messages, and context
      computeOptions     -- WAT compute options object providing DSS filename,
                            run time window, and run directory

    Output:
      Returns True on successful completion.
      Writes pre-processed DSS records to the shared pre-process DSS file.
    """
    
    # Retrieve the primary DSS output file path and run time window from compute options
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    # Retrieve run and project directories for resolving shared file paths
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    
    # Log resolved directories for traceability during compute
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    
    # Retrieve the balance computation time step string from the alternative
    balance_period = currentAlternative.getTimeStep()
    
    # Build the path to the shared data directory
    shared_dir = os.path.join(project_dir, 'shared')

    # Full path to the pre-process DSS output file
    output_dss_file = os.path.join(shared_dir,'DMS_Stanislaus_ResSim_Pre-Process.dss')

    # Full path to the DMS hydrology time-series DSS file
    hydro_dss = os.path.join(shared_dir, 'DMS_StanislausHydroTS.dss')
    
    # Fix data types and units in the hydrology DSS file before processing
    fix_DMS_types_units(hydro_dss)
    
    # Full path to the DMS meteorology DSS file
    met_dss_file = os.path.join(shared_dir,'DMS_StanislausMet.dss')
    
    # Fix data types and units in the meteorology DSS file before processing
    fix_DMS_types_units(met_dss_file)

    # Create a tiny constant flow record at daily resolution (used as a near-zero placeholder)
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    
    # Create a tiny constant flow record at hourly resolution
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    
    # Create a constant 10-degree water temperature record at daily resolution
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    
    # Create zero-flow placeholder records at daily and hourly resolution
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
                        
    # Compute and write New Melones combined inflow and generation release records
    compute_new_melones_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    
    # Compute and write Tulloch outflow and derived records
    compute_tulloch_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    
    # Compute additional plotting records (currently a no-op placeholder)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    return True


def preprocess_ResSim_Stanislaus(currentAlternative, computeOptions):
    """
    Orchestrates the pre-processing workflow for the HEC-ResSim Stanislaus model.

    Steps:
      1) Resolves file paths for the shared DSS files and project directory
      2) Fixes DMS data types and units in both hydrology and meteorology DSS files
      3) Creates constant DSS records (zero flow, zero temperature, gate open/closed placeholders)
      4) Computes New Melones inflow and release records
      5) Computes Tulloch outflow and derived records
      6) Computes Goodwin Dam balance flow records
      7) Computes any additional plotting records (currently a no-op)

    Inputs:
      currentAlternative -- WAT scripting alternative object providing time step,
                            compute messages, and context
      computeOptions     -- WAT compute options object providing DSS filename,
                            run time window, and run directory

    Output:
      Returns True on successful completion.
      Writes pre-processed DSS records to the shared pre-process DSS file.
    """
    # Retrieve the primary DSS output file path and run time window from compute options
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    # Retrieve run and project directories for resolving shared file paths
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    
    # Log resolved directories for traceability during compute
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    
    # Retrieve the balance computation time step string from the alternative
    balance_period = currentAlternative.getTimeStep()
    
    # Build the path to the shared data directory
    shared_dir = os.path.join(project_dir, 'shared')

    # Full path to the pre-process DSS output file
    output_dss_file = os.path.join(shared_dir,'DMS_Stanislaus_ResSim_Pre-Process.dss')

    # Full path to the DMS hydrology time-series DSS file
    hydro_dss = os.path.join(shared_dir, 'DMS_StanislausHydroTS.dss')
    
    # Fix data types and units in the hydrology DSS file before processing
    fix_DMS_types_units(hydro_dss)
    
    # Full path to the DMS meteorology DSS file
    met_dss_file = os.path.join(shared_dir,'DMS_StanislausMet.dss')
    
    # Fix data types and units in the meteorology DSS file before processing
    fix_DMS_types_units(met_dss_file)

    # Create zero-flow placeholder records at daily and hourly resolution
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    
    # Create zero water temperature placeholder records at daily and hourly resolution
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    
    # Create a gate-closed (0) constant record at hourly resolution (INST-VAL type)
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    
    # Create a gate-open (1) constant record at hourly resolution (INST-VAL type)
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')

    # Compute and write New Melones combined inflow and generation release records
    compute_new_melones_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    
    # Compute and write Tulloch outflow and derived records
    compute_tulloch_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    
    # Compute and write Goodwin Dam balance flow and Stanislaus river routing records
    compute_stanislaus_flows(currentAlternative, rtw, hydro_dss, output_dss_file)    
    
    # Compute additional plotting records (currently a no-op placeholder)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    return True

