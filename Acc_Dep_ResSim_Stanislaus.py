'''
Created on 8/7/2023
@note:
'''

# Import the balance flow computation module and reload to pick up any changes made since initial import
import create_balance_flow_jython as cbfj
reload(cbfj)

#  Import the Project accessor to retrieve the current project's directory and metadata
from com.rma.model import Project
import os
import Simple_DSS_Functions as sdf
reload(sdf)

def computeAlternative(currentAlternative, computeOptions):
    """
    Computes water balance flows for New Melones and Tulloch reservoirs
    in the Stanislaus River system.

    For each reservoir, this function:
      - Resolves DSS input/output file paths and the run time window
      - Resamples daily time-series records to hourly resolution as needed
      - Assembles inflow, outflow, stage, and evaporation record references
      - Reads the elevation-storage-area table from a CSV file
      - Calls create_balance_flows() to compute and write the balance flow DSS records

    Inputs:
      currentAlternative  -- the WAT scripting alternative object; provides the
                             compute time step, name, and message logging interface
      computeOptions      -- WAT compute options object; provides the DSS filename,
                             run time window, and run directory

    Output:
      Returns True on successful completion of both reservoir balance flow computations.
    """
    
    
    # Log the name of the alternative currently being computed
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    
    # Retrieve the primary DSS output file path from the compute options
    dss_file = computeOptions.getDssFilename()
    
    # Retrieve the run time window (start/end dates) for filtering DSS records
    rtw = computeOptions.getRunTimeWindow()

    # New Melones Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    # Retrieve the run directory from compute options (used for locating intermediate files)
    run_dir = computeOptions.getRunDirectory()
    
    # Retrieve the WAT project root directory for building absolute shared-data paths
    project_dir = Project.getCurrentProject().getProjectDirectory()
    
    # Log resolved directories for traceability during compute
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    
    # Retrieve the balance computation time step string from the alternative (e.g., '1Hour')
    balance_period_str = currentAlternative.getTimeStep()
    
    # Build the path to the shared data directory where DSS and CSV files reside
    shared_dir = os.path.join(project_dir, 'shared')

    # Full path to the DMS Stanislaus hydrology time-series DSS file (source hydro data)
    DMS_hydro_dss_file = os.path.join(shared_dir, "DMS_StanislausHydroTS.dss")
    
    # Full path to the pre-process DSS output file (written and read during this compute)
    output_dss_file = os.path.join(shared_dir,'DMS_Stanislaus_ResSim_Pre-Process.dss')
    
    # Full path to the historical fallback DSS file (used when observed data is unavailable)
    fallback_dss_file = os.path.join(shared_dir,'WTMP_Stanislaus_Historical.dss')
    
    # Resample New Melones combined inflow from daily to hourly resolution in the output DSS file
    sdf.resample_dss_ts(output_dss_file,'/MR Stan.-New Melones/Combined Inflow/Flow//1Day/ResSim_PreProcess/',rtw,output_dss_file,'1HOUR')

    # Define the New Melones inflow DSS record (combined inflow resampled to 1-hour)
    inflow_records = ['::'.join([output_dss_file,'/MR Stan.-New Melones/Combined Inflow/Flow//1Hour/ResSim_PreProcess/']),]

    # Define New Melones outflow DSS records: generation and outlet releases at 1-hour resolution
    outflow_records = ['/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/ResSim_PreProcess/',  # not in pre-process file ...
                       '/MR Stan.-New Melones/NML-Outlet Release/Flow//1Hour/240.1.125.3.1/']

    # DSS record path for New Melones reservoir elevation (instantaneous, 1-hour)
    stage_record = '/MR Stan.-New Melones/NML-Elevation/Elev//1Hour/240.1.145.1.1/'
    
    # Use a synthetic zero-flow record as the evaporation input for New Melones
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    # Read the elevation-storage-area lookup table from the New Melones CSV file
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_new_melones.csv'), 'New Melones') #TODO: check this

    # Control flags: disable conic interpolation, evap writing, and storage writing for New Melones
    use_conic = False
    write_evap = False
    write_storage = False

    # DSS record names for derived evaporation, storage, and balance flow outputs
    evap_dss_record_name = "/NEW MELONES/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/NEW MELONES/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/NEW MELONES/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    # Override output record name when conic interpolation is enabled
    if use_conic:
        output_dss_record_name = "/NEW MELONES/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/NEW MELONES/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # Compute and write the New Melones water balance flow to the output DSS file
    # alt_period=1440 minutes corresponds to a 1-day aggregation period for balance computation
    cbfj.create_balance_flows(currentAlternative, rtw, 'New Melones', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')
                                #alt_period=1440*7, alt_period_string='1Week'), was 1week in calibration; hard to make it work in scripting


    # Tulloch Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    # Tulloch inflows are the New Melones generation and outlet releases (reservoir-to-reservoir routing)
    inflow_records = ['/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/ResSim_PreProcess/',  # not in pre-process file ...
                       '/MR Stan.-New Melones/NML-Outlet Release/Flow//1Hour/240.1.125.3.1/']

    # Define Tulloch outflow DSS records: generation release, regulating flow, and spillway release
    outflow_records = ['::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Generation Release/Flow//1HOUR/ResSim_PreProcess/']),
                       '::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1HOUR/241.1.125.4.1/']),
                       '::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Spillway Release/Flow//1HOUR/241.1.125.3.1/'])]

    # Resample Tulloch reservoir elevation from daily to hourly; prepend first value for continuity
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Stan.-Tulloch/TUL-Reservoir Elevation/Elev//1Day/241.1.145.1.1/',rtw,
                        output_dss_file,'1HOUR',prepend_first_value=True,inst_val=True)
    
    # DSS record path for Tulloch reservoir elevation (resampled to 1-hour, written to output DSS)
    stage_record = '::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Reservoir Elevation/Elev//1Hour/241.1.145.1.1/'])
    
    # Use a synthetic zero-flow record as the evaporation input for Tulloch
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    # Read the elevation-storage-area lookup table from the Tulloch CSV file
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_tulloch.csv'), 'Tulloch') #TODO: check this

    # Control flags: disable conic interpolation, evap writing, and storage writing for Tulloch
    use_conic = False
    write_evap = False
    write_storage = False

    # DSS record names for derived evaporation, storage, and balance flow outputs
    evap_dss_record_name = "/TULLOCH/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/TULLOCH/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/TULLOCH/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    # Override output record name when conic interpolation is enabled
    if use_conic:
        output_dss_record_name = "/TULLOCH/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/TULLOCH/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # Compute and write the Tulloch water balance flow to the output DSS file
    # alt_period=1440 minutes corresponds to a 1-day aggregation period for balance computation
    cbfj.create_balance_flows(currentAlternative, rtw, 'Tulloch', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')
                                #alt_period=1440*7, alt_period_string='1Week'), was 1week in calibration; hard to make it work in scripting

    # Both reservoir balance flow computations completed successfully
    return True
