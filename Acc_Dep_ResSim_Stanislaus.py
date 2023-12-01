'''
Created on 8/7/2023
@author: Scott Burdick-Yahya, and Ben Saenz
@organization: Resource Management Associates
@contact: scott@rmanet.com
@note:
'''

import create_balance_flow_jython as cbfj
reload(cbfj)
from com.rma.model import Project
import os
import Simple_DSS_Functions as sdf
reload(sdf)

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    # New Melones Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period_str = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    DMS_hydro_dss_file = os.path.join(shared_dir, "DMS_StanislausHydroTS.dss")
    output_dss_file = os.path.join(shared_dir,'DMS_Stanislaus_ResSim_Pre-Process.dss')
    fallback_dss_file = os.path.join(shared_dir,'WTMP_Stanislaus_Historical.dss')

    sdf.resample_dss_ts(output_dss_file,'/MR Stan.-New Melones/Combined Inflow/Flow//1Day/ResSim_PreProcess/',rtw,output_dss_file,'1HOUR')

    inflow_records = ['::'.join([output_dss_file,'/MR Stan.-New Melones/Combined Inflow/Flow//1Hour/ResSim_PreProcess/']),]

    outflow_records = ['/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/ResSim_PreProcess/',  # not in pre-process file ...
                       '/MR Stan.-New Melones/NML-Outlet Release/Flow//1Hour/240.1.125.3.1/']

    stage_record = '/MR Stan.-New Melones/NML-Elevation/Elev//1Hour/240.1.145.1.1/'
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_new_melones.csv'), 'New Melones') #TODO: check this

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/NEW MELONES/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/NEW MELONES/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/NEW MELONES/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/NEW MELONES/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/NEW MELONES/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

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

    inflow_records = ['/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/ResSim_PreProcess/',  # not in pre-process file ...
                       '/MR Stan.-New Melones/NML-Outlet Release/Flow//1Hour/240.1.125.3.1/']

    outflow_records = ['::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Generation Release/Flow//1HOUR/ResSim_PreProcess/']),
                       '::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1HOUR/241.1.125.4.1/']),
                       '::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Spillway Release/Flow//1HOUR/241.1.125.3.1/'])]

    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Stan.-Tulloch/TUL-Reservoir Elevation/Elev//1Day/241.1.145.1.1/',rtw,
                        output_dss_file,'1HOUR',prepend_first_value=True,inst_val=True)
    stage_record = '::'.join([output_dss_file,'/MR Stan.-Tulloch/TUL-Reservoir Elevation/Elev//1Hour/241.1.145.1.1/'])
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_tulloch.csv'), 'Tulloch') #TODO: check this

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/TULLOCH/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/TULLOCH/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/TULLOCH/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/TULLOCH/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/TULLOCH/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'Tulloch', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')
                                #alt_period=1440*7, alt_period_string='1Week'), was 1week in calibration; hard to make it work in scripting


    return True
