
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
#import hec.hecmath.TimeSeriesMath as tsmath

sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

import DSS_Tools
reload(DSS_Tools)
import Simple_DSS_Functions as sdf
reload(sdf)

units_need_fixing = ['tenths','m/s','deg',] #'radians',]

def fix_DMS_types_units(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    dss = HecDss.open(dss_file)
    recs = dss.getPathnameList()
    for r in recs:
        tsm = dss.read(r)
        rlow = r.lower()
        if "/flow" in rlow or "/1day/" in rlow:
            if not "/elev" in rlow:
                tsm.setType('PER-AVER')
                dss.write(tsm)
        
        if tsm.getUnits() in units_need_fixing:
            if tsm.getUnits() == 'tenths':
                # save off a copy of cloud record in 0-1 for ResSim
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-FRAC'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'FRAC'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 10.0                
                dss.write(tsc)
            if tsm.getUnits() == 'radians':
                # save off a copy in deg
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-DEG'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'deg'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0                
                dss.write(tsc)
            if tsm.getUnits() == 'deg':
                # save off a copy in redians
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-RADIANS'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'radians'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)                
                dss.write(tsc)
            if tsm.getUnits() == 'm/s':
                # make a copy divied by kph conversion as a hack to get W2 linking the wind speed correctly 
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                if not "w2link" in rec_parts[3].lower():
                    rec_parts[3] += '-W2link'
                    tsc.fullName = '/'.join(rec_parts)
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    dss.write(tsc)
    dss.close()

def DMS_fix_units_types(hydro_dss,met_dss_file):
    fix_DMS_types_units(hydro_dss)
    fix_DMS_types_units(met_dss_file)

def compute_new_melones_flows(currentAlternative, rtw, hydro_dss, output_dss_file):
    # Sum Stanislaus tribs to get total New Melones inflow
    inflow_records = ['/MR Stan.-New Melones/11293200 MF Stan. R BL Sanbar-Flow/Flow//1Day/240.62.125.1.1/',
                      '/MR Stan.-New Melones/11295250 Colliverville PP NR Murphys-Flow/Flow//1Day/240.6.125.1.1/',
                      '/MR Stan.-New Melones/11295300 NF Stan. R BL Beaver Creek-Flow/Flow//1Day/240.65.125.1.1/',
                      '/MR Stan.-New Melones/11295505 Stan. PP NR Hathaway Pines-Flow/Flow//1Day/240.67.125.1.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Stan.-New Melones/Combined Inflow/Flow//1Day/ResSim_PreProcess/', output_dss_file)

	# add min 4 cfs to generation flow from New Melones (to prevent zero flow)
    DSS_Tools.min_ts(hydro_dss,'/MR Stan.-New Melones/NML-Generation Release/Flow//1Hour/240.1.125.1.1/', 4.0, output_dss_file, 'ResSim_PreProcess')

def compute_tulloch_stan_flows(currentAlternative, rtw, hydro_dss, output_dss_file):

	# Sum Tulloch outflows for Goodwin balance	
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Tulloch/TUL-Generation Release/Flow//1Day/241.1.125.2.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1Day/241.1.125.4.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Tulloch/TUL-Spillway Release/Flow//1Day/241.1.125.3.1/',rtw,output_dss_file,'1HOUR')
    outflow_records = ['/MR Stan.-Tulloch/TUL-Generation Release/Flow//1HOUR/241.1.125.2.1/',
                       '/MR Stan.-Tulloch/TUL-Ctrl Regulating Flow/Flow//1HOUR/241.1.125.4.1/',
                       '/MR Stan.-Tulloch/TUL-Spillway Release/Flow//1HOUR/241.1.125.3.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, outflow_records, output_dss_file,
              '/MR Stan.-Tulloch/Combined Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Goodwin Dam balance Flow
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Release to River/Flow//1Day/242.1.125.3.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Joint Canal Diversion/Flow//1Day/242.1.125.6.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-South Canal Diversion/Flow//1Day/242.1.125.5.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Spillway Flow/Flow//1Day/242.1.125.2.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(hydro_dss,'/MR Stan.-Goodwin/GDW-Outlet Flow/Flow//1Day/242.1.125.1.1/',rtw,output_dss_file,'1HOUR')
    flow_records = ['/MR Stan.-Goodwin/GDW-Release to River/Flow//1HOUR/242.1.125.3.1/',
	                '/MR Stan.-Goodwin/GDW-Joint Canal Diversion/Flow//1HOUR/242.1.125.6.1/',
	                '/MR Stan.-Goodwin/GDW-South Canal Diversion/Flow//1HOUR/242.1.125.5.1/',
                    '/MR Stan.-Tulloch/Combined Outflow/Flow//1Hour/ResSim_PreProcess/',]  
    out_rec = '/MR Stan.-Goodwin/Goodwin Dam Balance Flow/Flow//1HOUR/ResSim_PreProcess/'
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_records, output_dss_file, 
                                    [None,True,True,False], out_rec, output_dss_file)

	# Ripon Balance - shift goodwin release 1 day back before balance to account for some travel time
	#ripon_tsm = ripon_tsm_shift.shiftInTime(int timeShiftMinutes)

def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):
    pass
    

def preprocess_W2_Stanislaus(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_Stanislaus_ResSim_Pre-Process.dss')

    hydro_dss = os.path.join(shared_dir, 'DMS_StanislausHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_StanislausMet.dss')
    fix_DMS_types_units(met_dss_file)

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')

    compute_new_melones_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    return True


def preprocess_ResSim_Stanislaus(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_Stanislaus_ResSim_Pre-Process.dss')

    hydro_dss = os.path.join(shared_dir, 'DMS_StanislausHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_StanislausMet.dss')
    fix_DMS_types_units(met_dss_file)

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')

    compute_new_melones_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_tulloch_stan_flows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    return True

