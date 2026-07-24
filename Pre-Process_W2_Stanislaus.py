
#from com.rma.io import DssFileManagerImpl
import os,time,sys#,shutil
from distutils.dir_util import copy_tree
from com.rma.model import Project

# -----------------------------------------------------------------------
# sys.path cleanup: remove project-specific script folders from other
# watersheds to prevent module name collisions across WAT alternatives
# -----------------------------------------------------------------------

# print current path
print("Current paths: ", sys.path)

# create list of unwanted folders in sys.path
# These watershed-specific folder names are removed to avoid importing
# modules from a different project's scripts directory
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]

# initialize and search for unwanted paths
matching_paths = []
for p in sys.path:
    # Collect any sys.path entry that contains one of the unwanted watershed phrases
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")
for path in matching_paths:
    print(path)

# remove matching paths from sys.path
for path in matching_paths:
    # Only remove if still present (guards against duplicates already removed)
    if path in sys.path:
        sys.path.remove(path)

# append path
# Add the current project's scripts directory so the correct watershed modules are found
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

# Import the DMS pre-processing module for the W2 Stanislaus workflow
# and reload to ensure the latest version is used after sys.path has been corrected
import DMS_preprocess
reload(DMS_preprocess)

#List of CE-QUAL-W2 model alternative names whose input files will be configured
# before the W2 plugin executes the simulations
W2_models_for_input_copy = ['W2 New Melones Prescribed','W2 Tulloch Prescribed']

def backdate_W2_files_to_skip_compute(run_dir):
    """
    Sets the modification timestamps of all W2 model files to 7 days in the past,
    causing the WAT W2 plugin to believe the outputs are already up-to-date and
    skip re-running the simulation. Useful for debugging or bypassing long W2 runs.

    Inputs:
      run_dir -- full path to the WAT run directory, used to resolve the study directory

    Output:
      No return value. Modifies file timestamps in-place under the study's cequal-w2 folder.
    """
    
    # Resolve the study directory from the run directory path
    study_dir = study_dir_from_run_dir(run_dir)
    
    # Compute the target modification time as 7 days before the current time
    current_time = time.time()
    modification_time = current_time - 3600*24*7  # Subtract 7 days (in seconds)
    
    # Walk all files under the cequal-w2 study folder and update their modification times
    for root, dirs, files in os.walk(os.path.join(study_dir,'cequal-w2')):
        for file in files:
            # Skip hidden files (e.g., .DS_Store on macOS)
            if not file.startswith('.'):
                file_path = os.path.join(root,file)
                print("Changing modified time: ",file_path)
                # Preserve the original creation time; only change the modification time
                creation_time = os.path.getctime(file_path)
                os.utime(file_path,(creation_time,modification_time))

def study_dir_from_run_dir(run_dir):
    """
    Resolves the WAT study root directory from a given WAT run directory.

    The expected directory hierarchy is:
      study_dir / runs_dir / w2sim / run_dir

    Inputs:
      run_dir -- full path to the WAT run directory

    Output:
      Returns the full path to the study root directory (three levels up from run_dir).
    """
    
    # Navigate up three directory levels: run_dir -> w2sim -> runs_dir -> study_dir
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def annual_config_dirs_from_run_dir(run_dir,model_name,startyear_str):
    """
    Resolves the three key directory paths needed to configure a W2 model alternative
    for a specific simulation year.

    The three directories are:
      model_dir         -- the W2 model alternative's active input file directory
      annual_config_dir -- the year-specific W2 input configuration directory
      base_dir          -- the backup directory where original model files are preserved

    Inputs:
      run_dir       -- full path to the WAT run directory
      model_name    -- W2 model alternative name (e.g., 'W2 New Melones Prescribed')
      startyear_str -- four-character string of the simulation start year (e.g., '2020')

    Output:
      Returns a tuple (model_dir, annual_config_dir, base_dir) of full directory path strings.
    """
    
    # Resolve the study root directory from the run directory
    study_dir = study_dir_from_run_dir(run_dir)
    
    # Strip 'W2' and 'Prescribed' tokens from the model name to get the bare reservoir name
    model_name_no_w2_prescribed = model_name.replace('W2','').replace('Prescribed','').strip()    
    
    # Active W2 model alternative input directory under the study's cequal-w2 folder
    model_dir = os.path.join(study_dir,'cequal-w2',model_name_no_w2_prescribed,model_name)  
    
    # Year-specific W2 configuration directory containing pre-built annual input files
    annual_config_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,startyear_str)
    
    # Backup directory for preserving original model alternative files before overwriting
    base_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,'base')
    
    return model_dir,annual_config_dir,base_dir

def computeAlternative(currentAlternative, computeOptions):
    """
    Entry point for the WAT scripting alternative compute for the CE-QUAL-W2 Stanislaus workflow.

    For each W2 model listed in W2_models_for_input_copy, this function:
      1) Validates that the simulation spans a single calendar year (W2 requires single-year runs)
      2) Resolves the active model, annual config, and base backup directories
      3) Backs up original model input files to the base directory on first run
      4) Replaces all model input files (except the .w2Alt file) with year-specific config files
      5) Runs the DMS pre-processing step to prepare W2 input DSS records

    Inputs:
      currentAlternative -- WAT scripting alternative object providing the alternative name,
                            time step, and compute message logging interface
      computeOptions     -- WAT compute options object providing the DSS filename,
                            run time window, and run directory

    Output:
      No explicit return value. Side effects include:
        - W2 model alternative directories configured with year-specific input files
        - DMS pre-processing DSS records written to the shared pre-process DSS file
      NOTE: contains a pre-existing tab-indented block at the end of the function
      (the DMS_preprocess call and the commented-out backdate call) which may cause
      an IndentationError in Python 3 due to mixed tabs and spaces.
    """
    
    # Log the name of the alternative currently being computed
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    
    # Retrieve the run time window and extract start and end time strings
    rtw = computeOptions.getRunTimeWindow()
    
    # Extract the four-digit start year from the HEC time string (characters 5-8)
    starttime_str = rtw.getStartTimeString()
    startyear_str = starttime_str[5:9]
    currentAlternative.addComputeMessage('Found start year for W2 simulations:'+startyear_str)
    
    # Extract the four-digit end year for cross-year validation
    endtime_str = rtw.getEndTimeString()
    endyear_str = endtime_str[5:9]
    
    # Warn if the simulation spans more than one calendar year, as W2 requires single-year runs
    if startyear_str != endyear_str:
        currentAlternative.addComputeMessage('WARNING: Start year ({0}) is different from end year ({1}); W2 simulations will likely fail.'.format(starttime_str, endtime_str))    

    # Retrieve the run directory path for resolving study and model directories
    run_dir = computeOptions.getRunDirectory()
    
    # Iterate over each W2 model alternative that requires annual input file configuration
    for W2_model in W2_models_for_input_copy:
        
        # Resolve the active model, annual config, and base backup directories for this model and year
        model_dir,annual_config_dir,base_dir = annual_config_dirs_from_run_dir(run_dir,W2_model,startyear_str)
        
        # Log resolved directory paths for traceability
        currentAlternative.addComputeMessage('model_dir: '+model_dir)
        currentAlternative.addComputeMessage('annual_config_dir: '+annual_config_dir)
        currentAlternative.addComputeMessage('base_dir: '+base_dir)
        
        # Warn if no annual config exists for this model and year; W2 will likely be misconfigured
        if not os.path.exists(annual_config_dir):
            currentAlternative.addComputeMessage(W2_model+'- annual config not found; W2 may be configured incorrectly for this time window.')
        else:        
            # copy original W2 model alternative files to 'base' directory for safekeeping/later returning

            # Create the base backup directory if it does not yet exist
            if not os.path.exists(base_dir):
                os.mkdir(base_dir)
            
            # Check whether the base directory already contains backup files
            base_files = os.listdir(base_dir)
            
            # If the base directory is empty, copy the current model files there as a backup
            if len(base_files) == 0:
                currentAlternative.addComputeMessage(W2_model+'- base files not found; copying from model folder')
                copy_tree(model_dir,base_dir)

            # remove all W2 model input files EXCEPT the .w2alt file
            # Delete all model input files except the WAT alternative descriptor (.w2Alt and its backup)
            for mfile in os.listdir(model_dir):
                if not mfile.endswith('.w2Alt') and not mfile.endswith('.w2Alt.bak'):
                    os.remove(os.path.join(model_dir,mfile))

            # copy over annual config input files 
            # Replace deleted files with the year-specific W2 input configuration files
            copy_tree(annual_config_dir,model_dir)
            currentAlternative.addComputeMessage('Copied W2 inputs file for '+startyear_str+' to '+W2_model+' model alternative folder')
            
            # now, W2 model alternative directory is configured for startyear, ready for the W2 plugin to
            # work it's magic on the W2_con file and execute simulation    

    # Run the DMS pre-processing step to prepare W2 input DSS records for the Stanislaus watershed
    DMS_preprocess.preprocess_W2_Stanislaus(currentAlternative, computeOptions)

    # remove me most of the time
    #backdate_W2_files_to_skip_compute(run_dir)
