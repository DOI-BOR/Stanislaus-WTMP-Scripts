
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

# Java TimeZone class for timezone-aware time operations
from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

#  Import the accumulation/depletion compute module for the Stanislaus ResSim workflow
# and reload to ensure the latest version is used after sys.path has been corrected
import Acc_Dep_ResSim_Stanislaus
reload(Acc_Dep_ResSim_Stanislaus)

#  Import the DMS pre-processing module for the Stanislaus ResSim workflow
# and reload to ensure the latest version is used
import DMS_preprocess
reload(DMS_preprocess)


def computeAlternative(currentAlternative, computeOptions):
    """
    Entry point for the WAT scripting alternative compute for the ResSim Stanislaus workflow.

    Orchestrates two sequential sub-computes:
      1) DMS pre-processing  -- fixes DMS data types/units and computes derived DSS records
                                needed as ResSim inputs (via DMS_preprocess.preprocess_ResSim_Stanislaus)
      2) Accumulation/Depletion compute -- runs the Acc/Dep ResSim Stanislaus calculation
                                           (via Acc_Dep_ResSim_Stanislaus.computeAlternative)

    Both sub-computes must return True for the overall alternative compute to succeed.

    Inputs:
      currentAlternative -- WAT scripting alternative object providing the alternative name,
                            time step, and compute message logging interface
      computeOptions     -- WAT compute options object providing the DSS filename,
                            run time window, and run directory

    Output:
      Returns True if both sub-computes succeed.
      Returns None implicitly if either sub-compute returns a falsy value.
    """
    
    # Log the name of the alternative currently being computed
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # Run the DMS pre-processing step to prepare input DSS records
    data_preprocess = DMS_preprocess.preprocess_ResSim_Stanislaus(currentAlternative, computeOptions)

    # Run the accumulation/depletion compute for the Stanislaus ResSim alternative
    acc_dep = Acc_Dep_ResSim_Stanislaus.computeAlternative(currentAlternative, computeOptions)

    # Return True only if both sub-computes completed successfully
    if data_preprocess and acc_dep:
        return True

