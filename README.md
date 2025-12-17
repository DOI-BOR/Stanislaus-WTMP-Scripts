# Stanislaus WTMP Scripts
The HEC-WAT scripts folder for the Stanislaus River within the Bureau of Reclamation (USBR) Water Temperature Modeling Platform (WTMP).
This repository is a dependency of the Stanislaus WTMP Study repository.
These scripts are called at different points in the WTMP workflow.
These scripts are all in Jython, the implementation of Python in Java.

## W2 files
* Pre-Process_W2_Stanislaus.py
### Forecast
### Hindcast
### Planning
### Shared or not specified

## ResSim files
* Pre-Process_ResSim_StanPO.py
* Pre-Process_ResSim_StanPO_.py
### Forecast
### Hindcast
* Acc_Dep_ResSim_Stanislaus.py
### Planning
### Shared or not specified


## Shared files between W2 and ResSim
These files contain various function that are shared by different models.
* create_balance_flow_jython.py
* DMS_preprocess.py
* DSS_Tools.py
* Simple_DSS_Functions.py

## Dependencies
There are dependencies that come from the WAT and from external locations that must be brought in. 

- hec
  - hec.heclib
  - hec.io
  - hec.hecmath
- rma
  - com.rma
  - rma.util.RMAConst

## Usage
### Usage withing the build process
To be added later.

### Post build implementation
After making any desired changes to the code, files must be replaced in the scripts folder in a WTMP Stanislaus Study folder. 
Scripts will be triggered as various points of the modeling workflow.
