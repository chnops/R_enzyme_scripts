R_enzyme_scripts
================

These are scripts that assist with the processing of extracellular enzyme data

|_enzyme_processing_functions.R|
================================

This script includes functions that should be run first.  These are functions with equations that are in a seperate script
for ease of use.

|_enzyme_data_processing.R_|
============================

This script is used to take raw enzyme data and convert it to rates (nmol/g/h).
This is based on German et al., 2011 (DOI: 10.1016/j.soilbio.2011.03.017); note the corrigendum with changes
to equations.  These have been taken into account.
