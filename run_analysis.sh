#!/bin/bash

#Rscript ./prepare_brms_data.R &&
#Rscript ./gep_analysis_to_run_in_background_exclude_April_May.R
Rscript ./gep_analysis_to_run_in_background_include_April_May.R &&
Rscript ./GEP_examine_final.R

