#!/bin/bash

# Download necessary data from Nexus
python3 download_inputs_from_Nexus.py

mkdir output

# Compute the ion channel expression profiles for each me-type+t-type combination
python3 ic_selector_script.py

# Upload the final output file (met-types ion channel expression profiles) to Nexus
#python3 upload_output_Nexus.py
