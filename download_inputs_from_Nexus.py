# Author: Yann Roussel <yann.roussel@epfl.ch>
#
# License:
"""
This script download from NEXUS all the required input files for the ic_selector_script.py.
Dependencies: pandas, getpass and kgforge
"""

import getpass
import pandas as pd

from kgforge.core import KnowledgeGraphForge
from kgforge.core import Resource
from kgforge.specializations.resources import Dataset


### TO DO: change the way of getting TOKEN
print("enter your Nexus password")
TOKEN = getpass.getpass()
print("password taken")


nexus_endpoint = "https://bbp.epfl.ch/nexus/v1" # production environment

ORG = "bbp"
PROJECT = "ncmv3"

forge = KnowledgeGraphForge("https://raw.githubusercontent.com/BlueBrain/nexus-forge/master/examples/notebooks/use-cases/prod-forge-nexus.yml",
                           endpoint=nexus_endpoint,
                           bucket=f"{ORG}/{PROJECT}",
                           token= TOKEN,
                           debug=True
                           )
name_list = ["BBP_mtype_list", "mouse-whole-cortex-and-hippocampus-smart-seq",
             "P(marker_BBPmetype)_L1", "P(marker_BBPmetype)_L23_L6"]

for name in name_list:
    print(name)
    filters = {"type":"Dataset", "name":name}
    results = forge.search(filters, limit=3)
    print(f"{len(results)} results found")

    forge.download(results, "distribution.contentUrl", path="./input/")
    print("________")