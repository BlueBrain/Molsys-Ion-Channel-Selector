# Author: Yann Roussel <yann.roussel@epfl.ch>
#
# License:
"""
This script upload the output files from the ic_selector_script.py to NEXUS.
Dependencies: pandas, getpass and kgforge
"""

import getpass

from kgforge.core import KnowledgeGraphForge
from kgforge.core import Resource
from kgforge.specializations.resources import Dataset


### TO DO: change the way of getting TOKEN
TOKEN = getpass.getpass()
# TOKEN = ""


nexus_endpoint = "https://bbp.epfl.ch/nexus/v1" # production environment

ORG = "bbp"
PROJECT = "ncmv3"

forge = KnowledgeGraphForge("https://raw.githubusercontent.com/BlueBrain/nexus-forge/master/examples/notebooks/use-cases/prod-forge-nexus.yml",
                           endpoint=nexus_endpoint,
                           bucket=f"{ORG}/{PROJECT}",
                           token= TOKEN,
                           debug=True
                           )


my_data_distribution = forge.attach("./output/met_type_ion_channel_gene_expression.csv")

# brainLocation as a Resource => so that you can do #my_dataset.brainLocation.brainRegion
brainRegion = Resource(label="Isocortex")
brainLocation = Resource(brainRegion=brainRegion)
 
my_dataset = Dataset(forge, type=["Entity","Dataset", "RNASequencing"],
                     name="Mouse_met_types_ion_channel_expression",
                     brainLocation = brainLocation,
                     description="Output from IC_selector module"
                    )
my_dataset.add_distribution("./output/met_type_ion_channel_gene_expression.csv")
forge.register(my_dataset)