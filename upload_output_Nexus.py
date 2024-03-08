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
# TOKEN = getpass.getpass()
TOKEN = "eyJhbGciOiJSUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICI5T0R3Z1JSTFVsTTJHbFphVDZjVklnenJsb0lzUWJmbTBDck1icXNjNHQ4In0.eyJleHAiOjE3MDk5MTg2NTUsImlhdCI6MTcwOTg4OTg1NSwiYXV0aF90aW1lIjoxNzA5ODg5ODU0LCJqdGkiOiIwNDRlNmE1Ni02N2ExLTQzNGQtYmI1Zi1mMTMwYTliMDE5N2EiLCJpc3MiOiJodHRwczovL2JicGF1dGguZXBmbC5jaC9hdXRoL3JlYWxtcy9CQlAiLCJhdWQiOlsiaHR0cHM6Ly9zbGFjay5jb20iLCJjb3Jlc2VydmljZXMtZ2l0bGFiIiwiYWNjb3VudCJdLCJzdWIiOiJmOjBmZGFkZWY3LWIyYjktNDkyYi1hZjQ2LWM2NTQ5MmQ0NTljMjp5cm91c3NlbCIsInR5cCI6IkJlYXJlciIsImF6cCI6ImJicC1uaXNlLW5leHVzLWZ1c2lvbiIsIm5vbmNlIjoiNTY5MTgzNjYyYjNkNDA4ZmJlMTY0MWZlMWY3YTY0MDIiLCJzZXNzaW9uX3N0YXRlIjoiZTJmY2YyY2ItM2I2Ni00YzhmLWJlZmEtMmM5MjkyZDVhZTRlIiwicmVhbG1fYWNjZXNzIjp7InJvbGVzIjpbImJicC1wYW0tYXV0aGVudGljYXRpb24iLCJvZmZsaW5lX2FjY2VzcyIsInVtYV9hdXRob3JpemF0aW9uIiwiZGVmYXVsdC1yb2xlcy1iYnAiXX0sInJlc291cmNlX2FjY2VzcyI6eyJodHRwczovL3NsYWNrLmNvbSI6eyJyb2xlcyI6WyJyZXN0cmljdGVkLWFjY2VzcyJdfSwiY29yZXNlcnZpY2VzLWdpdGxhYiI6eyJyb2xlcyI6WyJyZXN0cmljdGVkLWFjY2VzcyJdfSwiYWNjb3VudCI6eyJyb2xlcyI6WyJtYW5hZ2UtYWNjb3VudCIsIm1hbmFnZS1hY2NvdW50LWxpbmtzIiwidmlldy1wcm9maWxlIl19fSwic2NvcGUiOiJvcGVuaWQgbmV4dXMgcHJvZmlsZSBsb2NhdGlvbiBlbWFpbCIsInNpZCI6ImUyZmNmMmNiLTNiNjYtNGM4Zi1iZWZhLTJjOTI5MmQ1YWU0ZSIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJuYW1lIjoiWWFubiBSb3Vzc2VsIiwibG9jYXRpb24iOiJCMSA1IDI1OC4wNDAiLCJwcmVmZXJyZWRfdXNlcm5hbWUiOiJ5cm91c3NlbCIsImdpdmVuX25hbWUiOiJZYW5uIiwiZmFtaWx5X25hbWUiOiJSb3Vzc2VsIiwiZW1haWwiOiJ5YW5uLnJvdXNzZWxAZXBmbC5jaCJ9.irrP9LpintK_M3GUg2ujF2iF9PL4pYe98e4EUijFmMWIcUyJUPIwqGxTxJxvSEEnnc3UdJxhvecaEGQER-bjZwZJR8Zepe-Y0gPt5mvirZwRYBAMC5xrVywyim8rIdK8pxxpvsdpiVbPnjVsXf6zbLE0ClV8H14qP7ejxroQuE_OYlxPeatVy3H_zXuVfhc6ZAdYzE5LTEbCMMDovwlqXNTglZRz1UQBnlQ7qWwq6eLFlUoo2yP0QoRIp778N_9GXpdX1hwUu0TpaZphNnoZl21_6c9cMRY9PvvAWaoLbjTfkb-R_2C9Dqv8adHFMi5nr150BqMoUTNL0Aq3ghM1gw"


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
