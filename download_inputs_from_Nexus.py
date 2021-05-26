from kgforge.core import KnowledgeGraphForge
from kgforge.core import Resource
from kgforge.specializations.resources import Dataset

import pandas as pd

import getpass

### TO DO: change the way of getting TOKEN
print("enter your Nexus password")
TOKEN = "eyJhbGciOiJSUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICI5T0R3Z1JSTFVsTTJHbFphVDZjVklnenJsb0lzUWJmbTBDck1icXNjNHQ4In0.eyJleHAiOjE2MjIwNDczNDIsImlhdCI6MTYyMjAxODU0NiwiYXV0aF90aW1lIjoxNjIyMDE4NTQyLCJqdGkiOiJkZWQ4M2UyZi02YjE1LTQ3ZWYtOGE4OS0wNTNkYmYyMTc5MGQiLCJpc3MiOiJodHRwczovL2JicGF1dGguZXBmbC5jaC9hdXRoL3JlYWxtcy9CQlAiLCJzdWIiOiJmOjBmZGFkZWY3LWIyYjktNDkyYi1hZjQ2LWM2NTQ5MmQ0NTljMjp5cm91c3NlbCIsInR5cCI6IkJlYXJlciIsImF6cCI6Im5leHVzLXdlYiIsIm5vbmNlIjoiNGU4MzIzYTkxZTgyNDYyZmI5YzdlNTliNDEwMjdmNWMiLCJzZXNzaW9uX3N0YXRlIjoiMTI4NDE2ZTktMzgwNi00MjIxLThlMjQtYWFhYmFlYzViNzkzIiwiYWNyIjoiMCIsImFsbG93ZWQtb3JpZ2lucyI6WyJodHRwczovL2Rldi5uZXh1cy5vY3AuYmJwLmVwZmwuY2giLCJodHRwczovL2JicC5lcGZsLmNoIiwiaHR0cDovL2Rldi5uZXh1cy5vY3AuYmJwLmVwZmwuY2giLCJodHRwczovL3N0YWdpbmcubmV4dXMub2NwLmJicC5lcGZsLmNoIiwiaHR0cHM6Ly9iYnAtbmV4dXMuZXBmbC5jaCIsImh0dHBzOi8vYmJwdGVhbS5lcGZsLmNoIiwiaHR0cDovL3N0YWdpbmcubmV4dXMub2NwLmJicC5lcGZsLmNoIl0sInNjb3BlIjoib3BlbmlkIHByb2ZpbGUgZ3JvdXBzIGVtYWlsIiwiZW1haWxfdmVyaWZpZWQiOnRydWUsIm5hbWUiOiJZYW5uIFJvdXNzZWwiLCJwcmVmZXJyZWRfdXNlcm5hbWUiOiJ5cm91c3NlbCIsImdpdmVuX25hbWUiOiJZYW5uIiwiZmFtaWx5X25hbWUiOiJSb3Vzc2VsIiwiZW1haWwiOiJ5YW5uLnJvdXNzZWxAZXBmbC5jaCJ9.iPrk2dczy-9ZS0CGHifXlfajmz6Tu-2hdZvk5ggP_W-J_0nnj2a8avrSzcb0-rDDPnm_yjdOez5LIhuW74rhAIAWVOlVZh8dvqc4Lnl1NC5W9zaeuVylh_Ruvdr4VsR6_o02fHq5R-yFuadCYXScmH1i41TBH50FJ6Rofxn39XLOmQIO9GQ7t5-GUdM4RlIPUenOQXe7VAIejjX6b3CtfycJwvia5ZGNFByYd9coPk4XuYPqwjj1eAqpzJvo-YUdIIvZuus7IEn1XDsN7V5QnwYPvvufmmI_etx5nLkHrlaBIgV6r3xmu5OqlubKV3nxXDHghJPyAsUh1g_KG1uG8w" #getpass.getpass()
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