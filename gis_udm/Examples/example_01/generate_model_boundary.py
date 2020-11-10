# import sys
# import importlib
# importlib.reload(sys.modules['gis_udm_tools.land_cover'])

#%%
# Modules
import os
import geopandas as gpd
from gis_udm import utils

#%%
# Input/Output database (Geopackage) 
# use os.path to create paths that are platform independent
input_gpkg = os.path.join(os.getcwd(), 'input_dbase.gpkg')
output_gpkg = os.path.join(os.getcwd(), 'output_dbase.gpkg')

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
#%%
# Create a boundary, based on the system's conduits, for the generation of land cover topology
conduits = gpd.read_file(input_gpkg, layer='conduits', driver='GPKG')
# conduits.crs
# conduits.plot()

# boundary
input_bnd = utils.boundary_from_sewer_network(
    network=conduits,
    minlen=5,
    buffer_network=60,
    aux_buffer=20,
    fill_area=1e06,
    simplify=5
)
# input_bnd.crs
# input_bnd.plot()
input_bnd.to_file(input_gpkg, layer='boundary', driver='GPKG')
