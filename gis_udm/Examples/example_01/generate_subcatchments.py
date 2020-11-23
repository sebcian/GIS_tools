# import sys
# import importlib
# importlib.reload(sys.modules['gis_udm_tools.land_cover'])

#%%
# Modules
import os
import geopandas as gpd
from gis_udm import utils, subcatchments

#%%
# Input/Output database (Geopackage) 
# use os.path to create paths that are platform independent
input_gpkg = os.path.join(os.getcwd(), 'input_dbase.gpkg')
output_gpkg = os.path.join(os.getcwd(), 'output_dbase.gpkg')

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Subcatchments generation
#%%
# Generate mask
conduits = gpd.read_file(input_gpkg, layer='conduits', driver='GPKG')[['drainage_id', 'system_type','geometry']]
conduits_combined = conduits.loc[conduits['system_type'] == 'combined', :]

boundary_combined = utils.boundary_from_sewer_network(
    network=conduits_combined,
    minlen=5,
    buffer_network=60,
    aux_buffer=20,
    fill_area=1e06,
    simplify=5
)
# boundary_combined.crs
# boundary_combined.plot()
boundary_combined = gpd.GeoDataFrame(geometry=boundary_combined)
# boundary_combined.crs

subcatchments_mask = gpd.read_file(input_gpkg, layer='subcatchments_mask', driver='GPKG')
# subcatchments_mask.crs
# subcatchments_mask.plot()

mask = gpd.overlay(subcatchments_mask, boundary_combined, how='intersection')
# mask.crs
# mask.plot()

#%%
# Subcatchments
sub = subcatchments.Subcatchments('epsg:31370')
thiessen = sub.voronoi_delineation_from_conduits(
    conduits=conduits_combined,
    conduits_id_col='drainage_id',
    split_len=2.5,
    area_thr=10,
    mask=mask,
    prefix='combined_subcatchment_'
)
# thiessen.plot()
# thiessen.to_file(output_gpkg, layer='combined_thiessen', driver='GPKG')

#%%
# Purge sub
sub.purge()

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
