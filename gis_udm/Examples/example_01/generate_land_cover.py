# import sys
# import importlib
# importlib.reload(sys.modules['gis_udm_tools.land_cover'])

#%%
# Modules
import os
import geopandas as gpd
from gis_udm import utils, land_cover

#%%
# Input/Output database (Geopackage) 
# use os.path to create paths that are platform independent
input_gpkg = os.path.join(os.getcwd(), 'input_dbase.gpkg')
output_gpkg = os.path.join(os.getcwd(), 'output_dbase.gpkg')

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

# Generate Land Cover Topology
#%%
# Create land cover class
lc = land_cover.LandCover(input_gpkg, output_gpkg, crs='EPSG:31370')

#%%
# import boundary
lc.preprocess_boundary('boundary', simplify=1)

#%%
# import land cover layers
lc.preprocess_land_cover_layer(
    layer='building_blocks',
    priority=1,
    boundary='boundary',
    seglen=2,
    snap=1,
    simplify=1,
    area_thr=10
)

lc.preprocess_land_cover_layer(
    layer='permeable',
    priority=2,
    boundary='boundary',
    seglen=2,
    snap=1,
    simplify=1,
    area_thr=10
)       

#%%
# generate land cover interactive
layers = ['building_blocks', 'permeable']
priorities = [1, 2]
lc.generate_land_cover_interactive(
    layers=layers,
    priorities=priorities,
    boundary='boundary',
    seglen=1,
    area_thr=10,
    snap=1,
    simplify=1,
    output='land_cover'
)

#%%
# generate land cover static
# (this is an alternative method to generate the land cover polygons, however, land cover interactive is preferred.)
layers = ['building_blocks', 'permeable']
priorities = [1, 2]
lc.generate_land_cover_static(
    layers=layers,
    priorities=priorities,
    boundary='boundary',
    resolution=0.5,
    area_thr=10,
    snap=1,
    simplify=1,
    seglen=1,
    output='land_cover_static'
)

#%%
# consolidate topology land cover interactive
lc.consolidate_topology(
    boundary='boundary',
    land_cover='land_cover',
    seglen=2,
    snap=1,
    simplify=1,
    area_thr=10,
    output_boundary='boundary',
    output_land_cover='land_cover'
)

# consolidate topology land cover static
lc.consolidate_topology(
    boundary='boundary',
    land_cover='land_cover_static',
    seglen=2,
    snap=1,
    simplify=1,
    area_thr=10,
    output_boundary='boundary',
    output_land_cover='land_cover_static'
)

#%%
# purge land cover topology class
lc.purge()

# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
