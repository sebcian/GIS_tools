# Modules
import os
import shutil
import tempfile
from shapely import wkt
from shapely.geometry import Point
from shapely.ops import split
import pandas as pd
import geopandas as gpd
import numpy as np

import wrappers.saga_session as saga
import wrappers.grass_session as grass


# Functions
def thiessen_points(linestring, split_len):
    segments = int(linestring.length // (split_len / 2))
    fractions = np.array(range(1, segments, 2)) / segments
    fractions = np.array([0.5]) if len(fractions) == 1 else fractions
    return [linestring.interpolate(i, normalized=True) for i in fractions]


# Classes
class Subcatchments(object):

    def __init__(self, crs):

        self.tmpdir = tempfile.mkdtemp(prefix='tmp_subcatchments_')
        if not os.path.exists(self.tmpdir):
            os.makedirs(self.tmpdir)

        print('GRASS GIS session:\n')
        gisdb = os.path.join(self.tmpdir, 'grassdata')
        self.crs = crs
        self.gs = grass.GrassSession(gisdb=gisdb, location='subcatchments', crs=self.crs)
        print('\n\n')
        print('SAGA GIS session:\n')
        saga.get_version()


    def voronoi_delineation_from_conduits(
        self, conduits, conduits_id_col, split_len, area_thr, mask, prefix
    ):

        # Create temporary directory
        dir_ = os.path.join(self.tmpdir, 'voronoi_delination')
        if not os.path.exists(dir_):
            os.makedirs(dir_)

        # Split linestrings with more than two vertices
        conduits = conduits.explode()
        conduits.reset_index(drop=True, inplace=True)
        
        pts_dict = {}
        for idx, geo in zip(conduits.index, conduits.geometry):
            if len(list(geo.coords)[1:-1]) >= 1:
                pts_dict[idx] = gpd.GeoSeries([Point(p) for p in list(geo.coords)[1:-1]]).unary_union
            
        for k, v in pts_dict.items():
            conduits.geometry[k] = split(conduits.geometry[k], v)
        
        conduits = conduits.explode()
        conduits.reset_index(drop=True, inplace=True)
        
        # Create points used for voronoi diagram
        points_lists = conduits.geometry.apply(lambda x: thiessen_points(x, split_len))
        voronoi_points = pd.DataFrame([(idx, point) for idx in points_lists.index for point in points_lists[idx]])
        voronoi_points.set_index(0, drop=True, inplace=True)
        voronoi_points.columns = ['geometry']
        voronoi_points.index.name = None
        voronoi_points = gpd.GeoDataFrame(voronoi_points)
        data = conduits.loc[voronoi_points.index, [i for i in conduits.columns if 'geometry' not in i]]
        data.reset_index(drop=True, inplace=True)
        voronoi_points.reset_index(drop=True, inplace=True)
        voronoi_points = data.merge(voronoi_points, left_index=True, right_index=True)
        voronoi_points.crs = self.crs
        
        # Execute SAGA-GIS thiessen polygons
        # export points to temporary folder
        p_path = os.path.join(dir_, 'points_v.shp')
        voronoi_points.to_file(p_path)
        # thiessen path
        t_path = os.path.join(dir_, 'thiessen.shp')
        # execute
        saga.run_command('shapes_points', 16, points=p_path, polygons=t_path, frame=1000)

        # Create GeoDataFrame with thiessen
        thiessen = gpd.read_file(t_path)
        thiessen.columns = voronoi_points.columns.values

        # Aggregate by id
        thiessen = thiessen.dissolve(by=conduits_id_col)
        thiessen.reset_index(inplace=True)
        geom = thiessen['geometry'].explode()
        idx = geom.index.get_level_values(0)
        thiessen = thiessen.loc[idx, [i for i in thiessen.columns if 'geometry' not in i]]
        thiessen['geometry'] = geom.values
        thiessen = gpd.GeoDataFrame(thiessen)
        thiessen.crs = self.crs

        # Mask thiessen
        masked = gpd.overlay(thiessen, mask, how='intersection')
        geom = masked['geometry'].explode()
        idx = geom.index.get_level_values(0)
        thiessen = masked.loc[idx, [i for i in masked.columns if 'geometry' not in i]]
        thiessen['geometry'] = geom.values
        thiessen = gpd.GeoDataFrame(thiessen)
        thiessen.crs = self.crs

        # Remove small areas
        temp_gpkg = os.path.join(dir_, 'thiessen.gpkg')
        thiessen.to_file(temp_gpkg, layer='thiessen', driver='GPKG')
        self.gs.run_command('v.in.ogr', input=temp_gpkg, layer='thiessen')
        self.gs.run_command('v.clean', input='thiessen', output='cleaned', tool='rmarea', threshold=area_thr, overwrite=True)
        flags = 'su' if os.path.exists(temp_gpkg) else 's'
        self.gs.run_command(
            'v.out.ogr', flags=flags, input='cleaned', output=temp_gpkg, output_layer='thiessen',
            format='GPKG', overwrite=True
        )

        thiessen = gpd.read_file(temp_gpkg, layer='thiessen', driver='GPKG')

        # Create subcatchments ids
        thiessen.reset_index(inplace=True, drop=True) # reset index
        n = len(str(len(thiessen.index)))
        idxs = pd.Series(thiessen.index)
        thiessen.insert(0, 'id', prefix + idxs.apply(lambda x: '{{:0{}}}'.format(n).format(x+1)))

        # Remove temporary
        if os.path.exists(dir_):
            shutil.rmtree(dir_)

        return thiessen


    def purge(self):
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)
