import pandas as pd
import geopandas as gpd
from osgeo import ogr
import shapely


def segmentize(geom, maxlength):
    wkt = geom.wkt  # shapely Polygon to wkt
    geom = ogr.CreateGeometryFromWkt(wkt)  # create ogr geometry
    geom.Segmentize(maxlength)  # densify geometry
    wkt2 = geom.ExportToWkt()  # ogr geometry to wkt
    new = shapely.wkt.loads(wkt2)  # wkt to shapely Polygon
    return new
