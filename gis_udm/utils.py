'''
DESCRIPTION / INSTRUCTIONS

'''
# Modules
import shapely
from shapely import wkt
from shapely.geometry import Point
from shapely.ops import split
import geopandas as gpd


# Functions
def boundary_from_sewer_network(network, minlen, buffer_network, aux_buffer, fill_area, simplify):
    '''
    network: GeoDataFrame
    minlen: Minimum length of conduits to be considered
    buffer_network: Buffer distance from conduit
    aux_buffer: Auxiliar buffer distance
    fill_area: Maximum area considered to fill holes in polygons
    simplify: Simplify geometries parameter

    '''
    # Function
    crs = network.crs

    points = gpd.GeoSeries([Point(i) for j in network.geometry for i in list(j.coords)[1:-1]]).unary_union
    network = network.unary_union
    network_exploded = gpd.GeoSeries(split(network, points)).explode()

    filtered = network_exploded[network_exploded.geometry.length >= minlen]
    filtered_buffer = filtered.apply(lambda x: x.buffer(buffer_network, cap_style=2, join_style=2)).unary_union
    filtered_aux_buffer = filtered_buffer.buffer(aux_buffer, cap_style=2, join_style=2).buffer(-aux_buffer, cap_style=2, join_style=2)
    raw_extent = gpd.GeoSeries(filtered_aux_buffer).explode()

    new_geom = []
    for pol in raw_extent:
        rings = [i.replace('POLYGON ((', '').replace('))', '') for i in pol.wkt.split('), (')]
        rings = ['POLYGON ((' + i + '))' for i in rings]
        rings = [wkt.loads(i) for i in rings]
        outer = rings[0]
        holes = rings[1::]
        keep_holes = [i for i in holes if i.area > fill_area]
        aux_ = []
        for i in keep_holes:
            aux_.append(i.wkt.replace('POLYGON ((', '(').replace('))', ')'))
        aux_.insert(0, outer.wkt.replace('POLYGON ((', '(').replace('))', ')'))
        new_geom.append(shapely.wkt.loads('POLYGON (' + ', '.join(aux_) + ')'))

    return gpd.GeoDataFrame(geometry=new_geom, crs=crs).simplify(simplify)
