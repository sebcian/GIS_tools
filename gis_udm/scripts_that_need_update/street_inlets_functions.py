# General Modules
import os
import sys
import subprocess
import yaml
import numpy as np
import pandas as pd

# General modules
import general_functions as gn
import sqlite_functions as sq

# Qgis modules
from initialize_qgis import *


# tmp directory
tmp_dir = os.path.join(os.getcwd(), 'tmp_dir')
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)


def connections_hub(si={}, cs={}, en={}, min_separation=float()):

    si_layer = si['layer']
    si_id_col = si['id_col']

    cs_layer = cs['layer']
    cs_id_col = cs['id_col']
    cs_systype_col = cs['systype_col']

    en_layer = en['layer']
    en_id_col = en['id_col']
    en_node_type_col = en['node_type_col']
    en_node_type_not = en['node_type_not_connect']
    en_node_systype_col = en['systype_col']

    # Buffer to existing network nodes
    tool = 'native:buffer'
    parameters = {'INPUT': en_layer,
                  'DISTANCE': min_separation,
                  'SEGMENTS': 5,
                  'END_CAP_STYLE': 0,
                  'JOIN_STYLE': 0,
                  'MITER_LIMIT': 2,
                  'DISSOLVE': True,
                  'OUTPUT': 'memory:'}

    e_nodes_buffer = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(e_nodes_buffer, 'test.sqlite', 'buffer', 'pk', 'geom')

    # Buffer to collection systems
    tool = 'native:buffer'
    parameters = {'INPUT': cs_layer,
                  'DISTANCE': 0.001,
                  'SEGMENTS': 5,
                  'END_CAP_STYLE': 1,
                  'JOIN_STYLE': 0,
                  'MITER_LIMIT': 2,
                  'DISSOLVE': True,
                  'OUTPUT': 'memory:'}

    cs_buffer = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(cs_buffer, 'test.sqlite', 'cs_buffer', 'pk', 'geom')

    # Points along collection_system vector layer
    tool = 'qgis:pointsalonglines'
    parameters = {'INPUT': cs_layer,
                  'DISTANCE': min_separation,
                  'START_OFFSET': 0,
                  'END_OFFSET': 0,
                  'OUTPUT': 'memory:'}
    p_along = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(p_along, 'test.sqlite', 'palong', 'pk', 'geom')

    # Remove duplicate geometries
    # gn.algorithms_list('delete')
    # gn.algorithms_help('qgis:deleteduplicategeometries')
    tool = 'qgis:deleteduplicategeometries'
    parameters = {'INPUT': p_along,
                  'OUTPUT': 'memory:'}
    p_along_cleaned = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(p_along, 'test.sqlite', 'palong_clean', 'pk', 'geom')

    # Select "p_along_cleaned" that disjoint e_nodes_buffer
    tool = 'native:selectbylocation'
    parameters = {'INPUT': p_along_cleaned,
                  'PREDICATE': 2,
                  'INTERSECT': e_nodes_buffer,
                  'METHOD': 0}

    select_ = processing.run(tool, parameters)['OUTPUT']
    # Save selected features
    tool = 'qgis:saveselectedfeatures'
    parameters = {'INPUT': select_,
                  'OUTPUT': 'memory:'}
    p_along_selection = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(p_along_selection, 'test.sqlite', 'palong_sel', 'pk', 'geom')

    # Manage attribute table of p_along_selection
    p_along_selection_pr = p_along_selection.dataProvider()
    del_fields = [p_along_selection.fields().indexFromName(i) for i in p_along_selection.fields().names() if i not in [cs_id_col, cs_systype_col]]
    p_along_selection_pr.deleteAttributes(del_fields)
    p_along_selection.updateFields()

    add_ = [QgsField(en_id_col, QVariant.String, typeName='string')]
    p_along_selection_pr.addAttributes(add_)
    p_along_selection.updateFields()
    n_digits = len(str(p_along_selection.featureCount()))
    idx_ = p_along_selection.fields().indexFromName(en_id_col)
    fid_dict = {}
    for i in p_along_selection.getFeatures():
        id_val = 'p_along_{{:0{}}}'.format(n_digits).format(i.id())
        attrs = {idx_: id_val}
        fid_dict[i.id()] = attrs
    p_along_selection_pr.changeAttributeValues(fid_dict)
    p_along_selection.updateFields()
    # sq.layer_to_spatialite(p_along_selection, 'test.sqlite', 'fields', 'pk', 'geom')

    mapping = []
    for i in p_along_selection.fields():
        if i.name() == cs_id_col:
            name_ = 'origin'
        else:
            name_ = i.name()
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': name_,
                        'precision': i.precision(),
                        'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': p_along_selection,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    p_along_refactored = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(p_along_refactored, 'test.sqlite', 'refact', 'pk', 'geom')

    # Select nodes from network nodes that intersect cs_buffer
    tool = 'native:selectbylocation'
    parameters = {'INPUT': en_layer,
                  'PREDICATE': 0,
                  'INTERSECT': cs_buffer,
                  'METHOD': 0}

    select_ = processing.run(tool, parameters)['OUTPUT']
    # Save selected features
    tool = 'qgis:saveselectedfeatures'
    parameters = {'INPUT': select_,
                  'OUTPUT': 'memory:'}
    en_selection = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(en_selection, 'test.sqlite', 'sel1', 'pk', 'geom')

    # Select nodes from "en_selection" that are not of specified node type (node_type_not_connect)
    # gn.algorithms_list('select')
    # gn.algorithms_help('qgis:selectbyexpression')
    tool = 'qgis:selectbyexpression'
    expr = ['{} <> \'{}\''.format(en_node_type_col, i) for i in en_node_type_not]
    expr = 'AND '.join(expr)
    parameters = {'INPUT': en_selection,
                  'EXPRESSION': expr,
                  'METHOD': 0}

    select_ = processing.run(tool, parameters)['OUTPUT']
    # Save selected features
    tool = 'qgis:saveselectedfeatures'
    parameters = {'INPUT': select_,
                  'OUTPUT': 'memory:'}
    en_selection2 = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(en_selection2, 'test.sqlite', 'sel', 'pk', 'geom')

    # Manage attribute table of en_selection2
    en_selection2_pr = en_selection2.dataProvider()
    del_fields = [en_selection2.fields().indexFromName(i) for i in en_selection2.fields().names() if i not in [en_id_col, en_node_type_col, en_node_systype_col]]
    en_selection2_pr.deleteAttributes(del_fields)
    en_selection2.updateFields()
    # sq.layer_to_spatialite(en_selection2, dbase='test.sqlite', name='en', pk='pk', geocol='geom')

    if not en_node_systype_col == cs_systype_col:
        mapping = []
        for i in en_selection2.fields():
            if i.name() == en_node_systype_col:
                name_ = cs_systype_col
            else:
                name_ = i.name()

            mapping.append({'expression': '\"{}\"'.format(i.name()),
                            'length': i.length(),
                            'name':name_,
                            'precision': i.precision(),
                            'type': i.type()})

        tool = 'qgis:refactorfields'
        parameters = {'INPUT': en_selection2,
                      'FIELDS_MAPPING': mapping,
                      'OUTPUT': 'memory:'}
        en_selection2 = processing.run(tool, parameters)['OUTPUT']
        # sq.layer_to_spatialite(si_refactored, 'test.sqlite', 'refact', 'pk', 'geom')

    # Merge selected nodes
    tool = 'native:mergevectorlayers'
    parameters = {'LAYERS': [en_selection2, p_along_refactored],
                  'CRS': en_layer.crs().authid(),
                  'OUTPUT': 'memory:'}

    n_along = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(layer=n_along, dbase='test.sqlite', name='merged', pk='pk', geocol='geom')

    # Nearest hub between street inlets and selected nodes
    # gn.algorithms_list('hub')
    # gn.algorithms_help('qgis:distancetonearesthublinetohub')
    tool = 'qgis:distancetonearesthublinetohub'
    parameters = {'INPUT': si_layer,
                  'HUBS': n_along,
                  'FIELD': en_id_col,
                  'UNIT': 0,
                  'OUTPUT': 'memory:'}

    hub = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(hub, 'test.sqlite', 'hub', 'pk', 'geom')

    # Managing attribute table
    hub_df = gn.layer_df(hub)
    n_along_df = gn.layer_df(n_along)
    hub_df['street_inlet_id'] = hub_df[si_id_col]
    aux = pd.DataFrame(hub_df['HubName'][hub_df['HubName'].str.startswith('p_along_')].unique(), columns=['HubName'])
    aux.index = range(1, len(aux.index) + 1)
    n_digits = len(str(len(aux.index)))
    aux['connection_node_id'] = ['connection_node_{{:0{}}}'.format(n_digits).format(i) for i in aux.index.values]
    hub_df = hub_df.merge(aux, how='left', on='HubName')
    idx = hub_df.index[pd.isnull(hub_df['connection_node_id'])]
    hub_df.loc[idx, 'connection_node_id'] = hub_df.loc[idx, 'HubName']
    hub_df['connection_node_origin'] = [n_along_df.loc[n_along_df[en_id_col] == i, 'origin'].values[0] for i in hub_df['HubName']]
    hub_df['system_type'] = [n_along_df.loc[n_along_df[en_id_col] == i, cs_systype_col].values[0] for i in hub_df['HubName']]

    # Add columns to attribute table
    hub_pr = hub.dataProvider()
    hub_pr.addAttributes([QgsField('street_inlet_id', QVariant.String, typeName='string'),
                          QgsField('connection_node_id', QVariant.String, typeName='string'),
                          QgsField('connection_node_origin', QVariant.String, typeName='string'),
                          QgsField('system_type', QVariant.String, typeName='string')])
    hub.updateFields()
    idx_si_id = hub.fields().indexFromName('street_inlet_id')
    idx_cn_id = hub.fields().indexFromName('connection_node_id')
    idx_origin = hub.fields().indexFromName('connection_node_origin')
    idx_systype = hub.fields().indexFromName('system_type')

    attrs_dict = {}
    for i, j, k, l, m in zip(hub.getFeatures(),
                             hub_df['street_inlet_id'], hub_df['connection_node_id'],
                             hub_df['connection_node_origin'], hub_df['system_type']):

        attrs_dict[i.id()] = {idx_si_id: j, idx_cn_id: k, idx_origin: l, idx_systype: m}
    hub_pr.changeAttributeValues(attrs_dict)
    hub.updateFields()
    # sq.layer_to_spatialite(hub, 'test.sqlite', 'hub', 'pk', 'geom')

    hub_dict = {'layer': hub,
                'street_inlet_id_col': 'street_inlet_id',
                'connection_node_id_col': 'connection_node_id',
                'connection_node_origin_col': 'connection_node_origin',
                'systype_col': 'system_type'}

    return hub_dict


def si_attribute_table(si={}, hub={}, dem=QgsRasterLayer()):

    si_layer = si['layer']
    si_id_col = si['id_col']

    hub_layer = hub['layer']
    hub_si_id_col = hub['street_inlet_id_col']
    hub_systype_col = hub['systype_col']

    # Street inlets
    # Create copy of layer
    si_layer.selectAll()
    tool = 'native:saveselectedfeatures'
    parameters = {'INPUT': si_layer,
                  'OUTPUT': 'memory:'}

    si_copy = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(hub, 'test.sqlite', 'hub', 'pk', 'geom')

    # Sytem type column
    si_copy_pr = si_copy.dataProvider()
    si_copy_pr.addAttributes([QgsField('system_type', QVariant.String, typeName='string')])
    si_copy.updateFields()
    idx = si_copy.fields().indexFromName('system_type')

    si_df = gn.layer_df(si_layer)
    hub_df = gn.layer_df(hub_layer)

    systype = [hub_df.loc[hub_df[hub_si_id_col] == i, hub_systype_col].values[0] for i in si_df[si_id_col]]
    systype_dict = {}
    for i, j in zip(si_layer.getFeatures(), systype):
        systype_dict[i.id()] = {idx: j}
    si_copy_pr.changeAttributeValues(systype_dict)
    si_copy.updateFields()
    # sq.layer_to_spatialite(si_refactored, 'test.sqlite', 'refact', 'pk', 'geom')

    # Geometry attributes
    tool = 'qgis:exportaddgeometrycolumns'
    parameters = {'INPUT': si_copy,
                  'CALC_METHOD': 0,
                  'OUTPUT': 'memory:'}

    si_geocols = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(hub, 'test.sqlite', 'hub', 'pk', 'geom')

    # Ground elevation
    tmp_shp = os.path.join(tmp_dir, 'tmp_shp.shp')
    tool = 'saga:addrastervaluestopoints'
    parameters = {'SHAPES': si_geocols,
                  'GRIDS': dem,
                  'RESAMPLING': 0,
                  'RESULT': tmp_shp}
    si_ground_elev = processing.run(tool, parameters)['RESULT']
    si_ground_elev = gn.read_shp(si_ground_elev)

    # Manage attribute table
    names = si_geocols.fields().names()
    names.append('ground_elevation')
    mapping = []
    for i, j in zip(si_ground_elev.fields(), names):
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': j,
                        'precision': i.precision(),
                        'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': si_ground_elev,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    si_refactored = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(si_refactored, 'test.sqlite', 'si', 'pk', 'geom')

    si_dict = {'layer': si_refactored,
               'id_col': si_id_col,
               'systype_col': 'system_type',
               'xcoord_col': 'xcoord',
               'ycoord_col': 'ycoord',
               'ground_elevation_col': 'ground_elevation'
               }

    return si_dict


def connection_nodes(hub={}, dem=QgsRasterLayer()):

    hub_layer = hub['layer']
    hub_conn_node_id_col = hub['connection_node_id_col']
    hub_systype_col = hub['systype_col']

    # Create connection nodes layer
    crs = hub_layer.crs().authid()
    conn_n = QgsVectorLayer('Point?crs={}&field=id:string&field=system_type:string'.format(crs), 'conn_n', 'memory')
    conn_n_pr = conn_n.dataProvider()

    # Get geometry and id and systype
    tuples_ = []
    for i in hub_layer.getFeatures():
        if i[hub_conn_node_id_col].startswith('connection_node_'):
            extract = i.geometry().asWkt().split(',')[-1].strip()
            geom_ = 'POINT(' + extract
            id_ = i[hub_conn_node_id_col]
            systype_ = i[hub_systype_col]
            tuples_.append((geom_, id_, systype_))

    # Remove duplicates (when street inlets connecting to the same connection node)
    tuples_ = list(set(tuples_))
    tuples_ = sorted(tuples_, key=lambda  tup: tup[1])

    # Create features
    for i, j, k in tuples_:
        feat = QgsFeature()
        geom_ = QgsGeometry.fromWkt(i)
        feat.setGeometry(geom_)
        feat.setAttributes([j, k])
        conn_n_pr.addFeatures([feat])

    conn_n_pr.updateExtents()
    # sq.layer_to_spatialite(conn_n, 'test.sqlite', 'conn_n', 'pk', 'geom')

    # Coordinates
    tool = 'qgis:exportaddgeometrycolumns'
    parameters = {'INPUT': conn_n,
                  'CALC_METHOD': 0,
                  'OUTPUT': 'memory:'}

    cn_geocols = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(hub, 'test.sqlite', 'hub', 'pk', 'geom')

    # Ground elevation
    tmp_shp = os.path.join(tmp_dir, 'tmp_shp.shp')
    tool = 'saga:addrastervaluestopoints'
    parameters = {'SHAPES': cn_geocols,
                  'GRIDS': dem,
                  'RESAMPLING': 0,
                  'RESULT': tmp_shp}
    cn_ground_elev = processing.run(tool, parameters)['RESULT']
    cn_ground_elev = gn.read_shp(cn_ground_elev)

    # Manage attribute table
    names = cn_geocols.fields().names()
    names.append('ground_elevation')
    mapping = []
    for i, j in zip(cn_ground_elev.fields(), names):
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': j,
                        'precision': i.precision(),
                        'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': cn_ground_elev,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    cn_refactored = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(cn_refactored, 'test.sqlite', 'cn', 'pk', 'geom')

    conn_n_dict = {'layer': cn_refactored,
                   'id_col': 'id',
                   'systype_col': 'system_type',
                   'xcoord_col': 'xcoord',
                   'ycoord_col': 'ycoord',
                   'ground_elevation_coord': 'ground_elevation'}

    return conn_n_dict


def connection_pipes(hub={}, si={}, cs={}, en={},
                     cp_us_offset=float(), cp_ds_offset=float(), cp_min_gradient=float()):

    hub_layer = hub['layer']
    hub_si_id_col = hub['street_inlet_id_col']
    hub_conn_node_id_col = hub['connection_node_id_col']
    hub_conn_node_origin_col = hub['connection_node_origin_col']
    hub_systype_col = hub['systype_col']

    si_layer = si['layer']
    si_id_col = si['id_col']
    si_ground_elev_col = si['ground_elevation_col']

    cs_layer = cs['layer']
    cs_id_col = cs['id_col']
    cs_us_node_id_col = cs['us_node_id_col']
    cs_ds_node_id_col = cs['ds_node_id_col']
    cs_us_invert_col = cs['us_invert_col']
    cs_ds_invert_col = cs['ds_invert_col']
    cs_conduit_height_col = cs['conduit_height_col']

    en_layer = en['layer']
    en_id_col = en['id_col']

    # Create connection pipes layer
    crs = hub_layer.crs().authid()
    conn_p = QgsVectorLayer('LineString?crs={}&field=id:string&field=system_type:string&field=street_inlet_id:string&field=connection_node_id:string'.format(crs), 'conn_p', 'memory')
    conn_p_pr = conn_p.dataProvider()

    n_digits = len(str(hub_layer.featureCount()))
    for i in hub_layer.getFeatures():
        geom_ = i.geometry()
        id_val_ = 'connection_pipe_{{:0{}}}'.format(n_digits).format(i.id())
        systype_ = i[hub_systype_col]
        si_id_ = i[hub_si_id_col]
        cn_id_ = i[hub_conn_node_id_col]
        feat = QgsFeature()
        feat.setGeometry(geom_)
        feat.setAttributes([id_val_, systype_, si_id_, cn_id_])
        conn_p_pr.addFeatures([feat])

    conn_p_pr.updateExtents()
    # sq.layer_to_spatialite(conn_p, 'test.sqlite', 'conn_p', 'pk', 'geom')

    # Additional attributes
    hub_df = gn.layer_df(hub_layer)

    cp_df = gn.layer_df(conn_p)
    si_df = gn.layer_df(si_layer)

    cs_df = gn.layer_df(cs_layer)
    # en_df_ = gn.layer_df(en_layer)

    df = pd.DataFrame(index=cp_df.index)

    # us invert levels
    si_elev = np.array([si_df.loc[si_df[si_id_col] == i, si_ground_elev_col].values[0] for i in cp_df['street_inlet_id']])
    cp_us_invert = si_elev - cp_us_offset

    df['us_invert'] = cp_us_invert

    # ds invert levels
    tuples_ = []
    for i in en_layer.getFeatures():
        id_ = i[en_id_col]
        geom_ = i.geometry().asWkt().replace('Point (', '').replace(')', '')
        x_, y_ = geom_.split(' ')
        tuples_.append((id_, x_, y_))
    en_coords = pd.DataFrame(tuples_, columns=['id', 'x', 'y'])

    tuples_ = []
    for i in hub_layer.getFeatures():
        id_ = i[hub_conn_node_origin_col]
        len_ = i.geometry().length()
        ds_point_ = i.geometry().asWkt().split(', ')[-1].replace(')', '')
        x_, y_ = ds_point_.split(' ')
        tuples_.append((id_, len_, x_, y_))
    cp_data = pd.DataFrame(tuples_, columns=['id', 'length', 'x', 'y'])

    ####################################################################################################################
    idx = hub_df.index[~pd.isnull(hub_df[hub_conn_node_origin_col])]

    with_origin = hub_df.loc[idx, hub_conn_node_origin_col]
    us_node = [cs_df.loc[cs_df[cs_id_col] == i, cs_us_node_id_col].values[0] for i in with_origin]
    us_x = np.array([en_coords.loc[en_coords['id'] == i, 'x'].values[0] for i in us_node], dtype=float)
    us_y = np.array([en_coords.loc[en_coords['id'] == i, 'y'].values[0] for i in us_node], dtype=float)
    us_invert = np.array([cs_df.loc[cs_df[cs_id_col] == i, cs_us_invert_col].values[0] for i in with_origin], dtype=float)

    ds_node = [cs_df.loc[cs_df[cs_id_col] == i, cs_ds_node_id_col].values[0] for i in with_origin]
    ds_x = np.array([en_coords.loc[en_coords['id'] == i, 'x'].values[0] for i in ds_node], dtype=float)
    ds_y = np.array([en_coords.loc[en_coords['id'] == i, 'y'].values[0] for i in ds_node], dtype=float)
    ds_invert = np.array([cs_df.loc[cs_df[cs_id_col] == i, cs_ds_invert_col].values[0] for i in with_origin], dtype=float)

    gradient = (us_invert - ds_invert) / ((ds_y - us_y) ** 2 + (ds_x - us_x) ** 2) ** 0.5
    conduit_height = np.array([cs_df.loc[cs_df[cs_id_col] == i, cs_conduit_height_col].values[0] for i in with_origin], dtype=float) / 1000.0 # conduit heigth given in mm in layer

    length = cp_data.loc[idx, 'length'].values.astype(float)
    cn_x = cp_data.loc[idx, 'x'].values.astype(float)
    cn_y = cp_data.loc[idx, 'y'].values.astype(float)

    dist_cn_csds = ((cn_y - ds_y) ** 2 + (cn_x - ds_x) ** 2) ** 0.5

    # Calculation of invert levels at downstream
    cp_ds_min_invert = ds_invert + gradient * dist_cn_csds  # ds in connection pipe
    cp_ds_max_invert = cp_ds_min_invert + conduit_height - cp_ds_offset
    cp_ds_invert = cp_us_invert[idx] - cp_min_gradient * length

    aux_idx_max = (cp_ds_max_invert - cp_ds_invert) < 0
    cp_ds_invert[aux_idx_max] = cp_ds_max_invert[aux_idx_max]

    aux_idx_min = (cp_ds_invert - cp_ds_min_invert) < 0
    cp_ds_invert[aux_idx_min] = cp_ds_min_invert[aux_idx_min]

    df.loc[idx, 'ds_invert'] = cp_ds_invert
    ####################################################################################################################

    ####################################################################################################################
    idx = hub_df.index[pd.isnull(hub_df[hub_conn_node_origin_col])]
    length = cp_data.loc[idx, 'length'].values.astype(float)
    cp_ds_invert = cp_us_invert[idx] - cp_min_gradient * length
    df.loc[idx, 'ds_invert'] = cp_ds_invert
    ####################################################################################################################

    df['length'] = cp_data['length'].values

    # Add to layer attribute table
    conn_p_pr = conn_p.dataProvider()
    conn_p_pr.addAttributes([QgsField('us_invert', QVariant.Double, typeName='double'),
                             QgsField('ds_invert', QVariant.Double, typeName='double'),
                             QgsField('length', QVariant.Double, typeName='double')])
    conn_p.updateFields()
    idx_us_invert = conn_p.fields().indexFromName('us_invert')
    idx_ds_invert = conn_p.fields().indexFromName('ds_invert')
    idx_length = conn_p.fields().indexFromName('length')

    attrs_dict = {}
    for i, j, k, l in zip(conn_p.getFeatures(),
                          df['us_invert'],
                          df['ds_invert'],
                          df['length']):

        attrs_dict[i.id()] = {idx_us_invert: j, idx_ds_invert: k, idx_length: l}
    conn_p_pr.changeAttributeValues(attrs_dict)
    conn_p.updateFields()
    # sq.layer_to_spatialite(conn_p, 'test.sqlite', 'conn_p', 'pk', 'geom')

    conn_p_dict = {'layer': conn_p,
                   'id_col': 'id',
                   'systype_col': 'system_type',
                   'street_inlet_id_col': 'street_inlet_id',
                   'connection_node_id_col': 'connection_node_id',
                   'us_invert_col': 'us_invert',
                   'ds_invert_col': 'ds_invert',
                   'length_col': 'length'}

    return conn_p_dict



def split_cs(cs, cn, en, fc):

    cs_layer = cs['layer']
    cs_id_col = cs['id_col']
    cs_us_node_id_col = cs['us_node_id_col']
    cs_ds_node_id_col = cs['ds_node_id_col']
    cs_us_invert_col = cs['us_invert_col']
    cs_ds_invert_col = cs['ds_invert_col']
    cs_link_suffix_col = cs['link_suffix_col']
    cs_gradient_col = cs['gradient_col']

    cn_layer = cn['layer']
    cn_id_col = cn['id_col']

    en_layer = en['layer']
    en_id_col = en['id_col']

    fc_layer = fc['layer']
    fc_id_col = fc['id_col']
    fc_us_node_id_col = fc['us_node_id_col']
    fc_link_suffix_col = fc['link_suffix_col']

    # Split collection system at connection nodes
    tmp_shp = os.path.join(tmp_dir, 'tmp_shp.shp')
    tool = 'saga:splitlinesatpoints'
    parameters = {'LINES': cs_layer,
                  'SPLIT': cn_layer,
                  'INTERSECT': tmp_shp,
                  'OUTPUT': 1,
                  'EPSILON': 0.001}
    cs_split = processing.run(tool, parameters)['INTERSECT']
    cs_split = gn.read_shp(cs_split)
    # sq.layer_to_spatialite(cn_refactored, 'test.sqlite', 'cn', 'pk', 'geom')

    # Manage attribute table
    names = cs_layer.fields().names()
    mapping = []
    for i, j in zip(cs_split.fields(), names):
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': j,
                        'precision': i.precision(),
                        'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': cs_split,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    cs_refactored = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(si_refactored, 'test.sqlite', 'si', 'pk', 'geom')

    # Multipart to single part
    cs_single = processing.run('native:multiparttosingleparts', parameters={'INPUT': cs_refactored,
                                                                            'OUTPUT': 'memory:'})['OUTPUT']

    # Add information to attribute table
    tuples_ = []
    for i in cs_single.getFeatures():
        fid_ = i.id()
        id_ = i[cs_id_col]
        us_node_ = i[cs_us_node_id_col]
        ds_node_ = i[cs_ds_node_id_col]
        us_invert_ = i[cs_us_invert_col]
        ds_invert_ = i[cs_ds_invert_col]
        gradient_ = i[cs_gradient_col]
        geom_ = i.geometry().asWkt()
        x1_ = float(geom_.split(', ')[0].replace('LineString (', '').split(' ')[0])
        y1_ = float(geom_.split(', ')[0].replace('LineString (', '').split(' ')[-1])

        x2_ = float(geom_.split(', ')[-1].replace(')', '').split(' ')[0])
        y2_ = float(geom_.split(', ')[-1].replace(')', '').split(' ')[-1])

        tuples_.append((fid_, id_, us_node_, ds_node_, us_invert_, ds_invert_, gradient_, x1_, y1_, x2_, y2_))

    df_cp = pd.DataFrame(tuples_, columns=['fid', 'id', 'us_node', 'ds_node', 'us_invert', 'ds_invert', 'gradient', 'x1', 'y1', 'x2', 'y2'])
    df_cp.index = df_cp['fid']

    tuples_ = []
    for i in en_layer.getFeatures():
        id_ = i[en_id_col]
        geom_ = i.geometry().asWkt()
        x_ = float(geom_.replace('Point (', '').split(' ')[0])
        y_ = float(geom_.replace(')', '').split(' ')[-1])
        tuples_.append((id_, x_, y_))

    for i in cn_layer.getFeatures():
        id_ = i[cn_id_col]
        geom_ = i.geometry().asWkt()
        x_ = float(geom_.replace('Point (', '').split(' ')[0])
        y_ = float(geom_.replace(')', '').split(' ')[-1])
        tuples_.append((id_, x_, y_))

    df_coord_n = pd.DataFrame(tuples_, columns=['id', 'x', 'y'])
    df_coord_n.index = df_coord_n['id']

    id1 = np.array([df_coord_n.loc[((abs(df_coord_n['x'] - i) < 0.001) & (abs(df_coord_n['y'] - j) < 0.001)), 'id'].values[0] for i, j in zip(df_cp['x1'], df_cp['y1'])])
    id2 = np.array([df_coord_n.loc[((abs(df_coord_n['x'] - i) < 0.001) & (abs(df_coord_n['y'] - j) < 0.001)), 'id'].values[0] for i, j in zip(df_cp['x2'], df_cp['y2'])])

    # Distances between original us_node and new nodes
    us_node_x = np.array([df_coord_n.loc[df_coord_n['id'] == i, 'x'].values[0] for i in df_cp['us_node']],dtype=float)
    us_node_y = np.array([df_coord_n.loc[df_coord_n['id'] == i, 'y'].values[0] for i in df_cp['us_node']],dtype=float)
    d1 = np.array(((df_cp['y1'] - us_node_y) ** 2 + (df_cp['x1'] - us_node_x) ** 2) ** 0.5, dtype=float)
    d2 = np.array(((df_cp['y2'] - us_node_y) ** 2 + (df_cp['x2'] - us_node_x) ** 2) ** 0.5, dtype=float)
    
    # Assigning new nodes and invert levels
    new_df = pd.DataFrame(index=df_cp.index)
    idx = d1 < d2
    new_df.loc[idx, 'us_node_id'] = id1[idx]
    new_df.loc[idx, 'ds_node_id'] = id2[idx]
    new_df.loc[idx, 'us_invert'] = df_cp.loc[idx, 'us_invert'] - df_cp.loc[idx, 'gradient'] * d1[idx]
    new_df.loc[idx, 'ds_invert'] = df_cp.loc[idx, 'us_invert'] - df_cp.loc[idx, 'gradient'] * d2[idx]

    idx = d1 > d2
    new_df.loc[idx, 'us_node_id'] = id2[idx]
    new_df.loc[idx, 'ds_node_id'] = id1[idx]
    new_df.loc[idx, 'us_invert'] = df_cp.loc[idx, 'us_invert'] - df_cp.loc[idx, 'gradient'] * d2[idx]
    new_df.loc[idx, 'ds_invert'] = df_cp.loc[idx, 'us_invert'] - df_cp.loc[idx, 'gradient'] * d1[idx]

    # Link suffix and id
    fc_df = gn.layer_df(fc_layer)
    us_nodes_unique = new_df['us_node_id'].unique()
    for n in us_nodes_unique:
        idx = new_df.index[new_df['us_node_id'] == n].values
        if n in fc_df[fc_us_node_id_col].values:
            not_available = fc_df.loc[fc_df[fc_us_node_id_col] == n, fc_link_suffix_col].astype(int).values
            c = 0
            i = 1
            new_values = []
            while c < len(idx):
                if i not in not_available:
                    new_values.append(i)
                    c += 1
                i += 1
        else:
            new_values = [int(i + 1) for i in range(len(idx))]

        new_df.loc[idx, cs_link_suffix_col] = new_values

    new_df[cs_link_suffix_col] = new_df[cs_link_suffix_col].astype(int)

    new_df[cs_id_col] = new_df[cs_us_node_id_col].astype(str) + '.' + new_df[cs_link_suffix_col].astype(str)
    # new_df = new_df.astype(str)

    # Update attribute table with new values
    cs_single_pr = cs_single.dataProvider()
    tuples_ = [(cs_single.fields().indexFromName(i), i) for i in new_df.columns.values]
    tuples_ = sorted(tuples_, key=lambda tup: tup[0])
    attrs_dict = {}
    for i in cs_single.getFeatures():
        attrs = {}
        for j in tuples_:
            attrs[j[0]] = new_df.loc[i.id(), j[1]]
        attrs_dict[i.id()] = attrs

    cs_single_pr.changeAttributeValues(attrs_dict)
    cs_single.updateFields()
    # sq.layer_to_spatialite(cs_single, 'test3.sqlite', 'scs', 'pk', 'geom')

    cs_split_dict = {'layer': cs_single,
                     'id_col': cs_id_col}
    
    return cs_split_dict
