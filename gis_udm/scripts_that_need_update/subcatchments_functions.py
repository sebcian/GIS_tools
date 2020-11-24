# General Modules
import os
import sys
import shutil
import subprocess
import yaml
import random
import string
import numpy as np
import pandas as pd

# Qgis modules
from pyqgis.initialize_qgis import *
# from processing.tools import dataobjects

# General modules
import pyqgis.general_functions as gn
import pyqgis.sqlite_functions as sq



########################################################################################################################
# FUNCTIONS
########################################################################################################################

# Define extents
def subcatchments_extents(dn_in, systype_col, systype_where, min_length, buffer_distance, secondary_buffer, fill_area):

    # Selection of features to buffer
    selection_ids = []
    for i in dn_in.getFeatures():
        if i.geometry().length() >= min_length and i[systype_col].lower() == systype_where:
            selection_ids.append(i.id())

    if selection_ids:
        dn_in.selectByIds(selection_ids)
        print('Selection was succesfull!!')
    else:
        print('None features selected. Exiting..')
        exit(-1)

    # Apply buffer to selection
    buffered = []
    segments = 5
    end_cap_style = 2
    join_style = 0
    miter_limit = 2
    for i in dn_in.selectedFeatures():
        print('Buffering feature: {}. System type: {}. Length: {}'.format(i.id(), i[systype_col], i.geometry().length()))
        buffered.append(i.geometry().buffer(buffer_distance,
                                            segments,
                                            end_cap_style,
                                            join_style,
                                            miter_limit))

    # Create new memory layer with buffered features
    crs = dn_in.crs().authid()
    dn_buff = QgsVectorLayer('Polygon?crs={}'.format(crs), 'dn_buff', 'memory')
    dn_buff_pr = dn_buff.dataProvider()
    for i in buffered:
        feat = QgsFeature()
        feat.setGeometry(i)
        dn_buff_pr.addFeatures([feat])
    dn_buff.updateExtents()

    # Buffer out and dissolve new layer
    tool = 'native:buffer'
    parameters = {'INPUT': dn_buff,
                  'DISTANCE': secondary_buffer,
                  'SEGMENTS': 5,
                  'END_CAP_STYLE': 0,
                  'JOIN_STYLE': 0,
                  'MITER_LIMIT': 2,
                  'DISSOLVE': True,
                  'OUTPUT': 'memory:'}

    pol_buff_out = processing.run(tool, parameters)['OUTPUT']

    # Buffer in dissolved
    tool = 'native:buffer'
    parameters = {'INPUT': pol_buff_out,
                  'DISTANCE': secondary_buffer*-1.0,
                  'SEGMENTS': 5,
                  'END_CAP_STYLE': 0,
                  'JOIN_STYLE': 0,
                  'MITER_LIMIT': 2,
                  'DISSOLVE': True,
                  'OUTPUT': 'memory:'}

    pol_buff_in = processing.run(tool, parameters)['OUTPUT']

    # Delete holes from polygon
    tool = 'native:deleteholes'
    parameters = {'INPUT': pol_buff_in,
                  'MIN_AREA': fill_area,
                  'OUTPUT': 'memory:'}
    pol_holes = processing.run(tool, parameters)['OUTPUT']

    # Simplify geometries
    tool = 'native:simplifygeometries'
    parameters = {'INPUT': pol_holes,
                  'METHOD': 0,
                  'TOLERANCE': 5,
                  'OUTPUT': 'memory:'}
    pol_sim = processing.run(tool, parameters)['OUTPUT']

    # Multipart to single part
    single = processing.run('native:multiparttosingleparts', parameters={'INPUT': pol_sim,
                                                                         'OUTPUT': 'memory:'})['OUTPUT']
    # Dissolve
    tool = "native:dissolve"
    parameters = {'INPUT': single,
                  'FIELD': [],
                  'OUTPUT': 'memory:'}
    dissolved = processing.run(tool, parameters)['OUTPUT']

    # Multipart to single part
    single = processing.run('native:multiparttosingleparts', parameters={'INPUT': dissolved,
                                                                         'OUTPUT': 'memory:'})['OUTPUT']

    return single


def delineate_thiessen(dn_in, dn_in_id_col, systype_col, systype_where, min_length, split_length, frame_buffer, rmarea, thiessen_extent=None, thiessen_mask=None):

    # Setting
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    # tmp directory
    tmp_dir = os.path.join(os.getcwd(), 'tmp_dir')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Selection of features
    selection_ids = []
    for i in dn_in.getFeatures():
        if i.geometry().length() >= min_length and i[systype_col].lower() == systype_where:
            selection_ids.append(i.id())

    if selection_ids:
        dn_in.selectByIds(selection_ids)
        print('Selection was succesfull!!')
    else:
        print('None features selected. Exiting..')
        exit(-1)

    # Save selected features
    tool = 'qgis:saveselectedfeatures'
    parameters = {'INPUT': dn_in,
                  'OUTPUT': 'memory:'}
    dn_selected = processing.run(tool, parameters)['OUTPUT']
    # sq.layer_to_spatialite(p_along_selection, 'test.sqlite', 'palong_sel', 'pk', 'geom')

    # Split lines in sub-segments
    tool = 'native:splitlinesbylength'
    parameters = {'INPUT': dn_selected,
                  'LENGTH': split_length,
                  'OUTPUT': 'memory:'}
    dn_split = processing.run(tool, parameters)['OUTPUT']

    # Create mid points in split lines
    dn_mid_point = gn.line_mid_point(layer=dn_split)

    # Create Thiessen polygons
    tmp_shp = os.path.join(tmp_dir, 'tmp_shp.shp')
    tool = 'saga:thiessenpolygons'
    parameters = {'POINTS': dn_mid_point,
                  'POLYGONS': tmp_shp,
                  'FRAME': frame_buffer}
    thiessen = processing.run(tool, parameters)['POLYGONS']
    thiessen = gn.read_shp(thiessen)

    # Column names
    names = dn_mid_point.fields().names()
    mapping = []
    for i, j in zip(thiessen.fields(), names):
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': j,
                        'precision': i.precision(),
                        'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': thiessen,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    refactored = processing.run(tool, parameters)['OUTPUT']

    # Clipping
    if thiessen_extent and thiessen_mask:
        # Clip thiessen mask with thiessen extent
        tool = 'native:clip'
        parameters = {'INPUT': thiessen_mask,
                      'OVERLAY': thiessen_extent,
                      'OUTPUT': 'memory:'}
        clip = processing.run(tool, parameters, context=context)['OUTPUT']
        clip_key = True

    elif thiessen_extent and not thiessen_mask:
        clip = thiessen_extent
        clip_key = True

    else:
        clip_key = False

    if clip_key:
        print('Clipping thiessen polygons')
        # Clip thiessen polygons with clip1
        tool = 'native:clip'
        parameters = {'INPUT': refactored,
                      'OVERLAY': clip,
                      'OUTPUT': 'memory:'}
        thiessen_clip = processing.run(tool, parameters, context=context)['OUTPUT']
    else:
        thiessen_clip = thiessen

    # Dissolved by id
    tool = 'native:dissolve'
    parameters = {'INPUT': thiessen_clip,
                  'FIELD': dn_in_id_col,
                  'OUTPUT': 'memory:'}

    dissolved = processing.run(tool, parameters, context=context)['OUTPUT']

    # Cleaning
    cleaned = gn.clean_vector(dissolved, rmarea)

    # Multipart to single part
    single = processing.run('native:multiparttosingleparts', parameters={'INPUT': cleaned,
                                                                         'OUTPUT': 'memory:'})['OUTPUT']

    # Add subcatchment id to attribute table
    single_pr = single.dataProvider()
    single_pr.addAttributes([QgsField('subcatchment_id', QVariant.String, typeName='string')])
    single.updateFields()
    idx = single_pr.fieldNameIndex('subcatchment_id')
    n_digits = len(str(single.featureCount()))
    field_dict = {}
    for i in single.getFeatures():
        num = 'subcatchment_{{:0{}}}'.format(n_digits).format(i.id())
        field_dict[i.id()] = {idx: num}
    single_pr.changeAttributeValues(field_dict)
    single.updateFields()

    # Manage attribute table
    original_names = cleaned.fields().names()
    original_names.append('subcatchment_id')
    mapping = []
    for i, j in zip(single.fields(), original_names):
        if i.name() == dn_in_id_col:
            name = 'connection'
        else:
            name = j
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': name,
                        'precision': i.precision(),
                        'type': i.type()})

    mapping.insert(0, mapping[-1])
    mapping = mapping[0:-1]

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': single,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    refactored = processing.run(tool, parameters)['OUTPUT']

    return refactored


def pre_population(pop_layer, pop_id_col, pop_pop_col, pop_area_col, pop_density_col,
                   ext_layer,
                   sub_layers, sub_id_cols, sub_systype,
                   pop_fracs):

    # Subcatchments overlay
    if len(sub_layers) > 1:
        i = 0
        while i < len(sub_layers) - 1:
            if i == 0:
                a = sub_layers[i]
                a_id_col = [sub_id_cols[i]]
                a_systype = [sub_systype[i]]
                tup_a = [(x, y) for x, y in zip(a_id_col, a_systype)]
            else:
                a = union_single
                a_id_col = a.fields().names()
                a_systype = [x.split('_')[0] for x in a_id_col]
                tup_a = [(x, y) for x, y in zip(a_id_col, a_systype)]

            b = sub_layers[i + 1]
            b_id_col = [sub_id_cols[i + 1]]
            b_systype = [sub_systype[i + 1]]
            tup_b = [(x, y) for x, y in zip(b_id_col, b_systype)]

            # Refactor a
            mapping = []
            for z in tup_a:
                for f in a.fields():
                    if f.name() == z[0]:
                        name = f.name().replace('{}_'.format(z[1]), '')
                        name = '{}_'.format(z[1]) + name
                        mapping.append({'expression': '\"{}\"'.format(f.name()),
                                        'length': f.length(),
                                        'name': name,
                                        'precision': f.precision(),
                                        'type': f.type()})

            tool = 'qgis:refactorfields'
            parameters = {'INPUT': a,
                          'FIELDS_MAPPING': mapping,
                          'OUTPUT': 'memory:'}
            a_refac = processing.run(tool, parameters)['OUTPUT']

            # Refactor b
            mapping = []
            for z in tup_b:
                for f in b.fields():
                    if f.name() == z[0]:
                        name = f.name().replace('{}_'.format(z[1]), '')
                        name = '{}_'.format(z[1]) + name
                        mapping.append({'expression': '\"{}\"'.format(f.name()),
                                        'length': f.length(),
                                        'name': name,
                                        'precision': f.precision(),
                                        'type': f.type()})

            tool = 'qgis:refactorfields'
            parameters = {'INPUT': b,
                          'FIELDS_MAPPING': mapping,
                          'OUTPUT': 'memory:'}
            b_refac = processing.run(tool, parameters)['OUTPUT']

            # Union
            tool = 'native:union'
            parameters = {'INPUT': a_refac,
                          'OVERLAY': b_refac,
                          'OUTPUT': 'memory:'}
            overlay_union = processing.run(tool, parameters)['OUTPUT']

            # Multipart to single part
            union_single = processing.run('native:multiparttosingleparts', parameters={'INPUT': overlay_union,
                                                                                       'OUTPUT': 'memory:'})['OUTPUT']

            i += 1
    else:
        # Refactor
        a = sub_layers[0]
        a_id_col = [sub_id_cols[0]]
        a_systype = [sub_systype[0]]
        tup_a = [(x, y) for x, y in zip(a_id_col, a_systype)]

        mapping = []
        for z in tup_a:
            for f in a.fields():
                if f.name() == z[0]:
                    name = f.name().replace('{}_'.format(z[1]), '')
                    name = '{}_'.format(z[1]) + name
                    mapping.append({'expression': '\"{}\"'.format(f.name()),
                                    'length': f.length(),
                                    'name': name,
                                    'precision': f.precision(),
                                    'type': f.type()})

        tool = 'qgis:refactorfields'
        parameters = {'INPUT': a,
                      'FIELDS_MAPPING': mapping,
                      'OUTPUT': 'memory:'}
        refac = processing.run(tool, parameters)['OUTPUT']
        union_single = refac


    # Clip population map with total extents and updates area and population fields
    tool = 'native:clip'
    parameters = {'INPUT': pop_layer,
                  'OVERLAY': ext_layer,
                  'OUTPUT': 'memory:'}
    pop_clip = processing.run(tool, parameters)['OUTPUT']

    # Update area and population information on pop_clip
    pop_clip_pr = pop_clip.dataProvider()
    idx_area = pop_clip_pr.fieldNameIndex(pop_area_col)
    idx_pop = pop_clip_pr.fieldNameIndex(pop_pop_col)
    attrs = {}
    for i in pop_clip.getFeatures():
        new_area = i.geometry().area()
        pden = np.nan if str(i[pop_density_col]) == 'NULL' else float(i[pop_density_col])
        new_pop = new_area * pden
        attrs[i.id()] = {idx_area: new_area, idx_pop: new_pop}
    pop_clip_pr.changeAttributeValues(attrs)
    pop_clip.updateFields()

    # Refactor population
    mapping = []
    for i in pop_clip.fields():
        if i.name() in [pop_id_col, pop_pop_col]:
            mapping.append({'expression': '\"{}\"'.format(i.name()),
                            'length': i.length(),
                            'name': i.name(),
                            'precision': i.precision(),
                            'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': pop_clip,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    pop_refac = processing.run(tool, parameters)['OUTPUT']

    # Intersect subcatchments union with population layer
    tool = 'native:intersection'
    parameters = {'INPUT': union_single,
                  'OVERLAY': pop_refac,
                  'INPUT_FIELDS': '',
                  'OVERLAY_FIELDS': '',
                  'OUTPUT': 'memory:'}

    intersection = processing.run(tool, parameters)['OUTPUT']

    # Multipart to single part
    intersection_single = processing.run('native:multiparttosingleparts', parameters={'INPUT': intersection,
                                                                               'OUTPUT': 'memory:'})['OUTPUT']

    # Update area in intersection single
    intersection_single_pr = intersection_single.dataProvider()
    intersection_single_pr.addAttributes([QgsField('area', QVariant.Double, typeName='double')])
    intersection_single.updateFields()
    idx_area = intersection_single_pr.fieldNameIndex('area')

    attrs = {}
    for i in intersection_single.getFeatures():
        new_area = i.geometry().area()
        attrs[i.id()] = {idx_area: new_area}
    intersection_single_pr.changeAttributeValues(attrs)
    intersection_single.updateFields()

    # Dataframe
    df = gn.layer_df(intersection_single)
    grouper = df.groupby(pop_id_col)
    df2 = grouper.agg({pop_pop_col: 'mean', 'area': 'sum'})
    df2['density'] = df2[pop_pop_col] / df2['area']
    for idx in df2.index:
        df.loc[df[pop_id_col] == idx, 'density'] = df2.loc[idx, 'density']

    df[pop_pop_col] = df['area'] * df['density']

    # Add new density and new population
    intersection_single_pr = intersection_single.dataProvider()
    intersection_single_pr.addAttributes([QgsField('density', QVariant.Double, typeName='double')])
    idx_den = intersection_single_pr.fieldNameIndex('density')
    idx_pop = intersection_single_pr.fieldNameIndex(pop_pop_col)
    attrs = {}
    for i in intersection_single.getFeatures():
        den = df.loc[i.id() - 1, 'density']
        pop = df.loc[i.id() - 1, pop_pop_col]
        attrs[i.id()] = {idx_den: den, idx_pop: pop}
    intersection_single_pr.changeAttributeValues(attrs)
    intersection_single.updateFields()

    # Population fractions
    aux_ = [i for i in df.columns if i not in [pop_id_col, pop_pop_col, 'area', 'density']]
    fracs_cols = ['{}_frac'.format(i) for i in sub_systype]
    fracs = pd.DataFrame(index=df.index, columns=fracs_cols)
    for i, j, k in zip(aux_, fracs_cols, pop_fracs):
        fracs[j][~pd.isnull(df[i])] = k
    fracs['SUM'] = fracs.loc[:, fracs_cols].sum(axis=1)
    f_fracs = fracs[fracs_cols].divide(fracs['SUM'], axis=0)

    # Add population fractions to layer
    intersection_single_pr = intersection_single.dataProvider()
    for i in fracs_cols:
        intersection_single_pr.addAttributes([QgsField(i, QVariant.Double, typeName='double')])
    intersection_single.updateFields()

    idx = [intersection_single_pr.fieldNameIndex(i) for i in fracs_cols]
    attrs_dict = {}
    for i in intersection_single.getFeatures():
        attr = {}
        for j, k in zip(idx, fracs_cols):
            frac = f_fracs.loc[i.id() - 1, k]
            attr[j] = frac
        attrs_dict[i.id()] = attr
    intersection_single_pr.changeAttributeValues(attrs_dict)
    intersection_single.updateFields()

    pre_population_dict = {'layer': intersection_single,
                           'sub_systypes': sub_systype,
                           'sub_id_cols': aux_,
                           'population_id_col': pop_id_col,
                           'population_col': pop_pop_col,
                           'area_col': 'area',
                           'density_col': 'density',
                           'fraction_cols': fracs_cols}

    return pre_population_dict


def pre_ato(sub_layers, sub_id_cols, sub_systype,
            lc_layer, lc_value_col,
            ato_fracs, ato_ignore_values):

    # Setting
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    # tmp directory
    tmp_dir = os.path.join(os.getcwd(), 'tmp_dir')
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Subcatchments overlay
    if len(sub_layers) > 1:
        i = 0
        while i < len(sub_layers) - 1:
            if i == 0:
                a = sub_layers[i]
                a_id_col = [sub_id_cols[i]]
                a_systype = [sub_systype[i]]
                tup_a = [(x, y) for x, y in zip(a_id_col, a_systype)]
            else:
                a = union_single
                a_id_col = a.fields().names()
                a_systype = [x.split('_')[0] for x in a_id_col]
                tup_a = [(x, y) for x, y in zip(a_id_col, a_systype)]

            b = sub_layers[i + 1]
            b_id_col = [sub_id_cols[i + 1]]
            b_systype = [sub_systype[i + 1]]
            tup_b = [(x, y) for x, y in zip(b_id_col, b_systype)]

            # Refactor a
            mapping = []
            for z in tup_a:
                for f in a.fields():
                    if f.name() == z[0]:
                        name = f.name().replace('{}_'.format(z[1]), '')
                        name = '{}_'.format(z[1]) + name
                        mapping.append({'expression': '\"{}\"'.format(f.name()),
                                        'length': f.length(),
                                        'name': name,
                                        'precision': f.precision(),
                                        'type': f.type()})

            tool = 'qgis:refactorfields'
            parameters = {'INPUT': a,
                          'FIELDS_MAPPING': mapping,
                          'OUTPUT': 'memory:'}
            a_refac = processing.run(tool, parameters)['OUTPUT']

            # Refactor b
            mapping = []
            for z in tup_b:
                for f in b.fields():
                    if f.name() == z[0]:
                        name = f.name().replace('{}_'.format(z[1]), '')
                        name = '{}_'.format(z[1]) + name
                        mapping.append({'expression': '\"{}\"'.format(f.name()),
                                        'length': f.length(),
                                        'name': name,
                                        'precision': f.precision(),
                                        'type': f.type()})

            tool = 'qgis:refactorfields'
            parameters = {'INPUT': b,
                          'FIELDS_MAPPING': mapping,
                          'OUTPUT': 'memory:'}
            b_refac = processing.run(tool, parameters)['OUTPUT']

            # Union
            tool = 'native:union'
            parameters = {'INPUT': a_refac,
                          'OVERLAY': b_refac,
                          'OUTPUT': 'memory:'}
            overlay_union = processing.run(tool, parameters)['OUTPUT']

            # Multipart to single part
            union_single = processing.run('native:multiparttosingleparts', parameters={'INPUT': overlay_union,
                                                                                       'OUTPUT': 'memory:'})['OUTPUT']

            i += 1
    else:
        union_single = sub_layers[0]

    # Refactor land cover layer
    mapping = []
    for i in lc_layer.fields():
        if i.name() in [lc_value_col]:
            mapping.append({'expression': '\"{}\"'.format(i.name()),
                            'length': i.length(),
                            'name': i.name(),
                            'precision': i.precision(),
                            'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': lc_layer,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}

    lc_refac = processing.run(tool, parameters, context=context)['OUTPUT']

    # Intersect subcatchments union with land cover layer
    tmp_shp = os.path.join(tmp_dir, 'tmp_shp.shp')
    tool = 'saga:intersect'
    parameters = {'A': union_single,
                  'B': lc_refac,
                  'SPLIT': True,
                  'RESULT': tmp_shp}
    intersection = processing.run(tool, parameters, context=context)['RESULT']
    intersection = gn.read_shp(intersection)

    # Column names
    names = union_single.fields().names()
    for i in lc_refac.fields().names():
        names.append(i)
    mapping = []
    for i, j in zip(intersection.fields(), names):
        mapping.append({'expression': '\"{}\"'.format(i.name()),
                        'length': i.length(),
                        'name': j,
                        'precision': i.precision(),
                        'type': i.type()})

    tool = 'qgis:refactorfields'
    parameters = {'INPUT': intersection,
                  'FIELDS_MAPPING': mapping,
                  'OUTPUT': 'memory:'}
    int_refact = processing.run(tool, parameters, context=context)['OUTPUT']
    # sq.layer_to_spatialite(int_refact, 'test.sqlite', 'union', 'pk', 'geom')

    # Multipart to single part
    int_single = processing.run('native:multiparttosingleparts', parameters={'INPUT': int_refact,
                                                                             'OUTPUT': 'memory:'})['OUTPUT']
    # sq.layer_to_spatialite(int_refact, 'test.sqlite', 'int_single', 'pk', 'geom')

    # Remove features with value = "ignore_ato_value"
    del_ = []
    for i in int_single.getFeatures():
        if i[lc_value_col] in ato_ignore_values:
            del_.append(i.id())

    int_single_pr = int_single.dataProvider()
    int_single_pr.deleteFeatures(del_)
    int_single.updateExtents()

    # Redo layer (to get ordered feature ids)
    crs = int_single.crs().authid()
    int_single_attrs = int_single.dataProvider().fields().toList()
    int_single_df = gn.layer_df(int_single)

    pre_ato = QgsVectorLayer('Polygon?crs={}'.format(crs), 'pre_ato', 'memory')
    pre_ato_pr = pre_ato.dataProvider()

    # Fields
    pre_ato_pr.addAttributes(int_single_attrs)
    pre_ato.updateFields()

    # Features
    c = 0
    for i in int_single.getFeatures():
        feat = QgsFeature()
        feat.setGeometry(i.geometry())
        feat.setAttributes(int_single_df.loc[c, :].tolist())
        pre_ato_pr.addFeatures([feat])
        c += 1
    pre_ato.updateExtents()
    # sq.layer_to_spatialite(pre_ato, 'test.sqlite', 'test','pk', 'geom')

    # Update area in intersection single
    pre_ato_pr = pre_ato.dataProvider()
    pre_ato_pr.addAttributes([QgsField('area', QVariant.Double, typeName='double')])
    pre_ato.updateFields()
    idx_area = pre_ato_pr.fieldNameIndex('area')

    attrs = {}
    for i in pre_ato.getFeatures():
        new_area = i.geometry().area()
        attrs[i.id()] = {idx_area: new_area}
    pre_ato_pr.changeAttributeValues(attrs)
    pre_ato.updateFields()

    # Dataframe
    df = gn.layer_df(pre_ato)

    # ATO fractions
    aux_ = [i for i in df.columns if i not in [lc_value_col, 'area']]
    fracs_cols = ['{}_frac'.format(i) for i in sub_systype]
    fracs = pd.DataFrame(index=df.index, columns=fracs_cols)
    for i, j, k in zip(aux_, fracs_cols, ato_fracs):
        fracs[j][~pd.isnull(df[i])] = k
    fracs['SUM'] = fracs.loc[:, fracs_cols].sum(axis=1)
    f_fracs = fracs[fracs_cols].divide(fracs['SUM'], axis=0)

    # Add ato fractions to layer
    pre_ato_pr = pre_ato.dataProvider()
    for i in fracs_cols:
        pre_ato_pr.addAttributes([QgsField(i, QVariant.Double, typeName='double')])
    pre_ato.updateFields()

    idx = [pre_ato_pr.fieldNameIndex(i) for i in fracs_cols]
    attrs_dict = {}
    for i in pre_ato.getFeatures():
        attr = {}
        for j, k in zip(idx, fracs_cols):
            frac = f_fracs.loc[i.id() - 1, k]
            attr[j] = frac
        attrs_dict[i.id()] = attr
    pre_ato_pr.changeAttributeValues(attrs_dict)
    pre_ato.updateFields()
    # sq.layer_to_spatialite(pre_ato, 'test.sqlite', 'int', 'pk', 'geom')

    pre_ato_dict = {'layer': pre_ato,
                    'sub_systypes': sub_systype,
                    'sub_id_cols': aux_,
                    'land_cover_value_id_col': lc_value_col,
                    'fraction_cols': fracs_cols}

    return pre_ato_dict


def population(pre_pop_df, pre_pop_sub_id_cols, pre_pop_pop_col, pre_pop_fraction_cols, pre_pop_sub_systypes,
               sub_layers, sub_id_col):

    out = {}
    for i, j, k, l, m in zip(sub_layers, sub_id_col, pre_pop_sub_id_cols, pre_pop_fraction_cols, pre_pop_sub_systypes):


        # Get population per subcatchment id
        df = pre_pop_df.loc[~pd.isnull(pre_pop_df[k]), [k, pre_pop_pop_col, l]].copy()
        df['out'] = np.ceil(df[pre_pop_pop_col].astype(float) * df[l].astype(float))
        grouper = df.groupby(k)
        pop_df = grouper.agg({'out': 'sum'})

        # Add population to layer df
        lyr_df = gn.layer_df(i)
        lyr_df['population'] = np.nan
        for x in pop_df.index:
            lyr_df.loc[lyr_df[j] == x, 'population'] = pop_df.loc[x, 'out']

        # Update attribute table
        pr = i.dataProvider()
        if pr.fieldNameIndex('population') == -1:
            pr.addAttributes([QgsField('population', QVariant.Int, typeName='integer')])
        i.updateFields()

        idx = pr.fieldNameIndex('population')
        attrs_dict = {}
        for feat in i.getFeatures():
            attrs_dict[feat.id()] = {idx: lyr_df.loc[feat.id() - 1, 'population']}
        pr.changeAttributeValues(attrs_dict)
        i.updateFields()

        out[m] = i

    return out


def ato(pre_ato_df, pre_ato_sub_systypes, pre_ato_sub_id_cols, pre_ato_value_col, pre_ato_area_col, pre_ato_fraction_cols,
        sub_layers, sub_id_col):

    out = {}
    for i, j, k, l, m in zip(sub_layers, sub_id_col, pre_ato_sub_id_cols, pre_ato_fraction_cols, pre_ato_sub_systypes):

        # Multiply: area * fraction
        df = pre_ato_df.loc[~pd.isnull(pre_ato_df[k]), [k, pre_ato_value_col, pre_ato_area_col, l]].copy()
        df['out'] = df[pre_ato_area_col].astype(float) * df[l].astype(float) / 10000.0

        # Get individual runoff surfaces
        grouper = df.groupby([k, pre_ato_value_col])
        runoff_surf = grouper.agg({'out': 'sum'}).unstack()
        runoff_surf.columns = ['runoff_surface_{:0>2d}'.format(i) for i in runoff_surf.columns.get_level_values(1)]
        runoff_surf['contributing_area'] = runoff_surf.sum(axis=1)

        # Add new columns to layer df
        lyr_df = gn.layer_df(i)
        lyr_df = lyr_df.reindex(columns= lyr_df.columns.tolist() + runoff_surf.columns.tolist())
        for x in runoff_surf.index.values:
            lyr_df.loc[lyr_df[j] == x, runoff_surf.columns.tolist()] = runoff_surf.loc[x, runoff_surf.columns.tolist()].values

        # Update attribute table
        pr = i.dataProvider()
        for x in runoff_surf.columns:
            if pr.fieldNameIndex(x) == -1:
                pr.addAttributes([QgsField(x, QVariant.Double, typeName='double')])
        i.updateFields()

        idx = [pr.fieldNameIndex(i) for i in runoff_surf.columns]
        attrs_dict = {}
        for feat in i.getFeatures():
            attr = {}
            for x, y in zip(idx, runoff_surf.columns):
                val = lyr_df.loc[feat.id() - 1, y]
                attr[x] = val
            attrs_dict[feat.id()] = attr
        pr.changeAttributeValues(attrs_dict)
        i.updateFields()

        out[m] = i

    return out


def set_slope(layer, slope_rast):

    tool = 'qgis:zonalstatistics'
    parameters = {'INPUT_RASTER': slope_rast,
                  'RASTER_BAND': 1,
                  'INPUT_VECTOR': layer,
                  'COLUMN_PREFIX': 'slope_',
                  'STATS': 2}
    lyr_slope = processing.run(tool, parameters)['INPUT_VECTOR']

    lyr_slope_dict = {'layer': lyr_slope,
                      'slope_col': 'slope_mean'}

    return lyr_slope_dict
