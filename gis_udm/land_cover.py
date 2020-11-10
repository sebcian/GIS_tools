'''
DESCRIPTION / INSTRUCTIONS

'''
# Modules
import os
import shutil
import tempfile
import shapely
from shapely import wkt
import geopandas as gpd

import pygis
import wrappers.grass_session as grass
import wrappers.qgis_session as qgs

# Class
class LandCover(object):
    """
    INPUT:
    
    input_gpkg: Geopackage of input layers
    output_gpkg: Geopackage for output
    crs: Coordinate reference system
    
    """
    
    def __init__(self, input_gpkg, output_gpkg, crs):

        self.tmpdir = tempfile.mkdtemp(prefix='tmp_lc_')
        
        self.input_gpkg = input_gpkg
        self.output_gpkg = output_gpkg

        print('GRASS GIS session:\n')
        gisdb = os.path.join(self.tmpdir, 'grassdata')
        self.gs = grass.GrassSession(gisdb=gisdb, location='lct', crs=crs)
        print('\n\n')
        print('QGIS session:\n')
        qgs.qgis_version()
        

    def preprocess_boundary(self, layer, simplify):
        db = self.input_gpkg

        print('\n\n\nIMPORT BOUNDARY LAYER STARTED')

        print('\n\nImporting LAYER into GRASS GIS')
        self.gs.import_interactive(input_path=db, layer=layer, snap=-1, output='tmp_imp')

        print('\n\nRemoving any existing inner-ring in BOUNDARY')
        self.gs.run_command('v.centroids', input='tmp_imp', output='tmp_centroids', overwrite=True)
        self.gs.run_command('v.db.addcolumn', map='tmp_centroids', columns='\"diss INT\"')
        self.gs.run_command('v.db.update', map='tmp_centroids', column='diss', value=1)
        self.gs.run_command('v.dissolve', input='tmp_centroids', column='diss', output='tmp_diss',
                               overwrite=True)

        print('\n\nSimplifying BOUNDARY geometries')
        self.gs.run_command('v.generalize', overwrite=True, input='tmp_diss',
                               output='tmp_simp', method='douglas', threshold=simplify)

        print('\n\nManaging BOUNDARY attribute table')
        self.gs.run_command('v.category', overwrite=True, input='tmp_simp',
                               output='tmp_catdel', option='del', cat='-1', layer='1')
        self.gs.run_command('v.category', overwrite=True, input='tmp_catdel',
                               output=layer, option='add')
        # self.gs.run_command('v.db.droptable', flags='f', map=self.boundary_preprocessed)
        self.gs.run_command('v.db.addtable', map=layer, columns='\"id TEXT\"')
        self.gs.run_command('v.db.update', map=layer, column='id', value='boundary')

        print('\n\nRemoving temporal files')
        self.gs.run_command('g.remove', flags='f', type='vector', pattern='tmp_*')

    
    def preprocess_land_cover_layer(self, layer, priority, boundary, seglen, snap, simplify, area_thr):
        db = self.input_gpkg
        bnd = boundary
        tmp_db = os.path.join(self.tmpdir, 'tmp.gpkg')

        print('\n\tImporting {} into GRASS GIS'.format(layer))
        self.gs.import_interactive(input_path=db, layer=layer, snap=-1, output='tmp_imp')

        print('CLIPPING layer with boundary')
        self.gs.clip_layer(layer_a='tmp_imp', layer_b=bnd, operator='and', output='tmp_clip')

        print('DENSIFYING layer features')
        self.gs.run_command('v.split', flags='n', input='tmp_clip', output='tmp_split',
                            length=seglen, overwrite=True)

        print('EXPORTING layer to temporary database')
        flags = 'su' if os.path.exists(tmp_db) else 's'
        self.gs.run_command('v.out.ogr', flags=flags, input='tmp_split', output=tmp_db,
                            output_layer='tmp_split', format='GPKG', overwrite=True)

        print('SELF-SNAPPING layer')
        input_layer = '{}|layername=tmp_split'.format(tmp_db)
        reference_layer = '{}|layername=tmp_split'.format(tmp_db)
        out_ = 'ogr:dbname=\'{}\' table=\"tmp_snap\" (geom) sql='.format(tmp_db)
        qgs.snap_geometries(input_layer, reference_layer, snap, out_)
        
        print('IMPORTING self-snapped layer')
        self.gs.import_interactive(input_path=tmp_db, layer='tmp_snap', snap=-1, output='tmp_imp')

        print('SIMPLIFYING layer geometries')
        self.gs.run_command('v.generalize', overwrite=True, input='tmp_imp', output='tmp_simp',
                            method='douglas', threshold=simplify)

        print('CLEANING layer')
        self.gs.run_command('v.clean', flags='c', input='tmp_simp',
                            output='tmp_clean',
                            tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
                            threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True)

        print('Managing ATTRIBUTE TABLE')
        self.gs.run_command('v.category', overwrite=True, input='tmp_clean',
                            output='tmp_catdel', option='del', cat='-1', layer='1')
        self.gs.run_command('v.category', overwrite=True, input='tmp_catdel',
                            output=layer, option='add')
        self.gs.run_command('v.db.droptable', flags='f', map=layer)
        self.gs.run_command('v.db.addtable', map=layer, columns='\"priority INT, class TEXT\"')
        self.gs.run_command('v.db.update', map=layer, column='priority', value=priority)
        self.gs.run_command('v.db.update', map=layer, column='class', value=layer)

        print('\n\nRemoving temporal files')
        self.gs.run_command('g.remove', flags='f', type='vector', pattern='tmp_*')
        if os.path.exists(tmp_db):
            os.remove(tmp_db)
    

    def generate_land_cover_interactive(
        self, layers, priorities, boundary, seglen, area_thr, snap, simplify, output
    ):

        tmp_dbase = os.path.join(self.tmpdir, 'generate_land_cover_interactive.gpkg')

        # Sorting layers tuples by priority
        l_ = [(i, j) for i, j in zip(priorities, layers)]
        l_ = sorted(l_, key=lambda x: x[0])
        priorities = [i[0] for i in l_]
        layers = [i[1] for i in l_]

        print('\n\n\nGENERATE LAND COVER PROCESS STARTED')
        i = 0
        loop_n = 1
        ref = layers[i]
        snp = layers[i + 1]

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('\n\nEXPORTING layer to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command('v.out.ogr', flags=flags, input=boundary, output=tmp_dbase,
                            output_layer=boundary, format='GPKG', overwrite=True)

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('DENSIFYING initial reference layer')
        self.gs.run_command('v.split', flags='n', input=ref, output='tmp_split',
                            length=seglen, overwrite=True)

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('EXPORTING initial reference to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command('v.out.ogr', flags=flags, input='tmp_split', output=tmp_dbase,
                            output_layer='tmp_split', format='GPKG', overwrite=True)

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('SNAPPING initial reference to boundary')
        input_layer = '{}|layername={}'.format(tmp_dbase, 'tmp_split')
        reference_layer = '{}|layername={}'.format(tmp_dbase, boundary)
        out_ = 'ogr:dbname=\'{}\' table=\"tmp_snap\" (geom) sql='.format(tmp_dbase)
        qgs.snap_geometries(input_layer, reference_layer, snap, out_)

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('IMPORTING snapped initial reference into GRASS GIS')
        self.gs.import_interactive(input_path=tmp_dbase, layer='tmp_snap',
                                   snap=-1, output='tmp_snap')

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('CLEANING initial reference layer')
        self.gs.run_command('v.clean', flags='c', input='tmp_snap', output='tmp_ref_cleaned',
                            tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
                            threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True)

        print('\n\n\t\tWorking with <<[ {} ]>>'.format(ref))
        print('EXPORTING cleaned reference layer to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command('v.out.ogr', flags=flags, input='tmp_ref_cleaned', output=tmp_dbase,
                            output_layer='tmp_ref_cleaned', format='GPKG', overwrite=True)

        while i < len(layers):

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('CLIPPING layer with reference')
            self.gs.clip_layer(layer_a=snp, layer_b='tmp_ref_cleaned', operator='not', output='tmp_clip')

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('DENSIFYING layer')
            self.gs.run_command('v.split', flags='n', input='tmp_clip', output='tmp_split',
                                length=seglen, overwrite=True)

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('EXPORTING layer to auxiliary database')
            flags = 'su' if os.path.exists(tmp_dbase) else 's'
            self.gs.run_command(
                'v.out.ogr', flags=flags, input='tmp_split', output=tmp_dbase,
                output_layer='tmp_split', format='GPKG', overwrite=True
            )

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('SNAPPING layer to reference')
            input_layer = '{}|layername={}'.format(tmp_dbase, 'tmp_split')
            reference_layer = '{}|layername=tmp_ref_cleaned'.format(tmp_dbase)
            out_ = 'ogr:dbname=\'{}\' table=\"tmp_snp1\" (geom) sql='.format(tmp_dbase)
            qgs.snap_geometries(input_layer, reference_layer, snap, out_)

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('SNAPPING layer to boundary')
            input_layer = '{}|layername=tmp_snp1'.format(tmp_dbase)
            reference_layer = '{}|layername={}'.format(tmp_dbase, boundary)
            out_ = 'ogr:dbname=\'{}\' table=\"tmp_snp2\" (geom) sql='.format(tmp_dbase)
            qgs.snap_geometries(input_layer, reference_layer, snap, out_)

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('IMPORTING layer into GRASS GIS')
            self.gs.import_interactive(input_path=tmp_dbase, layer='tmp_snp2', snap=-1, output='tmp_snp')

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('CLEANING layer')
            self.gs.run_command(
                'v.clean', flags='c', input='tmp_snp', output='tmp_snp_cleaned',
                tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
                threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True
            )

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('MERGING layer with reference')
            self.gs.merge_layers(layers=['tmp_ref_cleaned', 'tmp_snp_cleaned'], output='tmp_merged')

            print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(snp, loop_n, len(layers)-1))
            print('CLEANING (export/import) layer')
            flags = 'su' if os.path.exists(tmp_dbase) else 's'
            self.gs.run_command(
                'v.out.ogr', flags=flags, input='tmp_merged', output=tmp_dbase,
                output_layer='tmp_merged', format='GPKG', overwrite=True
            )
            self.gs.import_interactive(input_path=tmp_dbase, layer='tmp_merged', snap=-1, output='tmp_ref0')

            print('\n\n\t\tWorking with <<[ merged layer ]>>\t[loop {}/{}]'.format(loop_n, len(layers)-1))
            print('CLEANING layer')
            self.gs.run_command(
                'v.clean', flags='c', input='tmp_ref0', output='tmp_ref_cleaned',
                tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
                threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True
            )

            new_ref_name = '-'.join(layers[0:loop_n+1])
            print('\n\n\n\t\tNew REFERENCE layer <<[ {} ]>>\t[loop {}/{}]'.format(new_ref_name, loop_n, len(layers)-1))

            i += 2 if i == 0 else +1
            loop_n += 1
            if i < len(layers):
                snp = layers[i]

                print('\n\n\t\tWorking with <<[ {} ]>>\t[loop {}/{}]'.format(new_ref_name, loop_n, len(layers)-1))
                print('EXPORTING new reference to auxiliary database')
                flags = 'su' if os.path.exists(tmp_dbase) else 's'
                self.gs.run_command(
                    'v.out.ogr', flags=flags, input='tmp_ref_cleaned', output=tmp_dbase,
                    output_layer='tmp_ref_cleaned', format='GPKG', overwrite=True
                )

        print('\n\n\t\tLAND COVER LAYER CREATED')

        print('\n\nCLIPPING land cover with boundary')
        self.gs.clip_layer(layer_a='tmp_ref_cleaned', layer_b=boundary, operator='and', output='tmp_clip')

        print('\n\nSIMPLIFYING land cover layer')
        self.gs.run_command(
            'v.generalize', overwrite=True, input='tmp_clip', output='tmp_simplified',
            method='douglas', threshold=simplify
        )

        print('\n\nCLEANING lad cover geometries')
        self.gs.run_command(
            'v.clean', flags='c', input='tmp_simplified', output=output,
            tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
            threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True
        )

        print('\n\nMANAGING land cover attribute table')
        self.gs.update_table(output)

        print('\n\nREMOVING temporal files')
        self.gs.run_command('g.remove', flags='f', type='vector,raster', pattern='tmp_*')
        if os.path.exists(tmp_dbase):
            os.remove(tmp_dbase)


    def generate_land_cover_static(
        self, layers, priorities, boundary,
        resolution, area_thr, snap, simplify, seglen, output
    ):

        tmp_dbase = os.path.join(self.tmpdir, 'generate_land_cover_static.gpkg')

        # Sorting layers tuples by priority
        l_ = [(i, j) for i, j in zip(priorities, layers)]
        l_ = sorted(l_, key=lambda x: x[0])
        priorities = [i[0] for i in l_]
        layers = [i[1] for i in l_]

        print('\n\n\nGENERATE LAND COVER PROCESS STARTED')

        print('\n\nGetting and Setting REGION and RESOLUTION')
        extent = {'north': [], 'south': [], 'east': [], 'west': [], 'top': [], 'bottom': []}
        for lyr in layers:
            info = self.gs.read_command('v.info', flags='g', map=lyr).strip().split(os.linesep)
            ext = [i.split('=')[1] for i in info]
            for i, j in zip(extent.keys(), ext):
                extent[i].append(j)

        print(
            self.gs.read_command(
                'g.region', flags='p',
                n=max(extent['north']), s=min(extent['south']),
                e=max(extent['east']), w=min(extent['west']),
                res=resolution
            )
        )

        print('\n\nConverting LAYERS into RASTER')
        for lyr in layers:
            print('\n\n\t\t << [{}] >> to RASTER'.format(lyr))
            self.gs.run_command(
                'v.to.rast', overwrite=True, input=lyr, output='tmp_{}'.format(lyr),
                use='attr', attribute_column='priority', label_column='class'
            )

        print('\n\nPatching RASTER LAYERS')
        r_layers = ['tmp_{}'.format(i) for i in layers]
        self.gs.run_command('r.patch', input=','.join(r_layers), output='tmp_patched', overwrite=True)

        print('\n\nConverting PATCHED RASTER into VECTOR')
        self.gs.run_command(
            'r.to.vect', flags='s', overwrite=True, input='tmp_patched',
            output='tmp_patched', type='area', column='priority'
        )

        print('\n\nRemoving PATCHED VECTOR areas smaller than \t\t<< [{}] >> units'.format(area_thr))
        areas = self.gs.read_command(
            'v.to.db', flags='p', map='tmp_patched', option='area',
            type='centroid,boundary'
        )
        areas_list = [i.split('|') for i in areas.strip().split(os.linesep)[1::]]
        cat_keep = [i[0] for i in areas_list if float(i[1]) > area_thr]
        cat_file = os.path.join(self.tmpdir, 'extract_cats.txt')
        with open(cat_file, 'w') as f:
            f.writelines('\n'.join(cat_keep))

        self.gs.run_command(
            'v.extract', input='tmp_patched', output='tmp_extract',
            file=cat_file, overwrite=True
        )

        self.gs.run_command(
            'v.clean', flags='c', input='tmp_extract', output='tmp_rmarea',
            tool='rmarea', threshold=area_thr*1.05, overwrite=True
        )

        print('\n\nSimplifying PATCHED VECTOR geometries')
        self.gs.run_command(
            'v.generalize', overwrite=True, input='tmp_rmarea', output='tmp_simplified',
            method='douglas', threshold=simplify
        )

        print('\n\nExporting PATCHED VECTOR to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command(
            'v.out.ogr', flags=flags, input='tmp_simplified', output=tmp_dbase,
            output_layer='tmp_simplified', format='GPKG', overwrite=True
        )

        print('\n\nExporting BOUNDARY to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command(
            'v.out.ogr', flags=flags, input=boundary, output=tmp_dbase,
            output_layer=boundary, format='GPKG', overwrite=True
        )

        # Snapping process
        print('\n\nReading PATCHED VECTOR with Geopandas')
        patch = gpd.read_file(tmp_dbase, layer='tmp_simplified', driver='GPKG')

        i = 0
        loop_n = 1

        ref_name = layers[i]
        print('\n\n\t\t<<[ {} ]>>'.format(ref_name))
        print('Extracting INITIAL REFERENCE LAYER')
        ref = patch.loc[patch['priority'] == priorities[i], :].reset_index(drop=True).copy()
        ref['geometry'] = ref['geometry'].map(lambda x: pygis.segmentize(x, seglen))
        ref.to_file(tmp_dbase, layer=ref_name, driver='GPKG')

        print('\n\n\t\t<<[ {} ]>>'.format(ref_name))
        print('Snapping INITIAL REFERENCE to BOUNDARY')
        input_layer = '{}|layername={}'.format(tmp_dbase, ref_name)
        reference_layer = '{}|layername={}'.format(tmp_dbase, boundary)
        out_ = 'ogr:dbname=\'{}\' table=\"tmp_ref\" (geom) sql='.format(tmp_dbase)
        qgs.snap_geometries(input_layer, reference_layer, snap, out_)

        print('\n\n\t\t<<[ {} ]>>'.format(ref_name))
        print('Importing INITIAL REFERENCE into GRASS GIS')
        self.gs.run_command('v.in.ogr', overwrite=True, input=tmp_dbase, layer='tmp_ref', output='tmp_ref')

        print('\n\n\t\t<<[ {} ]>>'.format(ref_name))
        print('Converting INITIAL REFERENCE to RASTER')
        self.gs.run_command(
            'v.to.rast', overwrite=True, input='tmp_ref', output='tmp_ref',
            use='attr', attribute_column='priority'
        )

        ref_name = 'tmp_ref'
        snp_name = layers[i + 1]

        print('\n\n\t\t<<[ {} ]>>'.format(snp_name))
        print('Extracting INITIAL SNAP LAYER')
        snp = patch.loc[patch['priority'] == priorities[i+1], :].reset_index(drop=True).copy()
        snp['geometry'] = snp['geometry'].map(lambda x: pygis.segmentize(x, seglen))
        snp.to_file(tmp_dbase, layer=snp_name, driver='GPKG')

        while i < len(layers):

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(snp_name, loop_n, len(layers)-1))
            print('Snapping LAYER to REFERENCE')
            input_layer = '{}|layername={}'.format(tmp_dbase, snp_name)
            reference_layer = '{}|layername={}'.format(tmp_dbase, ref_name)
            out_ = 'ogr:dbname=\'{}\' table=\"tmp_snp\" (geom) sql='.format(tmp_dbase)
            qgs.snap_geometries(input_layer, reference_layer, snap, out_)

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(snp_name, loop_n, len(layers)-1))
            print('Snapping to BOUNDARY')
            input_layer = '{}|layername=tmp_snp'.format(tmp_dbase)
            reference_layer = '{}|layername={}'.format(tmp_dbase, boundary)
            out_ = 'ogr:dbname=\'{}\' table=\"tmp_snp2\" (geom) sql='.format(tmp_dbase)
            qgs.snap_geometries(input_layer, reference_layer, snap, out_)

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(snp_name, loop_n, len(layers)-1))
            print('Importing SNAPPED LAYER to GRASS GIS')
            self.gs.run_command(
                'v.in.ogr', overwrite=True, input=tmp_dbase, layer='tmp_snp2',
                output='tmp_snp'
            )

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(snp_name, loop_n, len(layers)-1))
            print('Converting SNAPPED LAYER into RASTER')
            self.gs.run_command(
                'v.to.rast', overwrite=True, input='tmp_snp', output='tmp_snp',
                use='attr', attribute_column='priority'
            )

            if len(layers) > 2:
                ref_name = 'tmp_ref' if i < len(layers)-1 else output
            else:
                ref_name = output
                
            new_ref_label = '-'.join(layers[0:loop_n + 1])

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Patching REFERENCE and SNAPPED LAYERS')
            self.gs.run_command('r.patch', input='tmp_ref,tmp_snp', output='tmp_patch', overwrite=True)

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Converting PATCHED RASTER into VECTOR')
            self.gs.run_command(
                'r.to.vect', flags='s', overwrite=True, input='tmp_patch',
                output='tmp_patch', type='area', column='priority'
            )

            # Remove features with area greater then threshold
            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Removing areas smaller than \t<< [{}] >> units'.format(area_thr))
            self.gs.run_command(
                'v.clean', flags='c', input='tmp_patch', output='tmp_rmarea',
                tool='rmarea', threshold=area_thr, overwrite=True
            )

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Cleaning PATCHED VECTOR geometries')
            self.gs.run_command(
                'v.clean', flags='c', input='tmp_rmarea', output='tmp_cleaned',
                tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
                threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True
            )

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Simplifying PATCHED VECTOR geometries')
            self.gs.run_command(
                'v.generalize', overwrite=True, input='tmp_cleaned',
                output='tmp_simplified', method='douglas', threshold=simplify
            )

            # Export-Import cleaning
            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Cleaning (Export-Import) PATCHED VECTOR')
            flags = 'su' if os.path.exists(tmp_dbase) else 's'
            self.gs.run_command(
                'v.out.ogr', flags=flags, input='tmp_simplified', output=tmp_dbase,
                output_layer='tmp_simplified', format='GPKG', overwrite=True
            )
            self.gs.run_command(
                'v.in.ogr', overwrite=True, input=tmp_dbase,
                layer='tmp_simplified', output='tmp_out_in'
            )

            print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
            print('Simplifying PATCHED VECTOR geometries')
            self.gs.run_command(
                'v.generalize', overwrite=True, input='tmp_out_in',
                output=ref_name, method='douglas', threshold=simplify
            )

            # Next loop
            i += 2 if i == 0 else +1
            loop_n += 1
            if i < len(layers):

                print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
                print('Converting NEW REFERENCE into RASTER')
                self.gs.run_command(
                    'v.to.rast', overwrite=True, input=ref_name, output=ref_name,
                    use='attr', attribute_column='priority'
                )

                print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
                print('Exporting NEW REFERENCE to auxiliary database')
                flags = 'su' if os.path.exists(tmp_dbase) else 's'
                self.gs.run_command(
                    'v.out.ogr', flags=flags, input=ref_name, output=tmp_dbase,
                    output_layer=ref_name, format='GPKG', overwrite=True
                )

                print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(new_ref_label, loop_n, len(layers)-1))
                print('Reading NEW REFERENCE with Geopandas')
                ref = gpd.read_file(tmp_dbase, layer=ref_name, driver='GPKG')
                ref['geometry'] = ref['geometry'].map(lambda x: pygis.segmentize(x, seglen))
                ref.to_file(tmp_dbase, layer=ref_name, driver='GPKG')

                snp_name = layers[i]

                print('\n\n\t\t<<[ {} ]>>\t[loop {}/{}]'.format(snp_name, loop_n, len(layers)-1))
                print('Extracting NEW SNAP LAYER')
                snp = patch.loc[patch['priority'] == priorities[i], :].reset_index(drop=True).copy()
                snp['geometry'] = snp['geometry'].map(lambda x: pygis.segmentize(x, seglen))
                snp.to_file(tmp_dbase, layer=snp_name, driver='GPKG')

        print('\n\n\t\tLAND COVER LAYER CREATED')

        print('\n\nManage LAND COVER attribute table')
        self.gs.run_command('v.db.dropcolumn', map=output, columns='label')
        self.gs.run_command('v.db.addcolumn', map=output, columns='\"class TEXT\"')
        query_file = os.path.join(self.tmpdir, 'table_queries.sql')
        with open(query_file, 'w') as f:
            f.write('BEGIN TRANSACTION;\n')
            for n, cl in zip(priorities, layers):
                f.write('UPDATE {} SET class=\"{}\" where priority=\"{}\";\n'.format(output, cl, n))
            f.write('COMMIT;\n')
        self.gs.run_command('db.execute', input=query_file)

        print('\n\nRemoving temporal files')
        self.gs.run_command('g.remove', flags='f', type='vector,raster', pattern='tmp_*')
        if os.path.exists(tmp_dbase):
            os.remove(tmp_dbase)
        if os.path.exists(query_file):
            os.remove(query_file)
        if os.path.exists(cat_file):
            os.remove(cat_file)


    def consolidate_topology(
        self, boundary, land_cover,
        seglen, snap, simplify, area_thr,
        output_boundary, output_land_cover
    ):

        bnd = boundary
        lc = land_cover
        tmp_dbase = os.path.join(self.tmpdir, 'consolidate.gpkg')

        print('\n\nDISSOLVE land cover')
        self.gs.run_command('g.copy', vector='{},tmp_lc'.format(lc), overwrite=True)
        self.gs.run_command('v.db.addcolumn', map='tmp_lc', columns='\"diss INT\"')
        self.gs.run_command('v.db.update', map='tmp_lc', column='diss', value=1)
        self.gs.run_command('v.dissolve', input='tmp_lc', column='diss', output='tmp_diss', overwrite=True)

        print('\n\nEXPORT dissolved land cover to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command(
            'v.out.ogr', flags=flags, input='tmp_diss', output=tmp_dbase,
            output_layer='tmp_diss', format='GPKG', overwrite=True
        )

        print('\n\nDENSIFYING boundary')
        self.gs.run_command(
            'v.split', flags='n', input=bnd, output='tmp_bnd_split', length=seglen, overwrite=True
        )

        print('\n\nEXPORT boundary to auxiliary database')
        flags = 'su' if os.path.exists(tmp_dbase) else 's'
        self.gs.run_command(
            'v.out.ogr', flags=flags, input='tmp_bnd_split', output=tmp_dbase,
            output_layer='tmp_bnd_split', format='GPKG', overwrite=True
        )

        print('\n\nSNAP boundary to land cover')
        input_layer = '{}|layername=tmp_bnd_split'.format(tmp_dbase)
        reference_layer = '{}|layername=tmp_diss'.format(tmp_dbase)
        out_ = 'ogr:dbname=\'{}\' table=\"tmp_snap\" (geom) sql='.format(tmp_dbase)
        qgs.snap_geometries(input_layer, reference_layer, snap, out_)

        print('\n\nIMPORT snapped boundary into GRASS GIS')
        self.gs.import_interactive(input_path=tmp_dbase, layer='tmp_snap', snap=-1, output='tmp_snap')

        print('\n\nSIMPLIFY boundary geometries')
        self.gs.run_command(
            'v.generalize', overwrite=True, input='tmp_snap',
            output=output_boundary, method='douglas', threshold=simplify
        )

        print('\n\nCLIP land cover with boundary')
        self.gs.clip_layer(layer_a=lc, layer_b=output_boundary, operator='and', output='tmp_clip')

        print('\n\nCLEAN land cover geometries')
        self.gs.run_command(
            'v.clean', flags='c', input='tmp_clip', output=output_land_cover,
            tool='break,rmdupl,bpol,rmdupl,rmdangle,rmdac,rmarea',
            threshold='0,0,0,0,-1,0,{}'.format(area_thr), overwrite=True
        )

        print('\n\nRemoving temporal files')
        self.gs.run_command('g.remove', flags='f', type='vector,raster', pattern='tmp_*')
        if os.path.exists(tmp_dbase):
            os.remove(tmp_dbase)

        print('\n\nEXPORT consolidated Boundary and Land cover layers')
        self.gs.to_gpkg(layer=output_boundary, db=self.output_gpkg, out_layer=output_boundary)
        self.gs.to_gpkg(layer=output_land_cover, db=self.output_gpkg, out_layer=output_land_cover)
        print('\n\nLayers were exported to: {}'.format(self.output_gpkg))


    def purge(self):
        if os.path.exists(self.tmpdir):
            shutil.rmtree(self.tmpdir)
