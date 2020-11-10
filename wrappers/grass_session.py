# Modules
import os
import sys
import shutil
import tempfile
import subprocess
import pandas as pd
import geopandas as gpd


# Functions
def make_command(prog, flags=None, overwrite=None, verbose=None, **options):
    args = [prog]

    if overwrite:
        args.append("--overwrite")
    if verbose:
        args.append("--verbose")
    if flags:
        args.append("-{}".format(flags))

    for key, val in options.items():
        args.append('{}={}'.format(key, val))

    return ' '.join(args)


# Classes
class GrassSession(object):
    """
    EXPLANATION
    """

    def __init__(self, gisdb, location='LOCATION', mapset='PERMANENT', crs=None):

        # Initialize variables
        self.grassbin = os.environ['GRASSBIN']
        
        self.gisdb = gisdb
        self.location = location
        self.mapset = mapset
        self.crs = crs

        # Set GISDBASE environment variable
        os.environ['GISDBASE'] = self.gisdb

        # Create location and mapset
        self.create_location()
        if self.mapset != 'PERMANENT':
            self.create_mapset()
        
        print(self.read_command('g.gisenv', verbose=True))

    
    def create_location(self):

        if not os.path.exists(self.gisdb):
            os.makedirs(self.gisdb)

        location_path = os.path.join(self.gisdb, self.location)
        if os.path.exists(location_path):
            shutil.rmtree(location_path)

        # Create grass location
        if self.crs:
            cmd = '\"{}\" -c \"{}\" -e \"{}\"'.format(self.grassbin, self.crs, location_path)
        else:
            cmd = '\"{}\" -c -e \"{}\"'.format(self.grassbin, location_path)

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate(input=os.linesep.encode())

        if 'Press any key to continue' in out.decode() or 'ERROR' in err.decode().upper():
            print(err.decode())
            raise RuntimeError('GRASS GIS location: {} not created'.format(location_path))


    def create_mapset(self):

        mapset_path = os.path.join(self.gisdb, self.location, self.mapset)
        if os.path.exists(mapset_path):
            shutil.rmtree(mapset_path)

        grass_command = 'g.mapset -c mapset=\"{}\" location=\"{}\" dbase=\"{}\"'.format(self.mapset, self.location, self.gisdb)
        mapset_permanent = os.path.join(self.gisdb, self.location, 'PERMANENT')
        cmd = '\"{}\" \"{}\" --exec \"{}\"'.format(self.grassbin, mapset_permanent, grass_command)

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate(input=os.linesep.encode())

        if 'Press any key to continue' in out.decode() or 'ERROR' in err.decode().upper():
            print(err.decode())
            raise RuntimeError('GRASS GIS mapset: {} not created'.format(mapset_path))


    def run_command(self, prog, flags=None, overwrite=None, verbose=None, **options):

        mapset = os.path.join(self.gisdb, self.location, self.mapset)
        cmd = make_command(prog, flags, overwrite, verbose, **options)
        command = '\"{}\" \"{}\" --exec {}'.format(self.grassbin, mapset, cmd)

        p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate(input=os.linesep.encode())

        if 'Press any key to continue' in out.decode() or 'ERROR' in err.decode().upper():
            print(err.decode())
            raise RuntimeError('The algorithm failed to execute. Check ERROR in GRASS output above')
        else:
            print(err.decode())


    def read_command(self, prog, flags=None, overwrite=None, verbose=None, **options):

        location = os.path.join(self.gisdb, self.location, self.mapset)
        cmd = make_command(prog, flags, overwrite, verbose, **options)
        command = '\"{}\" \"{}\" --exec {}'.format(self.grassbin, location, cmd)

        p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = p.communicate(input=os.linesep.encode())
        
        if 'Press any key to continue' in out.decode() or 'ERROR' in err.decode().upper():
            print(err.decode())
            raise RuntimeError('The algorithm failed to execute. Check ERROR in GRASS output above')
        else:
            print(err.decode())
            return out.decode()


    def import_interactive(self, input_path, layer, snap, output):

        location = os.path.join(self.gisdb, self.location, self.mapset)
        input_txt = 'To accept the imported layer press << a >>\n To import using a different snap value, insert the new value << value >>\n'
        x = ''

        while True:
            cmd = 'v.in.ogr --overwrite input={} layer={} snap={} output={}'.format(input_path, layer, snap, output)
            command = '\"{}\" \"{}\" --exec {}'.format(self.grassbin, location, cmd)
            p = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            out, err = p.communicate(input=os.linesep.encode())
    
            if 'Press any key to continue' in out.decode() or 'ERROR' in err.decode().upper():
                print(err.decode())
                raise RuntimeError('The algorithm failed to execute. Check ERROR in GRASS output above')
            else:
                print(err.decode())
                while not x.lower() == 'a' or isinstance(x, float):
                    x = input(input_txt)

                if x.lower() == 'a':
                    break
                else:
                    snap = x


    def clip_layer(self, layer_a, layer_b, operator, output):
        # Clip
        self.run_command('v.overlay', overwrite=True, verbose=True, ainput=layer_a, binput=layer_b, operator=operator, output=output)

        # Manage attribute table
        attrs_table = self.read_command('v.db.select', map=output)
        colnames = attrs_table.split(os.linesep)[0].split('|')
        drop_cols = [i for i in colnames if i.startswith('b_')]
        drop_cols.append(colnames[1])
        self.run_command('v.db.dropcolumn', map=output, columns=','.join(drop_cols))

        old_cols = [i for i in colnames if i not in drop_cols and i != 'cat']
        new_cols = [i[2::] for i in old_cols]

        for old, new in zip(old_cols, new_cols):
            self.run_command('v.db.renamecolumn', map=output, column='{},{}'.format(old, new))


    def merge_layers(self, layers, output):

        count_loops = 1
        c = 0
        la = layers[c]
        lb = layers[c + 1]

        # loops
        while c < len(layers):
            # Check columns names
            colsdata_la = self.read_command('db.describe', flags='c', table=la)
            colsdata_la = colsdata_la.strip().split(os.linesep)[2::]
            colsdata_la = ['{} {}'.format(i.split(':')[1].strip(), i.split(':')[2].strip()) for i in colsdata_la][1::]
            la_names = [i.split(' ')[0] for i in colsdata_la]

            colsdata_lb = self.read_command('db.describe', flags='c', table=lb)
            colsdata_lb = colsdata_lb.strip().split(os.linesep)[2::]
            colsdata_lb = ['{} {}'.format(i.split(':')[1].strip(), i.split(':')[2].strip()) for i in colsdata_lb][1::]
            lb_names = [i.split(' ')[0] for i in colsdata_lb]

            # Overlay
            if len(layers) > 2:
                patch = 'tmp_patch_{}'.format(count_loops) if c < len(layers) - 1 else output
            else:
                patch = output

            self.run_command('v.overlay', ainput=la, binput=lb, operator='or', output=patch, overwrite=True,
                             verbose=True)

            # Manage attribute table
            if sorted(la_names) == sorted(lb_names):
                tmp_file = os.path.join(tempfile.mkdtemp(prefix='tmp_query_merge_layers_'), 'query.sql')
                with open(tmp_file, 'w') as f:
                    f.write('BEGIN TRANSACTION;\n')
                    for col in la_names:
                        a_col = 'a_{}'.format(col)
                        b_col = 'b_{}'.format(col)
                        f.write('UPDATE {} SET {}={} WHERE {} IS NULL;\n'.format(patch, a_col, b_col, a_col))
                    f.write('COMMIT;\n')
                self.run_command('db.execute', input=tmp_file)

                for col in la_names:
                    a_col = 'a_{}'.format(col)
                    self.run_command('v.db.renamecolumn', map=patch, column='\"{},{}\"'.format(a_col, col))

                # Drop columns
                to_drop = self.read_command('db.describe', flags='c', table=patch).strip().split(os.linesep)[2::]
                to_drop = [i.split(':')[1].strip() for i in to_drop]
                to_drop = ','.join([i for i in to_drop if i.startswith('a_') or i.startswith('b_')])
                self.run_command('v.db.dropcolumn', map=patch, columns=to_drop)

            else:
                ov_cols = self.read_command('db.describe', flags='c', table=patch).strip().split(os.linesep)[2::]
                ov_cols = [i.split(':')[1].strip() for i in ov_cols]
                a_cat = [i for i in ov_cols if i.startswith('a_')][0]
                b_cat = [i for i in ov_cols if i.startswith('b_')][0]
                to_drop = ','.join([a_cat, b_cat])
                self.run_command('v.db.dropcolumn', map=patch, columns=to_drop)

                for col in la_names:
                    a_col = 'a_{}'.format(col)
                    self.run_command('v.db.renamecolumn', map=patch, column='\"{},{}\"'.format(a_col, col))

                for col in lb_names:
                    b_col = 'b_{}'.format(col)
                    self.run_command('v.db.renamecolumn', map=patch, column='\"{},{}\"'.format(b_col, col))

            # Counter increment
            c += 2 if c == 0 else +1
            count_loops += 1
            if c < len(layers):
                la = patch
                lb = layers[c]

        # Remove temporary
        self.run_command('g.remove', type='vector', pattern='tmp_patch*', flags='f')
        if os.path.exists(tmp_file):
            shutil.rmtree(os.path.dirname(tmp_file))


    def to_gpd(self, layer):
        tmp_file = os.path.join(tempfile.mkdtemp(prefix='tmp_gpd_'), layer + '.geojson')
        print('\n\nExporting {} to GeodataFrame'.format(layer))
        flags = 'su' if os.path.exists(tmp_file) else 's'
        self.run_command('v.out.ogr', flags=flags, input=layer, output=tmp_file,
                         format='GeoJSON', overwrite=True)

        out = gpd.read_file(tmp_file, layer=layer, driver='GeoJSON')

        if os.path.exists(tmp_file):
            shutil.rmtree(os.path.dirname(tmp_file))

        return out


    def to_gpkg(self, layer, db, out_layer):

        flags = 'su' if os.path.exists(db) else 's'
        self.run_command('v.out.ogr', flags=flags, input=layer, output=db,
                         output_layer=out_layer, format='GPKG', overwrite=True)


    def get_table(self, layer):
        colsdata = self.read_command('db.describe', flags='c', table=layer)
        colsdata = colsdata.strip().split(os.linesep)[2::]
        colsdata = ['{} {}'.format(i.split(':')[1].strip(), i.split(':')[2].strip()) for i in colsdata]
        cols_names = [i.split(' ')[0] for i in colsdata]

        data0 = self.read_command('v.db.select', map=layer)
        data = [i.split('|') for i in data0.strip().split(os.linesep)[1::]]

        df = pd.DataFrame.from_records(data, columns=cols_names).set_index(cols_names[0])
        df.index = df.index.astype(int)

        return df

    def update_table(self, layer):

        # Get DataFrame with attribute table
        colsdata = self.read_command('db.describe', flags='c', table=layer)
        colsdata = colsdata.strip().split(os.linesep)[2::]
        colsdata = ['{} {}'.format(i.split(':')[1].strip(), i.split(':')[2].strip()) for i in colsdata]
        cols_names = [i.split(' ')[0] for i in colsdata]

        data0 = self.read_command('v.db.select', map=layer)
        data = [i.split('|') for i in data0.strip().split(os.linesep)[1::]]

        df = pd.DataFrame.from_records(data, columns=cols_names).set_index(cols_names[0])
        df.index = df.index.astype(int)

        # Get cat values
        cat_val = self.read_command('v.category', input=layer, option='print').strip().split(os.linesep)
        cat_val = [int(i.split('/')[0]) for i in cat_val]

        self.run_command('v.category', overwrite=True, input=layer, output='tmp_catdel', option='del', cat='-1')
        self.run_command('v.category', overwrite=True, input='tmp_catdel', output=layer, option='add')
        self.run_command('v.db.droptable', flags='f', map=layer)
        self.run_command('v.db.addtable', map=layer)

        if colsdata != '':
            self.run_command('v.db.addcolumn', map=layer, columns='\"' + ','.join(colsdata[1::]) + '\"')
            new_df = df.loc[cat_val, :]
            query_values = new_df.columns.values + '=' + '\"' + new_df.iloc[:, :].astype(str) + '\"'
            query_values = query_values.apply(','.join, axis=1)
            query_values = query_values.str.replace('\"\"', 'null', regex=False)

            tmp_file = os.path.join(tempfile.mkdtemp(prefix='tmp_update_table_'), 'tmp_query.sql')

            with open(tmp_file, 'w') as f:
                f.write('BEGIN TRANSACTION;\n')
                for i, query in enumerate(query_values):
                    f.write('UPDATE {} SET {} WHERE cat={};\n'.format(layer,
                                                                      query,
                                                                      i + 1))
                f.write('COMMIT;\n')

            self.run_command('db.execute', input=tmp_file)

        self.run_command('g.remove', flags='f', type='vector', name='tmp_catdel')

        if os.path.exists(tmp_file):
            shutil.rmtree(os.path.dirname(tmp_file))
