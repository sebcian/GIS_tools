# Modules
import os
import sys
import shutil
import json
import subprocess
import tempfile


def init_qgis(filepath):
    with open(filepath, 'w') as f:
        f.writelines('import os\n')
        f.writelines('import sys\n')
        f.writelines('sys.path.append(r\'{}\')\n'.format(os.environ['QGISPLUGIN']))
        f.writelines('import qgis\n')
        f.writelines('from qgis.core import *\n')
        f.writelines('import processing\n')
        f.writelines('QgsApplication.setPrefixPath(os.environ[\'QGIS_PREFIX_PATH\'], True)\n')
        f.writelines('app = QgsApplication([], False)\n')
        f.writelines('app.initQgis()\n')
        f.writelines('QgsApplication.processingRegistry().addProvider(qgis.analysis.QgsNativeAlgorithms())\n')
        f.writelines('processing.core.Processing.Processing.initialize()\n')
        f.writelines('context = processing.tools.dataobjects.createContext()\n')
        f.writelines('context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)\n')


def run_script(script):

    cmds = [
        '"{}"'.format(os.environ['QGISBIN']),
        '"{}"'.format(script)
    ]
    cmd = ' '.join(cmds)
    
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    
    if p.returncode != 0 and not 'deprecationwarning' in err.decode().lower():
        print(err.decode())
        raise ValueError('There were errors executing QGIS Processing algorithm. Check algorithm parameters')
    else:
        print(out.decode('windows-1252'))

def qgis_version():
    tmp_dir = tempfile.mkdtemp(prefix='tmp_qgis_')
    tmp_file = os.path.join(tmp_dir, r'version.py')

    init_qgis(tmp_file)
    with open(tmp_file, 'a') as f:
        f.writelines('print(Qgis.QGIS_VERSION)\n')
    
    run_script(tmp_file)

    if os.path.exists(tmp_file):
        shutil.rmtree(tmp_dir)
    
    
def find_processing_tool(keyword='buffer'):
    
    tmp_dir = tempfile.mkdtemp(prefix='tmp_qgis_')
    tmp_file = os.path.join(tmp_dir, r'find_processing_tool.py')

    init_qgis(tmp_file)

    with open(tmp_file, 'a') as f:
        f.writelines('for alg in QgsApplication.processingRegistry().algorithms():\n')
        f.writelines('\tif \'{}\' in alg.id():\n'.format(keyword))
        f.writelines('\t\tprint(alg.id())\n')
    
    run_script(tmp_file)

    if os.path.exists(tmp_file):
        shutil.rmtree(tmp_dir)


def help_processing_tool(tool):

    tmp_dir = tempfile.mkdtemp(prefix='tmp_qgis_')
    tmp_file = os.path.join(tmp_dir, r'help_processing_tool.py')

    init_qgis(tmp_file)

    with open(tmp_file, 'a') as f:
        f.writelines('processing.algorithmHelp(\'{}\')\n'.format(tool))
    
    run_script(tmp_file)

    if os.path.exists(tmp_file):
        shutil.rmtree(tmp_dir)


def write_processing(tool, parameters, filepath):

    with open(filepath, 'a') as f:
        f.writelines('tool = \'{}\'\n'.format(tool))
        f.writelines('parameters = {}\n'.format(json.dumps(parameters)))
        f.writelines('print(processing.run(tool, parameters, context=context))\n')


def run_processing(tool, parameters):

    tmp_dir = tempfile.mkdtemp(prefix='tmp_qgis_')
    tmp_file = os.path.join(tmp_dir, r'run_processing_tool.py')

    init_qgis(tmp_file)
    write_processing(tool, parameters, tmp_file)
    run_script(tmp_file)

    if os.path.exists(tmp_file):
        shutil.rmtree(tmp_dir)
    

def snap_geometries(input_layer, reference_layer, tolerance, output):

    tmp_dir = tempfile.mkdtemp(prefix='tmp_qgis_')
    tmp_script = os.path.join(tmp_dir, r'snap_geometries.py')
    tmp_gpkg = os.path.join(tmp_dir, r'snap.gpkg')

    init_qgis(tmp_script)

    tool = 'qgis:snapgeometries'
    parameters = {
        'INPUT': input_layer,
        'REFERENCE_LAYER': reference_layer,
        'BEHAVIOR': 0,
        'TOLERANCE': tolerance,
        'OUTPUT': 'ogr:dbname=\'{}\' table=\"snapped\" (geom) sql='.format(tmp_gpkg)
    }
    write_processing(tool, parameters, tmp_script)

    tool = 'native:fixgeometries'
    parameters = {
        'INPUT': '{}|layername=snapped'.format(tmp_gpkg),
        'OUTPUT': 'ogr:dbname=\'{}\' table=\"snapped_fixed\" (geom) sql='.format(tmp_gpkg)
    }
    write_processing(tool, parameters, tmp_script)

    tool = 'native:multiparttosingleparts'
    parameters = {
        'INPUT': '{}|layername=snapped_fixed'.format(tmp_gpkg),
        'OUTPUT': output
    }
    write_processing(tool, parameters, tmp_script)

    run_script(tmp_script)
