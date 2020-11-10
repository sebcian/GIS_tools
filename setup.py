# Modules
import os
import sys
import yaml

# Read configuration file
with open('configuration_file.yml') as f:
    config_ = yaml.load(f, Loader=yaml.FullLoader)

# Path to conda environment folders
conda_env_activate = os.path.join(os.environ['CONDA_PREFIX'], 'etc', 'conda', 'activate.d')
conda_env_deactivate = os.path.join(os.environ['CONDA_PREFIX'], 'etc', 'conda', 'deactivate.d')

# Create dirs
if not os.path.exists(conda_env_activate):
    os.makedirs(conda_env_activate)

if not os.path.exists(conda_env_deactivate):
    os.makedirs(conda_env_deactivate)

# Create GIS tools environment files
grassbin = config_['grass_executable']

sagabin = config_['saga_cmd_executable']
sagatools = config_['saga_modules_tools_folder']

qgisbin = config_['qgis_python_executable']
qgisprefixpath = config_['qgis_prefix_path']
qgisplugin = config_['qgis_python_plugins_folder']


if 'PYTHONPATH' in os.environ.keys():
    pp = os.environ['PYTHONPATH']
    pp = pp.split(';')
    pp = list(set(pp))
    pp = ';'.join(pp)

else:
    pp = os.getcwd()

pythonpath = pp

if 'win' in sys.platform:
    # Bat file (Windows)
    
    with open(os.path.join(conda_env_activate, 'GIS_tools.bat'), 'w') as f:
        
        f.writelines('@echo off\n')
        
        # GRASS GIS
        f.writelines('set GRASSBIN={}\n'.format(grassbin))
        
        # SAGA GIS
        f.writelines('set SAGABIN={}\n'.format(sagabin))
        f.writelines('set SAGA_MLB={}\n'.format(sagatools))
        
        # QGIS
        f.writelines('set QGISBIN={}\n'.format(qgisbin))
        f.writelines('set QGISPREFIX_PATH={}\n'.format(qgisprefixpath))
        f.writelines('set QGISPLUGIN={}\n'.format(qgisplugin))
        
        # PYTHONPATH
        f.writelines('set PYTHONPATH={}\n'.format(pythonpath))


    with open(os.path.join(conda_env_deactivate, 'GIS_tools.bat'), 'w') as f:

        f.writelines('@echo off\n')
        
        # GRASS GIS
        f.writelines('set GRASSBIN=\n')
        
        # SAGA GIS
        f.writelines('set SAGABIN=\n')
        f.writelines('set SAGA_MLB=\n')
        
        # QGIS
        f.writelines('set QGISBIN=\n')
        f.writelines('set QGISPREFIX_PATH=\n')
        f.writelines('set QGISPLUGIN=\n')
        
        # PYTHONPATH
        if 'PYTHONPATH' in os.environ.keys():
            f.writelines('set PYTHONPATH={}\n'.format(os.environ['PYTHONPATH']))
        else:
            f.writelines('set PYTHONPATH=\n')


if 'linux' in sys.platform:
    # Bash file (Linux)
    
    with open(os.path.join(conda_env_activate, 'GIS_tools.sh'), 'w') as f:
        
        f.writelines('#!/bin/bash\n')
        
        # GRASS GIS
        f.writelines('export GRASSBIN={}\n'.format(grassbin))
        
        # SAGA GIS
        f.writelines('export SAGABIN={}\n'.format(sagabin))
        f.writelines('export SAGA_MLB={}\n'.format(sagatools))
        
        # QGIS
        f.writelines('export QGISBIN={}\n'.format(qgisbin))
        f.writelines('export QGISPREFIX_PATH={}\n'.format(qgisprefixpath))
        f.writelines('export QGISPLUGIN={}\n'.format(qgisplugin))
        
        # PYTHONPATH
        f.writelines('export PYTHONPATH={}\n'.format(pythonpath))


    with open(os.path.join(conda_env_deactivate, 'GIS_tools.sh'), 'w') as f:

        f.writelines('#!/bin/bash\n')
        
        # GRASS GIS
        f.writelines('unset GRASSBIN\n')
        
        # SAGA GIS
        f.writelines('unset SAGABIN\n')
        f.writelines('unset SAGA_MLB\n')
        
        # QGIS
        f.writelines('unset QGISBIN\n')
        f.writelines('unset QGISPREFIX_PATH\n')
        f.writelines('unset QGISPLUGIN\n')
        
        # PYTHONPATH
        if 'PYTHONPATH' in os.environ.keys():
            f.writelines('export PYTHONPATH={}\n'.format(os.environ['PYTHONPATH']))
        else:
            f.writelines('unset PYTHONPATH\n')
