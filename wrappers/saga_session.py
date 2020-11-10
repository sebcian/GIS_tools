# Modules
import os
import subprocess


# Functions
def make_command(library, module, **options):
    args = ['saga_cmd', library, str(module)]

    for key, val in options.items():
        args.append('-{} {}'.format(key.upper(), val))

    return ' '.join(args)


def get_version():
    env = os.environ
    env['PATH'] = ';'.join([
        os.path.dirname(env['SAGABIN']),
        env['PATH']
    ])

    cmd = 'saga_cmd --version'
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, env=env)
    out, err = p.communicate()

    if err:
        print(err.decode())
        raise RuntimeError('The algorithm failed to execute. Check SAGA-GIS output above')
    else:
        print(out.decode())


def run_command(library, module, **options):
    env = os.environ
    env['PATH'] = ';'.join([
        os.path.dirname(env['SAGABIN']),
        env['PATH']
    ])
    
    cmd = make_command(library, module, **options)
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, env=env)
    out, err = p.communicate()

    if err:
        print(err.decode())
        raise RuntimeError('The algorithm failed to execute. Check SAGA-GIS output above')
    else:
        print(out.decode())
