import json
import sys
import f90nml

import importlib.util
import subprocess

def install_module_if_not_present(module_name):
    if importlib.util.find_spec(module_name) is None:
        subprocess.check_call(["pip", "install", module_name])
        print(f"Installed {module_name}!")


install_module_if_not_present('f90nml')
input_str = sys.argv[1]
if 'INPUT' in input_str:
    nml = f90nml.read(input_str)
    print(json.dumps(nml, indent=4))
else:
    ValueError('Missing INPUT file path!')