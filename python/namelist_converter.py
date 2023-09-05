import json
import sys
import importlib.util
import subprocess

def install_module_if_not_present(module_name: str) -> None:
    """
    Installs the specified Python module if it's not already present.

    Args:
        module_name (str): The name of the module to check and install.

    Raises:
        subprocess.CalledProcessError: If there is an error installing the module.
    """
    if importlib.util.find_spec(module_name) is None:
        subprocess.check_call(["pip", "install", module_name])
        print(f"Installed {module_name}!")


# Ensure 'f90nml' is installed
install_module_if_not_present('f90nml')

import f90nml

# Get the input filename from command line arguments
input_str = sys.argv[1]

if 'INPUT' in input_str:
    nml = f90nml.read(input_str)
    print(json.dumps(nml, indent=4))
else:
    raise ValueError('Missing INPUT file path!')
