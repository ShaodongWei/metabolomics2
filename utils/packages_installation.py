
# List of required packages
required_modules = [
    "pandas==2.1.2", "numpy==1.26.1", "openpyxl==3.1.1", "PyYAML==6.0.1", 
    "python-dateutil==2.8.2", "matplotlib==3.8.0", "scikit-learn==1.3.2", 
    "seaborn==0.13.0", "statsmodels==0.14.0"
]

import subprocess
import importlib

for mudule in required_modules:
    subprocess.check_call(['pip', 'install', mudule])

print('\ninstallation finished')