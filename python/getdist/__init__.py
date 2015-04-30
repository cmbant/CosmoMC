__author__ = 'Antony Lewis'
__all__ = ['ResultObjs', 'paramNames', 'MCSamples', 'chains','inifile']
__version__ = "0.2.0"

from getdist.inifile import IniFile
from getdist.MCSamples import MCSamples, loadMCSamples
from getdist.chains import WeightedSamples

output_base_dir = None
cache_dir = None
use_plot_data = False

def get_defaults():
    import os
    return os.path.join(os.path.dirname(__file__), 'analysis_defaults.ini')

def set_logging(loglevel):
    import logging
    logging.basicConfig(level=loglevel)

def get_config():
    import os
    config_file = os.environ.get('GETDIST_CONFIG', None)
    if not config_file:
        config_file = os.path.join(os.path.dirname(__file__), 'config.ini')
    if os.path.exists(config_file):
        return IniFile(config_file)
    else:
        return IniFile()

config_ini = get_config()
default_grid_root = config_ini.string('default_grid_root', '')
output_base_dir = config_ini.string('output_base_dir', '')
cache_dir = config_ini.string('cache_dir', '')
default_getdist_settings = config_ini.string('default_getdist_settings', get_defaults())
use_plot_data = config_ini.bool('use_plot_data', use_plot_data)
default_plot_output = config_ini.string('default_plot_output', 'pdf')
loglevel = config_ini.string('logging', '')
if loglevel: set_logging(loglevel)
