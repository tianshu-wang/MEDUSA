# Test script to make sure all setups in setups compile

# Modules
import logging
import scripts.utils.fornax as fornax
import glob
import os

current_dir = os.getcwd()
test = current_dir[0:(len(current_dir) - 14)]
setup_directory = test + 'setups/'
logger = logging.getLogger('fornax' + __name__[7:])

# set pgen_choices to list of .cpp files in src/pgen/
setup_choices = [os.path.basename(i) for i in glob.glob(setup_directory + '*.mk')]
setup_choices = ['rad_ccsn_{}d.mk'.format(i) for i in '123']
setup_choices = sorted(setup_choices)


# Prepare Fornax
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    for setup in setup_choices:
        logger.debug('Compiling setup ' + setup)
        fornax.add_testing_mk()
        fornax.configure(setup, **kwargs)
        fornax.make(**kwargs)


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
