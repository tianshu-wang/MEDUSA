# Functions for interfacing with Fornax during testing

# Modules
from glob import glob
import logging
import os
import subprocess
from timeit import default_timer as timer
from .log_pipe import LogPipe

# Global variables
fornax_rel_path = '../../'
fornax_path = os.path.abspath(fornax_rel_path)
saved_filenames = ['machine.mk', 'setup.mk', 'testing.mk']
saved_files = []
global_config_args = []
global_run_args = []
global_test_name = None
global_coverage_cmd = None
global_silent = False


def find_machine_files(machine=None, include=None, base_dir=None):
    if base_dir is None:
        base_dir = fornax_path
    path = os.path.join(base_dir, 'machines')
    if machine is None:
        machine = subprocess.run(['hostname'], capture_output=True, text=True)
        machine = machine.stdout.strip().split('.')[0]
    mk = [i for i in glob(path + '/*.mk') if os.path.basename(i).startswith(machine)]
    env = [i for i in glob(path + '/*.env') if os.path.basename(i).startswith(machine)]
    if include is not None:
        mk = [i for i in mk if include.lower() in i.lower()]
        env = [i for i in env if include.lower() in i.lower()]
    if len(mk) == 0:
        mk = [path + '/default.mk']
    return mk, env


def add_testing_mk(warnings_to_errors=True, debug=False, base_dir=None, **kwargs):
    def _add(key, value):
        if key in kwargs:
            kwargs[key] = str(kwargs[key]) + ' ' + value
        else:
            kwargs[key] = value

    if base_dir is None:
        base_dir = fornax_path
    if warnings_to_errors:
        _add('CFLAGS', '-Werror')
    if debug:
        _add('CFLAGS', '-g')

    with open(os.path.join(base_dir, 'testing.mk'), 'w') as f:
        for key, value in kwargs.items():
            f.write('{}+={} '.format(key, value))
        f.write('\n')
    return


# Function for configuring Fornax
def configure(setup='ccsn_1d.mk', machine=None, env=False, **kwargs):
    current_dir = os.getcwd()
    os.chdir(fornax_rel_path)
    try:
        if machine is None or env is None:
            if machine is None and env is None:
                machine, env = find_machine_files(include='gnu')
            else:
                files = find_machine_files(machine)
                machine = files[0] if machine is None else machine
                env = files[1] if env is None else env
            machine = machine[0]
            try:
                env = env[0]
            except (IndexError, TypeError):
                env = None
        machine = os.path.abspath(machine)
        if not os.path.isfile(setup):
            setup = os.path.join('setups', setup)
        setup = os.path.abspath(setup)
        logger = logging.getLogger('fornax.make')
        logger.debug('Using {} as machine.mk'.format(machine))
        if os.path.exists('machine.mk'):
            os.remove('machine.mk')
        try:
            os.symlink(machine, 'machine.mk')
        except FileExistsError:
            pass
        if '/' not in setup and not os.path.isfile(setup):
            setup = 'setups/' + setup
        logger.debug('Using {} as setup.mk'.format(setup))
        if os.path.exists('setup.mk'):
            os.remove('setup.mk')
        try:
            os.symlink(setup, 'setup.mk')
        except FileExistsError:
            pass
        if env is not None and env is not False:
            logger.debug('Using {} as env.mk'.format(env))
            out_log = LogPipe('fornax.configure', logging.INFO)
            err_log = LogPipe('fornax.configure', logging.ERROR)
            try:
                subprocess.check_call(['source', env], stdout=out_log, stderr=err_log, shell=True)
            except subprocess.CalledProcessError as err:
                raise FornaxError('Return code {0} from command \'{1}\''
                                  .format(err.returncode, ' '.join(err.cmd)))
            finally:
                out_log.close()
                err_log.close()
    finally:
        os.chdir(current_dir)


# Function for compiling Fornax
def make(clean_first=True, obj_only=False, **kwargs):
    current_dir = os.getcwd()
    os.chdir(fornax_rel_path)
    logger = logging.getLogger('fornax.make')
    out_log = open(os.devnull, 'w') if global_silent else LogPipe('fornax.make',
                                                                  logging.INFO)
    try:
        clean_command = ['make', 'clean']
        if obj_only:
            # used in setup_compile.py to save expensive linking time for Intel Compiler:
            make_command = ['make', '-j8', 'objs']
        else:
            # disable parallel GNU Make execution for Lcov (issues w/ Jenkins filesystem)
            if global_coverage_cmd is not None:
                make_command = ['make']
            else:
                make_command = ['make', '-j8']
        try:
            if clean_first:
                logger.debug('Executing: ' + ' '.join(clean_command))
                subprocess.check_call(clean_command, stdout=out_log)
            logger.debug('Executing: ' + ' '.join(make_command))
            t0 = timer()
            subprocess.check_call(make_command, stdout=out_log)
            logger.debug('Compilation took {0:.3g} seconds.'.format(timer() - t0))
        except subprocess.CalledProcessError as err:
            logger.error("Something bad happened", exc_info=True)
            raise FornaxError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        out_log.close()
        os.chdir(current_dir)


# Functions for running Fornax
def run(input_filename, arguments, lcov_test_suffix=None):
    current_dir = os.getcwd()
    os.chdir('bin')
    out_log = LogPipe('fornax.run', logging.INFO)
    try:
        input_filename_full = '../' + fornax_rel_path + 'inputs/' + \
                              input_filename
        run_command = ['./fornax', input_filename_full]
        try:
            cmd = run_command + arguments + global_run_args
            logging.getLogger('fornax.run').debug('Executing: ' + ' '.join(cmd))
            subprocess.check_call(cmd, stdout=out_log)
        except subprocess.CalledProcessError as err:
            raise FornaxError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        else:
            os.chdir(current_dir)
            # (optional) if execution completes without error, and a lcov_test_suffix is
            # explicitly passed, process Lcov tracefile immediately after run_command
            analyze_code_coverage(global_test_name, lcov_test_suffix)
    finally:
        out_log.close()
        os.chdir(current_dir)


def restart(input_filename, arguments):
    current_dir = os.getcwd()
    os.chdir('bin')
    out_log = LogPipe('fornax.make', logging.INFO)
    try:
        run_command = ['./fornax', '-r', input_filename]
        try:
            cmd = run_command + arguments
            logging.getLogger('fornax.run').debug('Executing (restart): ' + ' '.join(cmd))
            subprocess.check_call(cmd, stdout=out_log)
        except subprocess.CalledProcessError as err:
            raise FornaxError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        out_log.close()
        os.chdir(current_dir)


def mpirun(mpirun_cmd, mpirun_opts, nproc, input_filename, arguments,
           lcov_test_suffix=None):
    current_dir = os.getcwd()
    os.chdir('bin')
    out_log = LogPipe('fornax.run', logging.INFO)
    try:
        input_filename_full = '../' + fornax_rel_path + 'inputs/' + \
                              input_filename
        run_command = [mpirun_cmd] + mpirun_opts + ['-n', str(nproc), './fornax', '-i',
                                                    input_filename_full]
        run_command = list(filter(None, run_command))  # remove any empty strings
        try:
            cmd = run_command + arguments + global_run_args
            logging.getLogger('fornax.run').debug('Executing (mpirun): ' + ' '.join(cmd))
            subprocess.check_call(cmd, stdout=out_log)
        except subprocess.CalledProcessError as err:
            raise FornaxError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        else:
            os.chdir(current_dir)
            # (optional) if execution completes without error, and a lcov_test_suffix is
            # explicitly passed, process Lcov tracefile immediately after run_command
            analyze_code_coverage(global_test_name, lcov_test_suffix)
    finally:
        out_log.close()
        os.chdir(current_dir)


# Function for saving configure-generated files that may already exist
def save_files():
    global saved_files
    for filename in saved_filenames:
        rel_path_to_file = fornax_rel_path + filename
        if os.path.isfile(rel_path_to_file):
            with open(rel_path_to_file, 'r') as current_file:
                saved_files.append(current_file.read())
        else:
            saved_files.append(None)


# Function for restoring configure-generated files that previously existed
def restore_files():
    for filename, saved_file in zip(saved_filenames, saved_files):
        rel_path_to_file = fornax_rel_path + filename
        if saved_file is None:
            os.system('rm -f ' + rel_path_to_file)
        else:
            with open(rel_path_to_file, 'w') as current_file:
                current_file.write(saved_file)


# Function for analyzing code coverage after fornax.run() or fornax.mpirun()
def analyze_code_coverage(test_name, lcov_test_suffix=None):
    # Only run Lcov if a string is passed to optional lcov_test_suffix argument (most
    # regression tests only need to run Lcov ONCE after test.analyze() is complete).

    # Regression tests with multiple fornax.make() calls and binaries require that Lcov is
    # executed in test.run() after intermediate fornax.mpi/run() calls, before the test is
    # complete. Test author must ensure that the obj/ directory contains the
    # correct/matching files when fornax.run() is called with lcov_test_suffix=string
    if lcov_test_suffix is not None and global_coverage_cmd is not None:
        # Empty string --> use base test name for output file
        if lcov_test_suffix == '':
            lcov_test_name = global_test_name
        # Append nonempty string suffix to base test name with an underscore
        else:
            lcov_test_name = '_'.join([global_test_name, lcov_test_suffix])
        # For now, assumes Lcov flags for adding test-dependent info (name, output file):
        test_lcov_cmd = (
            global_coverage_cmd
            + ' --test-name {0} -output-file {0}.info'.format(lcov_test_name)
        )
        # is this usage of os.system() safe?
        os.system(test_lcov_cmd)


# General exception class for these functions
class FornaxError(RuntimeError):
    pass
