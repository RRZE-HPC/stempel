import os
import errno
import math
import subprocess
import sys
import logging
import argparse
from shutil import copyfile
from glob import glob
import pathlib
from datetime import datetime

from ruamel import yaml
from kerncraft import kerncraft
from kerncraft.likwid_bench_auto import get_machine_topology

#from . import stempel
stempel_logpath = pathlib.Path(pathlib.Path.home(), '.stempel')
stempel_logpath.mkdir(exist_ok=True, parents=True)

current_time = datetime.now().strftime("%Y%m%d-%H%M")
log_name = 'myanalysis_{}.log'.format(current_time)
stempel_logfilepath = str(pathlib.Path(stempel_logpath, log_name))
logging.basicConfig(level=logging.INFO, filename=stempel_logfilepath,
                    format='%(asctime)s %(levelname)s %(message)s')


def create_parser():
    """This method creates a parser
    """
    example_gen = 'Example usage: analysis -w ~/workspace -i'
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-w', '--workspace', metavar=('WORKSAPCE'),
                        required=True, help='Path to the workspace')
    parser.add_argument('-i', '--iaca', action='store_true', default=False,
                        help='Defines wether to run with IACA or not')
    parser.add_argument('-p', '--prova', nargs=2, metavar=('PROVAPATH', 'PROVAWORKSPACE'),
                        help='Defines wether to run an experiment through PROVA! or not')
    parser.add_argument('-l', '--likwid', nargs=2, metavar=('LIKWID_INC_PATH', 'LIKWID_LIB_PATH'),
                        help='Specifies the path to the likwid include and library respectively.')
    parser.add_argument('-k', '--kind', metavar=('KIND'),
                        choices=['star', 'box'], help='Kind of stencil')
    parser.add_argument('-m', '--machinepath', metavar=('MACHINEPATH'),
                        required=True, help='Path to the machinefile(s)')
    parser.add_argument('-r', '--radius', type=int, metavar=('RADIUS'),
                        help='Radius of the stencil')
    parser.add_argument('-d', '--dimensions', metavar=('DIMENSIONS'),
                        type=int, help='Number of dimensions of the stencil')
    parser.add_argument('-e', '--executions', metavar=('EXECUTIONS'),
                        type=int, default=5, help='Number of executions of the code')
    parser.add_argument('-t', '--threads', metavar=('THREADS'),
                        nargs='*', help='Number of dimensions of the stencil')
    parser.add_argument('--method_type', metavar=('METHODTYPE'),
                        default='OpenMP-4.0-GCC-4.9.3-2.25', help='Method type to use for the actual implementation')
    parser.add_argument('-c', '--classification', metavar=('CLASSIFICATION'),
                        choices=['isotropic', 'heterogeneous', 'homogeneous', 'point-symmetric'], help='Classification of the stencil')
    parser.add_argument('-C', '--coefficients', metavar=('COEFFICIENTS'),
                        choices=['constant', 'variable'], help='Kind of the coefficients of the stencil')
    parser.set_defaults(func=run_gen)

    return parser


def getlast_dir(mypath):
    try:
        directory = max(glob(os.path.join(mypath, '*/')), key=os.path.getmtime)
    except ValueError as e:
        directory = ''
    return directory


def copyheaders(headers_path, destination_dir):
    headers = os.listdir(headers_path)
    for h in headers:
        filesrc = os.path.join(headers_path, h)
        filedest = os.path.join(destination_dir, h)
        try:
            copyfile(filesrc, filedest)
        except IOError as e:
            logging.error('Unable to copy file {} to {}. {}'.format(
                filesrc, filedest, e))


def basesetup(provapath, provaworkspace, likwid_inc, likwid_lib):
    os.environ['BENCHROOT'] = provapath
    os.environ['WORKSPACE'] = provaworkspace
    os.environ['PATH'] += ':' + provapath
    os.environ['LIKWID_LIB'] = likwid_lib
    os.environ['LIKWID_INC'] = likwid_inc
    if 'LD_LIBRARY_PATH' in os.environ:
        os.environ['LD_LIBRARY_PATH'] += ':' + likwid_lib
    else:
        os.environ['LD_LIBRARY_PATH'] = likwid_lib
    logging.info('LIKWID_INC: {}'.format(likwid_inc))

    method_avail = os.path.join(
        provaworkspace, 'methodType', 'methodType_avail')
    method_installed = os.path.join(
        provaworkspace, 'methodType', 'methodType_installed')
    pathlib.Path(method_avail).mkdir(parents=True, exist_ok=True)
    pathlib.Path(method_installed).mkdir(parents=True, exist_ok=True)


def create_project(provapath, provaworkspace, project, params, values, threads):

    # create a project
    mypathname = os.path.join(provaworkspace, project)
    if os.path.exists(mypathname):
        logging.error(
            "Run failed: The project {} already exists. Please choose a different name or delete it first.".format(project))
        print("Run failed: The project {} already exists. Please choose a different name or delete it first.".format(project))
        sys.exit(1)

    cmd = ['workflow', 'project', '-c', '-p', project, '--params',
           params, '--values', values, '--threads', str(threads)]
    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        out = subprocess.check_output(cmd)
        logging.info('Returned: {}'.format(out))
    except subprocess.CalledProcessError as e:
        logging.error("Run failed: {}".format(e))
        sys.exit(1)


def create_method(provapath, provaworkspace, project, method_type, method_name):

    mypathname = os.path.join(provaworkspace, project, method_name)
    if os.path.exists(mypathname):
        logging.error(
            "Run failed: The method {} already exists in this project. Please choose a different name or delete it first.".format(method_name))
        print("Run failed: The method {} already exists in this project. Please choose a different name or delete it first.".format(method_name))
        sys.exit(1)
    #create a method
    cmd = ['workflow', 'method', '-c', '-p', project,
           '-m', method_type, '-n', method_name]
    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        out = subprocess.check_output(cmd)
        logging.info('Returned: {}'.format(out))
    except subprocess.CalledProcessError as e:
        logging.error("Run failed: {}".format(e))
        sys.exit(1)


def run_exp(provapath, provaworkspace, project, executions, param_values, method_name, threads, pinning):
    # cmd_prefix = ['.', os.path.join(
    #     provapath, 'util', 'BaseSetup.sh'), provaworkspace, '&&']

    exp_threads = threads.split()
    # compile the code and run an experiment of the method
    cmd = ['workflow', 'run_exp', '-p', project, '-e', str(executions),
           '-d'] + param_values + ['-m', method_name, '-t'] + exp_threads + ['--pin', pinning]

    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        out = subprocess.check_output(cmd)
        logging.info('Returned: {}'.format(out))
    except subprocess.CalledProcessError as e:
        #print("Run failed:", e)
        logging.error("Run failed: {}".format(e))
        sys.exit(1)


def build_graph(project, exp_dir, param_values, method_name, threads, pinning):
    # cmd_prefix = ['.', os.path.join(
    #     provapath, 'util', 'BaseSetup.sh'), provaworkspace, '&&']
    exp_threads = threads.split()
    # compile the code and run an experiment of the method
    method_name += '_' + pinning
    cmd = ['workflow', 'build_graph', '-p', project, '-e', exp_dir,
           '-d'] + param_values + ['-m', method_name, '-t'] + exp_threads + ['-f', '0', '-T', 'stdev', '-M', 'MLUP/s']
    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        out = subprocess.check_output(cmd)
        logging.info('Returned: {}'.format(out))
    except subprocess.CalledProcessError as e:
        #print("Run failed:", e)
        logging.error("Run failed: {}".format(e))
        sys.exit(1)

def logspace(start, end, num):
    start = math.log10(start)
    end = math.log10(end)
    steplength = float((end-start))/float((num-1))
    i = 0
    while i < num:
        yield int(10**(start+i*steplength))
        i += 1

def run_gen(args, output_file=sys.stdout):

    if args.kind:
        kind = [args.kind]
    else:
        kind = ['star', 'box']

    workspace = args.workspace

    machinefilepath = os.path.dirname(args.machinepath)
    machinefiles = [os.path.basename(args.machinepath)]

    stencilfiles = os.path.join(workspace, 'stencils')
    radius_star = None
    radius_box = None

    if args.radius:
        radius_star = [args.radius]
        radius_box = [args.radius]
    else:
        radius_star = [1, 2, 3, 4]
        radius_box = [1, 2]
    if args.dimensions:
        dimensions = [args.dimensions]
    else:
        dimensions = [2, 3]

    if args.classification:
        classification = [args.classification]
    else:
        classification = ['isotropic', 'heterogeneous',
                          'homogeneous', 'point-symmetric']

    if args.coefficients:
        coefficients = [args.coefficients]
    else:
        coefficients = ['constant', 'variable']

    iaca = args.iaca
    if args.prova and args.likwid:
        provapath, provaworkspace = args.prova
        likwid_inc, likwid_lib = args.likwid
        withprova = True
    else:
        withprova = None
        logging.info(
            'Executing without PROVA because either PROVAPATH/PROVAWORKSPACE or LIKWID_INC/LIKWID_LIB is missing.')

    if args.threads:
        exp_threads = ''
        for t in args.threads:
            exp_threads += t + ' '
        exp_threads = exp_threads.rstrip()
    else:
        exp_threads = '2'

    executions = args.executions

    method_type = args.method_type
    method_name = 'openMP'
    pinning = 'node'

    for k in kind:
        if k == 'star':
            radius = radius_star
            cache_model = "LC"
        else:
            radius = radius_box
            cache_model = "SIM"
        for c in coefficients:
            for l in classification:
                for r in radius:
                    for d in dimensions:
                        stencil_name = str(d) + 'd-' + str(r) + 'r-' + l + \
                            '-' + c + '-' + k + '-stencil.c'
                        logging.info('Working on: {}'.format(stencil_name))
                        stencil_path = os.path.join(stencilfiles, str(
                            d) + 'D', str(r) + 'r', l, k, c)
                        try:
                            os.makedirs(stencil_path)
                        except OSError as e:
                            if e.errno != errno.EEXIST:
                                raise
                        cmd = ['stempel', 'gen', '-D', str(d), '-r', str(
                            r), '-k', k, '-C', c, '--' + l, '--store', os.path.join(stencil_path, stencil_name)]
                        try:
                            # print(cmd)
                            logging.info(
                                'Running command: {}'.format(' '.join(cmd)))
                            subprocess.check_output(cmd)
                        except subprocess.CalledProcessError as e:
                            #print("Run failed:", e)
                            logging.error(
                                'Failed to execute {}: {}'.format(cmd, e))
                            sys.exit(1)

                        if d == 2:
                            logging.info('Retrieving machinefile')
                            for machine in machinefiles:
                                with open(os.path.join(machinefilepath, machine)) as myFile:
                                    results = yaml.load(
                                        myFile, Loader=yaml.RoundTripLoader)
                                try:
                                    dimension = results['memory hierarchy'][
                                        2]['size per group']
                                    dimension = dimension.partition('.')[0]
                                except:
                                    line = results['memory hierarchy'][
                                        2]['cache per group']
                                    sets = line['sets']
                                    ways = line['ways']
                                    cl_size = line['cl_size']
                                    dimension = int(
                                        sets) * int(ways) * int(cl_size) / 1048576
                                    dimension = str(
                                        dimension).partition(' ')[0]
                                size = float(dimension) * 1000000 / 8
                                size = int(math.sqrt(size))
                                half_size = str(int(size / 2))
                                quarter_size = str(int(size / 4))
                                double_size = str(int(size * 2))
                                size = str(int(size))

                                sizes = [quarter_size, half_size, size, double_size]
                                param_values = []
                                for size in sizes:
                                    if iaca:
                                        ECM = 'ECM'
                                    else:
                                        ECM = 'ECMData'
                                    cmd = ['kerncraft', '-p', cache_model, '-p', 'Roofline', '-p', ECM, os.path.join(
                                        stencil_path, stencil_name), '-m', os.path.join(machinefilepath, machine), '-D', 'M', size, '-D', 'N', size, '-v']
                                    logging.info(
                                        'Running command: {}'.format(' '.join(cmd)))
                                    try:
                                        # print(cmd)
                                        out = subprocess.check_output(cmd)
                                        with open(os.path.join(stencil_path, stencil_name.split('.')[0] + '-' + size + '-'  + machine.split('.')[0] + '.txt'), 'wb') as f:
                                            f.write(out)
                                    except subprocess.CalledProcessError as e:
                                        #print("kerncraft failed:", e)
                                        logging.error(
                                            'Failed to execute {}: {}'.format(cmd, e))
                                        logging.error('It returned: {}'.format(out))
                                        sys.exit(1)
                                    # blocksize = 32

                                    param_values.append("{} {}".format(size, size))

                                # #run stempel bench to create actual C code
                                cmd = ['stempel', 'bench', os.path.join(stencil_path, stencil_name), '-m', os.path.join(
                                    machinefilepath, machine), '--store', '--nocli']  # , '-b', blocksize]
                                logging.info(
                                    'Running command: {}'.format(' '.join(cmd)))
                                try:
                                   # print(cmd)
                                    subprocess.check_output(cmd)
                                except subprocess.CalledProcessError as e:
                                    #print("Run failed:", e)
                                    logging.error(
                                        'Failed to execute {}: {}'.format(cmd, e))
                                    logging.error('It returned: {}'.format(out))
                                    sys.exit(1)
                                logging.info('Successfully created benchmark file: {}{}'.format(
                                    stencil_name.split('.')[0], '_compilable.c'))

                                if withprova:
                                    # run the code through prova!
                                    mystencilname = stencil_name.split('.')[0]
                                    project = mystencilname
                                    mystencilname = mystencilname + '_compilable.c'
                                    
                                    params = 'M_MAX N_MAX'
                                    values = '{} {}'.format(sizes[0], sizes[0])
                                    threads = 2
                                    logging.info('Using: {}'.format(mystencilname))
                                    logging.info('Running: {}'.format(project))
                                    logging.info('Running with sizes: {}'.format(param_values))

                                    run_prova(stencil_path, mystencilname, provapath, provaworkspace, likwid_inc, likwid_lib, project, params, values, threads,
                                        method_type, method_name, executions, param_values, exp_threads, pinning)
                        else:  # d == 3
                            for machine in machinefiles:
                                logging.info('Retrieving machinefile')
                                with open(os.path.join(machinefilepath, machine)) as myFile:
                                    results = yaml.load(
                                        myFile, Loader=yaml.RoundTripLoader)
                                try:
                                    dimension = results['memory hierarchy'][
                                        2]['size per group']
                                    dimension = dimension.partition('.')[0]
                                except:
                                    line = results['memory hierarchy'][
                                        2]['cache per group']
                                    sets = line['sets']
                                    ways = line['ways']
                                    cl_size = line['cl_size']
                                    dimension = int(
                                        sets) * int(ways) * int(cl_size) / 1048576
                                    dimension = str(
                                        dimension).partition(' ')[0]
                                size = float(dimension) * 1000000 / 8
                                size = int(round(size ** (1. / 3)))
                                half_size = str(int(size / 2))
                                quarter_size = str(int(size / 4))
                                double_size = str(int(size * 2))
                                size = str(int(size))
                                sizes = [quarter_size, half_size, size, double_size]
                                param_values = []
                                for size in sizes:
                                    if iaca:
                                        ECM = 'ECM'
                                    else:
                                        ECM = 'ECMData'
                                    cmd = ['kerncraft', '-p', cache_model, '-p', 'Roofline', '-p', ECM, os.path.join(stencil_path, stencil_name), '-m', os.path.join(
                                        machinefilepath, machine), '-D', 'M', size, '-D', 'N', size, '-D', 'P', size, '-v']
                                    logging.info(
                                        'Running command: {}'.format(' '.join(cmd)))
                                    try:
                                        # print(cmd)
                                        out = subprocess.check_output(cmd)
                                        with open(os.path.join(stencil_path, stencil_name.split('.')[0] + '-' + size + '-' + machine.split('.')[0] + '.txt'), 'wb') as f:
                                            f.write(out)
                                    except subprocess.CalledProcessError as e:
                                        #print("kerncraft failed:", e)
                                        logging.error(
                                            'Failed to execute {}: {}'.format(cmd, e))
                                        logging.error('It returned: {}'.format(out))
                                        sys.exit(1)

                                    #param_values = param_values + ' "{} {} {}"'.format(size, size, size)
                                    param_values.append("{} {} {}".format(size, size, size))

                                # #run stempel bench to create actual C code
                                cmd = ['stempel', 'bench', os.path.join(stencil_path, stencil_name), '-m', os.path.join(
                                    machinefilepath, machine), '--store', '--nocli']  # , '-b', blocksize]
                                logging.info(
                                    'Running command: {}'.format(' '.join(cmd)))
                                try:
                                    # print(cmd)
                                    subprocess.check_output(cmd)
                                except subprocess.CalledProcessError as e:
                                    #print("Run failed:", e)
                                    logging.error(
                                        'Failed to execute {}: {}'.format(cmd, e))
                                    logging.error('It returned: {}'.format(out))
                                    sys.exit(1)
                                logging.info('Successfully created benchmark file: {}{}'.format(
                                    stencil_name.split('.')[0], '_compilable.c'))


                                if withprova:
                                    # run the code through prova!
                                    mystencilname = stencil_name.split('.')[0]
                                    project = mystencilname
                                    mystencilname = mystencilname + '_compilable.c'

                                    params = 'M_MAX N_MAX P_MAX'
                                    values = '{} {} {}'.format(
                                        sizes[0], sizes[0], sizes[0])
                                    threads = 2

                                    logging.info('Using: {}'.format(mystencilname))
                                    logging.info('Running: {}'.format(project))
                                    logging.info('Running with sizes: {}'.format(param_values))
                                    # run the code through prova!
                                    run_prova(stencil_path, mystencilname, provapath, provaworkspace, likwid_inc, likwid_lib, project, params, values, threads,
                                              method_type, method_name, executions, param_values, exp_threads, pinning)


def run_prova(stencil_path, stencil_name, provapath, provaworkspace, likwid_inc, likwid_lib, project, params, values, threads, method_type, method_name, executions, param_values, exp_threads, pinning):
    filesrc = os.path.join(stencil_path, stencil_name)
    filekernelsrc = os.path.join(stencil_path, 'kernel.c')
    destination_dir = os.path.join(provaworkspace, project,
                                   method_name, 'src')
    filedest = os.path.join(destination_dir, 'example.c')
    filekerneldest = os.path.join(destination_dir, 'kernel.c')
    logging.info('Later {} will be copied to {}'.format(filesrc, filedest))

    basesetup(provapath, provaworkspace, likwid_inc, likwid_lib)
    create_project(provapath, provaworkspace, project, params, values, threads)
    create_method(provapath, provaworkspace, project, method_type, method_name)
    try:
        copyfile(filesrc, filedest)
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            filesrc, filedest, e))

    try:
        copyfile(filekernelsrc, filekerneldest)
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            filekernelsrc, filekerneldest, e))

    headers_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    headers_path = os.path.join(headers_path, 'headers')
    copyheaders(headers_path, destination_dir)

    run_exp(provapath, provaworkspace, project, executions, param_values,
            method_name, exp_threads, pinning)

    exp_dir = getlast_dir(os.path.join(provaworkspace, project, 'experiments'))
    exp_dir_name = os.path.basename(os.path.normpath(exp_dir))

    build_graph(project, exp_dir_name, param_values,
                method_name, exp_threads, pinning)

    try:
        outfile = os.path.join(exp_dir, 'results.json')
    except OSError as e:
        logging.error('Unable to correctly retrieve the last experiment.')

    try:
        logging.info('Copying file: results.json to {}'.format(stencil_path))
        copyfile(outfile, os.path.join(stencil_path, 'results.json'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, 'results.dat')
    try:
        logging.info('Copying file: results.dat')
        copyfile(outfile, os.path.join(stencil_path, 'results.dat'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, 'graph.svg')
    try:
        logging.info('Copying file: graph.svg')
        copyfile(outfile, os.path.join(stencil_path, 'graph.svg'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, 'gnuplot.gp')
    try:
        logging.info('Copying file: gnuplot.gp')
        copyfile(outfile, os.path.join(stencil_path, 'gnuplot.gp'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, '.experiment')
    try:
        logging.info('Copying file: .experiment')
        copyfile(outfile, os.path.join(stencil_path, '.experiment'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))


def main():

    logging.info('=== Started ===')
    """This method is the main, it creates a paerser, uses it and runs the
    business logic
    """

    # Create and populate parser
    parser = create_parser()

    # # Parse given arguments
    args = parser.parse_args()

    # # BUSINESS LOGIC
    args.func(args, parser)

    logging.info('Finished')


if __name__ == '__main__':
    main()
