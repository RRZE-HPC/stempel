import os
import errno
import math
import subprocess
import sys
import logging
import argparse
from shutil import copyfile
from glob import glob

from ruamel import yaml
from kerncraft import kerncraft
from kerncraft.likwid_bench_auto import get_machine_topology

#from . import stempel


logging.basicConfig(level=logging.INFO, filename='/tmp/myanalysis.log',
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
    parser.add_argument('-k', '--kind', metavar=('KIND'),
                        choices=['star', 'box'], help='Kind of stencil')
    parser.add_argument('-m', '--machinepath', metavar=('MACHINEPATH'),
                        required=True, help='Path to the machinefile(s)')
    parser.add_argument('-r', '--radius', type=int, metavar=('RADIUS'),
                        help='Radius of the stencil')
    parser.add_argument('-d', '--dimensions', metavar=('DIMENSIONS'),
                        type=int, help='Number of dimensions of the stencil')
    parser.add_argument('-e', '--executions', metavar=('EXECUTIONS'),
                        type=int, default=5, help='Number of dimensions of the stencil')
    parser.add_argument('-t', '--threads', metavar=('THREADS'),
                        nargs='*', help='Number of dimensions of the stencil')
    parser.add_argument('--method_type', metavar=('METHODTYPE'),
                        default='OpenMP-4.0-GCC-4.9.3-2.25', help='Number of dimensions of the stencil')
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

def create_project(provapath, project, params, values, threads):
    cmd_prefix = ['source', os.path.join(
        provapath, 'util', 'BaseSetup.sh'), '&&']

    # create a project
    cmd = cmd_prefix + ['workflow', 'project', '-c', '-p', project,
                        '--params', params, '--values', values, '--threads', str(threads)]
    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        #print("Run failed:", e)
        logging.error("Run failed: {}".format(e))
        sys.exit(1)


def create_method(provapath, project, method_type, method_name):
    cmd_prefix = ['source', os.path.join(
        provapath, 'util', 'BaseSetup.sh'), '&&']

    # create a method
    cmd = cmd_prefix + ['workflow', 'method', '-c', '-p', project,
                        '-m', method_type, '-n', method_name]
    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        #print("Run failed:", e)
        logging.error("Run failed: {}".format(e))
        sys.exit(1)


def run_exp(provapath, project, executions, param_values, method_name, threads, pinning):
    cmd_prefix = ['source', os.path.join(
        provapath, 'util', 'BaseSetup.sh'), '&&']

    # compile the code and run an experiment of the method
    cmd = cmd_prefix + ['workflow', 'run_exp', '-p', project, '-e', str(executions),
                        '-d', param_values, '-m', method_name, '-t', str(threads), '--pin', pinning]
    try:
        logging.info('Running command: {}'.format(' '.join(cmd)))
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError as e:
        #print("Run failed:", e)
        logging.error("Run failed: {}".format(e))
        sys.exit(1)


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
    if args.prova:
        provapath, provaworkspace = args.prova
        withprova = True
    else:
        withprova = None

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
        else:
            radius = radius_box
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
                            #print(cmd)
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
                                    results = yaml.load(myFile, Loader=yaml.RoundTripLoader)
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
                                size = str(int(math.sqrt(size)))
                                if iaca:
                                    ECM = 'ECM'
                                else:
                                    ECM = 'ECMData'
                                cmd = ['kerncraft', '-p', 'LC', '-p', 'Roofline', '-p', ECM, os.path.join(
                                    stencil_path, stencil_name), '-m', os.path.join(machinefilepath, machine), '-D', 'M', size, '-D', 'N', size]
                                logging.info(
                                    'Running command: {}'.format(' '.join(cmd)))
                                try:
                                    #print(cmd)
                                    out = subprocess.check_output(cmd)
                                    with open(os.path.join(stencil_path, stencil_name.split('.')[0] + '-' + machine.split('.')[0] + '.txt'), 'w') as f:
                                        f.write(out)
                                except subprocess.CalledProcessError as e:
                                    #print("kerncraft failed:", e)
                                    logging.error(
                                        'Failed to execute {}: {}'.format(cmd, e))
                                    sys.exit(1)
                                # blocksize = 32
                                # #run stempel bench to create actual C code
                                cmd = ['stempel', 'bench', os.path.join(stencil_path, stencil_name), '-m', os.path.join(
                                    machinefilepath, machine), '-D', size, '-D', size]  # , '-b', blocksize]
                                logging.info(
                                    'Running command: {}'.format(' '.join(cmd)))
                                try:
                                   #print(cmd)
                                    subprocess.check_output(cmd)
                                except subprocess.CalledProcessError as e:
                                    #print("Run failed:", e)
                                    logging.error(
                                        'Failed to execute {}: {}'.format(cmd, e))
                                    sys.exit(1)
                                logging.info('Successfully created benchmark file: {}{}'.format(
                                    stencil_name.split('.')[0], '_compilable.c'))
                                
                                if withprova:
                                    # run the code through prova!
                                    project = stencil_name.split('.')[0]
                                    params = '"M N"'
                                    values = '"{} {}"'.format(size, size)
                                    param_values = values
                                    threads = 2

                                    run_prova(stencil_path, stencil_name, provapath, provaworkspace, project, params, values, threads,
                                              method_type, method_name, executions, param_values, exp_threads, pinning)
                        else:  # d == 3
                            for machine in machinefiles:
                                logging.info('Retrieving machinefile')
                                with open(os.path.join(machinefilepath, machine)) as myFile:
                                    results = yaml.load(myFile, Loader=yaml.RoundTripLoader)
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
                                size = str(int(round(size ** (1. / 3))))
                                if iaca:
                                    ECM = 'ECM'
                                else:
                                    ECM = 'ECMData'
                                cmd = ['kerncraft', '-p', 'LC', '-p', 'Roofline', '-p', ECM, os.path.join(stencil_path, stencil_name), '-m', os.path.join(
                                    machinefilepath, machine), '-D', 'M', size, '-D', 'N', size, '-D', 'P', size]
                                logging.info(
                                    'Running command: {}'.format(' '.join(cmd)))
                                try:
                                    #print(cmd)
                                    out = subprocess.check_output(cmd)
                                    with open(os.path.join(stencil_path, stencil_name.split('.')[0] + '-' + machine.split('.')[0] + '.txt'), 'w') as f:
                                        f.write(out)
                                except subprocess.CalledProcessError as e:
                                    #print("kerncraft failed:", e)
                                    logging.error(
                                        'Failed to execute {}: {}'.format(cmd, e))
                                    sys.exit(1)
                                # blocksize = 32
                                # #run stempel bench to create actual C code
                                cmd = ['stempel', 'bench', os.path.join(stencil_path, stencil_name), '-m', os.path.join(
                                    machinefilepath, machine), '-D', 'M', size, '-D', 'N', size, '-D', 'P', size, '--store']  # , '-b', blocksize]
                                logging.info(
                                    'Running command: {}'.format(' '.join(cmd)))
                                try:
                                    #print(cmd)
                                    subprocess.check_output(cmd)
                                except subprocess.CalledProcessError as e:
                                    #print("Run failed:", e)
                                    logging.error(
                                        'Failed to execute {}: {}'.format(cmd, e))
                                    sys.exit(1)
                                logging.info('Successfully created benchmark file: {}{}'.format(
                                    stencil_name.split('.')[0], '_compilable.c'))

                                if withprova:
                                    # machine = get_machine_topology()
                                    # actual_machine_name = machine['model name'].replace(' ', '_').replace('(TM)', '') + '.yml'
                                    # if not actual_machine_name in machinefiles:
                                    #     cmd = ['python', 'likwid_bench_auto']
                                    # logging.info(
                                    #     'Running command: {}'.format(' '.join(cmd)))
                                    # try:
                                    #     out = subprocess.check_output(cmd)
                                    # except subprocess.CalledProcessError as e:
                                    #     #print("Run failed:", e)
                                    #     logging.error(
                                    #         'Failed to execute {}: {}'.format(cmd, e))
                                    #     sys.exit(1)
                                    # with open(os.path.join(machinepath, actual_machine_name), 'w', encoding='utf8') as outfile:
                                    #     yaml.dump(machine, outfile, default_flow_style=False, allow_unicode=True)
                                    # logging.info('Successfully created machinefile for the current architecture')
                                    

                                    project = stencil_name.split('.')[0]
                                    params = '"M N P"'
                                    values = '"{} {}"'.format(size, size, size)
                                    param_values = values
                                    threads = 2
                                    # run the code through prova!
                                    run_prova(stencil_path, stencil_name, provapath, provaworkspace, project, params, values, threads,
                                              method_type, method_name, executions, param_values, exp_threads, pinning)


def run_prova(stencil_path, stencil_name, provapath, provaworkspace, project, params, values, threads, method_type, method_name, executions, param_values, exp_threads, pinning):
    filesrc = os.path.join(stencil_path, stencil_name)
    filedest = os.path.join(provaworkspace, project,
                            method_name, 'src', 'example.c')
    logging.info('Later {} will be copied to {}'.format(filesrc, filedest))

    create_project(provapath, project, params, values, threads)
    create_method(provapath, project, method_type, method_name)
    #copyfile(filesrc, filedest)
    run_exp(provapath, project, executions, param_values,
            method_name, exp_threads, pinning)

    exp_dir = getlast_dir(os.path.join(provaworkspace, project, 'experiments'))
    try:
        outfile = os.path.join(exp_dir, 'results.json')
    except OSError as e:
        logging.error('Unable to correctly retrieve the last experiment.')

    try:
        copyfile(outfile, os.path.join(stencil_path, 'results.json'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, 'results.dat')
    try:
        copyfile(outfile, os.path.join(stencil_path, 'results.dat'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, 'graph.svg')
    try:
        copyfile(outfile, os.path.join(stencil_path, 'graph.svg'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, 'gnuplot.gp')
    try:
        copyfile(outfile, os.path.join(stencil_path, 'gnuplot.gp'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))

    outfile = os.path.join(exp_dir, '.experiment')
    try:
        copyfile(outfile, os.path.join(stencil_path, '.experiment'))
    except IOError as e:
        logging.error('Unable to copy file {} to {}. {}'.format(
            outfile, stencil_path, e))


def main():

    logging.info('=== Started ===')
    """This method is the main, it creates a paerser, uses it and runs the
    business logic
    """
    # provapath = '/export/hpwc/PROVA'
    # provaworkspace = '~/Desktop/workspace'
    # project = 'try1'
    # params = '"M N"'
    # values = '"200 200"'
    # threads = 2
    # method_type = 'OpenMP-4.0-GCC-4.9.3-2.25'
    # method_name = 'mytry'

    # executions = '5'
    # param_values = '"200 200"'
    # exp_threads = '2 4 8 16'
    # pinning = 'node'

    # d = 2
    # r = 2
    # l = 'isotropic'
    # k = 'star'
    # c = 'constant'
    # #name = str(d) + 'd-' + str(r) + 'r-' + l + '-' + c + '-' + k + '-stencil.c'
    # stencil_name = '2d-2r-isotropic-constant-star-stencil.c'
    # workspace = '~/Desktop/workspace'
    # stencil_path = os.path.join(workspace, 'stencils', str(
    #     d) + 'D', str(r) + 'r', l, k, c)

    # run_prova(stencil_path, stencil_name, provapath, provaworkspace, project, params, values, threads,
    #           method_type, method_name, executions, param_values, exp_threads, pinning)

    # exit(1)
    # Create and populate parser
    parser = create_parser()

    # # Parse given arguments
    args = parser.parse_args()

    # # BUSINESS LOGIC
    args.func(args, parser)

    # kind = ['star', 'box']

    # #path = '/Users/guerrera/Documents/Projects/Stempel/stempel'
    # workspace = '/Users/guerrera/Desktop/stempel_stencils'
    # machinefilepath = os.path.join(workspace, 'machine-files')
    # machinefiles = os.listdir(machinefilepath)
    # stencilfiles = os.path.join(workspace, 'stencils')
    # radius_star = [1, 2, 3, 4]
    # radius_box = [1, 2]
    # dimensions = [2, 3]
    # classification = ['isotropic', 'heterogeneous',
    #                   'homogeneous', 'point-symmetric']

    # coefficients = ['constant', 'variable']
    # iaca = None
    # stencils = []
    # for k in kind:
    #     if k == 'star':
    #         radius = radius_star
    #     else:
    #         radius = radius_box
    #     for c in coefficients:
    #         for l in classification:
    #             for r in radius:
    #                 for d in dimensions:

    #                     name = str(d) + 'd-' + str(r) + 'r-' + l + \
    #                         '-' + c + '-' + k + '-stencil.c'
    #                     logging.info('Working on: {}'.format(name))
    #                     mypath = os.path.join(stencilfiles, str(
    #                         d) + 'D', str(r) + 'r', l, k, c)
    #                     stencils.append(name)
    #                     try:
    #                         os.makedirs(mypath)
    #                     except OSError as e:
    #                         if e.errno != errno.EEXIST:
    #                             raise
    #                     cmd = ['stempel', 'gen', '-D', str(d), '-r', str(
    #                         r), '-k', k, '-C', c, '--' + l, '--store', os.path.join(mypath, name)]
    #                     try:
    #                         print(cmd)
    #                         # subprocess.check_output(cmd)
    #                     except subprocess.CalledProcessError as e:
    #                         print("Run failed:", e)
    #                         sys.exit(1)
    #                     if d == 2:
    #                         logging.info('Retrieving machinefile')
    #                         for machine in machinefiles:
    #                             with open(os.path.join(machinefilepath, machine)) as myFile:
    #                                 results = yaml.safe_load(myFile)
    #                             try:
    #                                 dimension = results['memory hierarchy'][
    #                                     2]['size per group']
    #                                 dimension = dimension.partition('.')[0]
    #                             except:
    #                                 line = results['memory hierarchy'][
    #                                     2]['cache per group']
    #                                 sets = line['sets']
    #                                 ways = line['ways']
    #                                 cl_size = line['cl_size']
    #                                 dimension = int(
    #                                     sets) * int(ways) * int(cl_size) / 1048576
    #                                 dimension = str(
    #                                     dimension).partition(' ')[0]
    #                             size = float(dimension) * 1000000 / 8
    #                             size = str(int(math.sqrt(size)))
    #                             if iaca:
    #                                 ECM = 'ECM'
    #                             else:
    #                                 ECM = 'ECMData'
    #                             cmd = ['kerncraft', '-p', 'LC', '-p', 'Roofline', '-p', ECM, os.path.join(
    #                                 mypath, name), '-m', os.path.join(machinefilepath, machine), '-D', 'M', size, '-D', 'N', size]
    #                             logging.info(
    #                                 'Running command: {}'.format(' '.join(cmd)))
    #                             try:
    #                                 print(cmd)
    #                                 out = subprocess.check_output('ls')
    #                                 with open(os.path.join(mypath, name.split('.')[0] + '-' + machine.split('.')[0] + '.txt'), 'w') as f:
    #                                     f.write(out)
    #                             except subprocess.CalledProcessError as e:
    #                                 print("kerncraft failed:", e)
    #                                 sys.exit(1)
    #                             # blocksize = 32
    #                             # #run stempel bench to create actual C code
    #                             # cmd = ['stempel', 'gen', os.path.join(mypath, name), '-m', os.path.join(machinefilepath, machine), '-D', size, '-D', size, '-b', blocksize]
    #                             # try:
    #                             #     print(cmd)
    #                             #     # subprocess.check_output(cmd)
    #                             # except subprocess.CalledProcessError as e:
    #                             #     print("Run failed:", e)
    #                             #     sys.exit(1)
    #                     else:  # d == 3
    #                         for machine in machinefiles:
    #                             logging.info('Retrieving machinefile')
    #                             with open(os.path.join(machinefilepath, machine)) as myFile:
    #                                 results = yaml.safe_load(myFile)
    #                             try:
    #                                 dimension = results['memory hierarchy'][
    #                                     2]['size per group']
    #                                 dimension = dimension.partition('.')[0]
    #                             except:
    #                                 line = results['memory hierarchy'][
    #                                     2]['cache per group']
    #                                 sets = line['sets']
    #                                 ways = line['ways']
    #                                 cl_size = line['cl_size']
    #                                 dimension = int(
    #                                     sets) * int(ways) * int(cl_size) / 1048576
    #                                 dimension = str(
    #                                     dimension).partition(' ')[0]
    #                             size = float(dimension) * 1000000 / 8
    #                             size = str(int(round(size ** (1. / 3))))
    #                             if iaca:
    #                                 ECM = 'ECM'
    #                             else:
    #                                 ECM = 'ECMData'
    #                             cmd = ['kerncraft', '-p', 'LC', '-p', 'Roofline', '-p', ECM, os.path.join(mypath, name), '-m', os.path.join(
    #                                 machinefilepath, machine), '-D', 'M', size, '-D', 'N', size, '-D', 'P', size]
    #                             logging.info(
    #                                 'Running command: {}'.format(' '.join(cmd)))
    #                             try:
    #                                 print(cmd)
    #                                 out = subprocess.check_output('ls')
    #                                 with open(os.path.join(mypath, name.split('.')[0] + '-' + machine.split('.')[0] + '.txt'), 'w') as f:
    #                                     f.write(out)
    #                             except subprocess.CalledProcessError as e:
    #                                 print("kerncraft failed:", e)
    #                                 sys.exit(1)
    logging.info('Finished')


if __name__ == '__main__':
    main()
