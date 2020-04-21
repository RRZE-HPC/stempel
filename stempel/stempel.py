#!/usr/bin/env python
"""This module analyzes and generates C code out of a stencil specification
done via command line
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division


import sys
import argparse
import importlib
import subprocess
import shutil
import os
import itertools

import sympy
from ruamel import yaml

from kerncraft.pycparser_utils import clean_code


from kerncraft.machinemodel import MachineModel
from kerncraft.kerncraft import AppendStringRange
from kerncraft.kernel import Kernel
from kerncraft.kernel import KernelCode

from .benchkernel import KernelBench
from .utilities import count_ops
from . import __version__
# Version check
if sys.version_info[0] == 2 and sys.version_info < (2, 7) or \
        sys.version_info[0] == 3 and sys.version_info < (3, 4):
    print("Must use python 2.7 or 3.4 and greater.", file=sys.stderr)
    sys.exit(1)


def class_for_name(module_name, class_name):
    """
    This method returns a class given its name and the module it belongs to
    """
    # load the module, will raise ImportError if module cannot be loaded
    mod = importlib.import_module(module_name)
    # get the class, will raise AttributeError if class cannot be found
    myclass = getattr(mod, class_name)
    return myclass


def print_header(args, output_file, stencil):
    """
    This method prints all the infos regarding how the tool has been called
    """
    # print header (taken from kerncraft)
    print('{:=^80}'.format(' stempel '), file=output_file)
    print('{:<40}{:>40}'.format('-D ', args.dimensions), file=output_file)
    print('{:<40}{:>40}'.format('-r ', args.radius), file=output_file)
    print('--{}'.format(args.classification), file=output_file)
    print('{:<40}{:>40}'.format('-k ', args.kind), file=output_file)
    print('{:<40}{:>40}'.format('-t ', args.datatype), file=output_file)
    print('{:<40}{:>40}'.format('-C ', args.coefficient), file=output_file)
    print('{:-^80}'.format(' ' + stencil.name + ' ' + args.classification
                           + ' stencil '), file=output_file)


def call_kerncraft(code='', verbosity=''):
    """
    this method tries to call kerncraft on the C like code generated in the
    first place. It calls it in order to get an analysis via ECM of the code
    provided as input on the machine. In order to make it simplier, at the
    moment is not possible to pass a machine file to it, but it will take a
    default one.
    """
    verb = ''
    if verbosity > 0:
        verb = ' -v'
    if verbosity > 1:
        verb += 'v'
    try:
        model = subprocess.Popen(
            ['kerncraft -p ECMData -m phinally.yaml {} -D N 10000 -D M '
             '10000{}'.format(code, verb)], stdout=subprocess.PIPE
        ).communicate()[0].decode("utf-8")
    except OSError:
        print('likwid-topology execution failed, is it installed and loaded?',
              file=sys.stderr)
        sys.exit(1)

    return model


def create_parser():
    """This method creates a parser
    """
    example_gen = 'Example usage: stempel gen -D 2 -r 1 -i -k star -C variable -d 2'
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))

    subparsers = parser.add_subparsers(dest='subparser_name',
                                       help='sub-command help')
    subparsers.required = True
    parser_gen = subparsers.add_parser('gen', help='Generate a C-like code for '
                                       'the stencil computation descrivbed by '
                                       'the command line parameters',
                                       epilog=example_gen)
    parser_gen.add_argument('-D', '--dimensions', metavar=('DIMENSIONS'),
                            type=int, default=2, required=True, help='Define '
                            'the number of dimensions to create in the final C '
                            'code. Value must be integer.')

    parser_gen.add_argument('-r', '--radius', metavar=('RADIUS'), type=int,
                            default=2, required=True, help='Define the radius of '
                            'the stencil in each dimension. Value must be integer.')

    classification = parser_gen.add_mutually_exclusive_group(required=True)
    classification.add_argument('-c', '--classification', dest='classification',
                                choices=['isotropic', 'heterogeneous', 'homogeneous',
                                         'point-symmetric'],
                                help="Define the coefficitions to be used.")
    classification_group = classification.add_mutually_exclusive_group()
    classification_group.add_argument('-i', '--isotropic', action='store_const',
                                      dest='classification', const='isotropic',
                                      help='Define if the coefficients are '\
                                      'isotropic (does not depend on the direction).')
    classification_group.add_argument('-e', '--heterogeneous', action='store_const',
                                      dest='classification', const='heterogeneous',
                                      help='Define if the weighting factors expose'\
                                      ' no symmetry, i.e. a different coefficient '\
                                      'for each direction')
    classification_group.add_argument('-o', '--homogeneous', action='store_const',
                                      dest='classification', const='homogeneous',
                                      help='Define if the stencil has a homogenuous'\
                                      ' coefficient (i.e. only a sincgle scalar)')
    classification_group.add_argument('-p', '--point-symmetric',
                                      action='store_const', dest='classification',
                                      const='point-symmetric', help='Define if the'\
                                      ' weighting factors are symmetric with respect'\
                                      ' to the origin, in which case, they have '\
                                      'the same value.')

    parser_gen.add_argument('-k', '--kind', choices=['star', 'box'],
                            type=str, default='star',
                            help='Kind of stencil to generate. Value must be star '
                            'or box')

    parser_gen.add_argument('-C', '--coefficient', type=str,
                            default='constant', choices=['constant', 'variable'],
                            help='Define if the stencil has a fixed coeffient or a'
                            ' matrix of coefficients. Value must be scalar or '
                            'matrix')

    parser_gen.add_argument('-t', '--datatype', type=str,
                            choices=['float', 'double', 'float _Complex', 'double _Complex'],
                            default='double',
                            help='Define the datatype of the grids used in the '
                            'stencil. Value must be double, float, "float _Complex" or "double '
                            '_Complex"')

    parser_gen.add_argument('-d', '--dimofcoeffs', type=int,
                            help='Variable coefficients are stored as an array. '
                            'Each variable coefficient results in an array of the '
                            'same size of the stencils\' grid. All variable '
                            'coefficients are packed into an array (e.g. W[M][N][2]'
                            '). This parameter defines which dimension (first, '
                            ' second, third, ...) will be used for discerning the '
                            'coefficients. Value must be an integer.')

    parser_gen.add_argument('--store', type=argparse.FileType('w'),
                            help='Addes results to a C file for later processing.')

    parser_gen.add_argument('--verbose', '-v', action='count', default=0,
                            help='Increases verbosity level.')

    parser_gen.set_defaults(func=run_gen)

    parser_bench = subparsers.add_parser('bench', help='Generate a benchmark '
                                         'code for the stencil kernel passed as'
                                         'input. It must be a C source code.')

    parser_bench.add_argument('code_file', metavar='FILE',
                              type=argparse.FileType(),
                              help='File with declarations and C loop kernel')

    parser_bench.add_argument('--machine', '-m', type=argparse.FileType('r'),
                              required=True, help='Path to machine description '
                              'yaml file.')
    parser_bench.add_argument('--block', '-b', nargs='?', type=int, const=1, default=0,
                              help='Blocking factor:\n'\
                              '0:  no blocking\n'\
                              '1 (default):  blocking for the middle (3D) or outermost (2D) loop\n'\
                              '>1: full blocking in all three dimensions')
    parser_bench.add_argument('-D', '--define', nargs=2, metavar=('KEY', 'VALUE'), default=[],
                              action=AppendStringRange,
                              help='Define constant to be used in C code. Values '\
                              'must be integer')
    parser_bench.add_argument('--nocli', action='store_true', default=False,
                                      help='Define wehther to generate the '\
                                      'version of the code to be run via '\
                                      'command line interface or via PROVA')
    parser_bench.add_argument('--store', action='store_true',
                              help='Addes results to a C file for later processing.')
    parser_bench.add_argument('--initwithrand', action='store_true', default=False,
                              help='Initialize the arrays with random numbers.')

    parser_bench.set_defaults(func=run_bench)
    # for s in stencils.__all__:
    #     ag = parser.add_argument_group('arguments for '+s+' stencil',
    #         getattr(stencils, s).name)
    #     getattr(stencils, s).configure_arggroup(ag)

    return parser


def check_arguments(args, parser):
    """This method checks that some of the arguments given to the parser
    respect our convention (are what they are supposed to be)
    """
    if args.datatype not in ['float', 'double', 'float _Complex', 'double _Complex']:
        parser.error('--coefficient can only be "float", "double", "float _Complex" or '
                     '"double _Complex".')

    if args.kind == 'box':
        if args.dimensions == 2 and args.radius > 6:
            parser.error('Please choose a different combination of values for '
                         'dimensions and radius, because this generates too many coefficients.')
        elif args.dimensions == 3 and args.radius > 5:
            parser.error('Please choose a different combination of values for '
                         'dimensions and radius, because this generates too many coefficients.')


def change_decl(decltoreplace, dimofcoeffs, stencil, code):
    """
    This method accept a declaration as a string, the stencil, a new position
    and the code. It returns the code after having rearranged the declaration:
    it shifts the dimension holding the coefficients.
    Example: W[i][j][3] --> W[i][3][j]
    """
    size = str(stencil.num_coefficients)
    newdecl = decltoreplace.replace('[' + size + ']', '')
    letter_length = 0
    for i in range(dimofcoeffs - 1):
        letter_length += len(stencil.dims[i])

    pos = 1 + 2 * (dimofcoeffs - 1) + letter_length

    newdecl = newdecl[:pos] + '[' + size + ']' + newdecl[pos:]

    code = code.replace(decltoreplace, newdecl)
    return code


def change_loop(decltoreplace, dimofcoeffs, stencil, code):
    """
    This method accept a declaration as a string, the stencil, a new position
    and the code. It returns the code after having rearranged the loop:
    it shifts the dimension holding the coefficients.
    Example: W[i][j][3] --> W[i][3][j]
    """
    assert isinstance(code, str) or isinstance(code, unicode), "The type of"\
        "the code variable must be either unicode or string"

    array_of_coeffs = [stencil.coefficients[0] + '[' +
                       str(i) + ']' for i in range(stencil.num_coefficients)]
    newarray = []
    for i in array_of_coeffs:
        # remove trailing dimension
        newi = i.rsplit('[', 1)[0]

        size = i.rsplit('[', 1)[1].rsplit(']', 1)[0]

        letter_length = 0
        for i in range(dimofcoeffs - 1):
            letter_length += len(stencil.loop_variables[i])
        pos = 1 + 2 * (dimofcoeffs - 1) + letter_length

        newi = newi[:pos] + '[' + size + ']' + newi[pos:]
        newarray.append(newi)

    for i in range(stencil.num_coefficients):
        code = code.replace(array_of_coeffs[i], newarray[i])

    return code


def run_gen(args, parser, output_file=sys.stdout):
    """This method creates an object of type Stencil and calls the appropriate
    methods to generate the C code
    """
    # Checking arguments
    check_arguments(args, parser)
    # Create a new Stencil
    # first we need to retrive the name of the stencil class out of the "kind"
    # passed via command line
    mykind = (args.kind).title()
    if(args.coefficient) == 'constant':
        mykind = mykind + 'Constant'
    elif(args.coefficient) == 'variable':
        mykind = mykind + 'Variable'

    # with tox
    stencil_class = class_for_name('stempel.stencils', mykind)
    # without tox
    #stencil_class = class_for_name('stencils', mykind)

    stencil = stencil_class(dimensions=args.dimensions, radius=args.radius,
                            classification=args.classification,
                            datatype=args.datatype)

    # create the declaration part of the final C code
    declaration = stencil.declaration()
    # create the loop part of the final C code
    loop, flop = stencil.loop()

    code = declaration + '\n'.join(loop)

    if args.coefficient == 'variable':

        toreplace = stencil.coefficients[0]

        decltoreplace = toreplace
        for i in range(args.dimensions):
            decltoreplace = decltoreplace.replace(
                stencil.loop_variables[i], stencil.dims[i])

        decltoreplace = decltoreplace + \
            '[' + str(stencil.num_coefficients) + ']'

        # assert args.dimofcoeffs < len(stencil.dims)+1, "The stencil provided"\
        # "does not have enough dimensions to place the index of the coefficients"\
        # "at position {}".format(args.dimofcoeffs)
        if args.dimofcoeffs and args.dimofcoeffs < len(stencil.dims) + 1:
            code = change_decl(decltoreplace, args.dimofcoeffs, stencil, code)
            code = change_loop(decltoreplace, args.dimofcoeffs, stencil, code)

    print_header(args, output_file, stencil)

    # Save storage to file or print to STDOUT
    if args.store:
        args.store.write(code)
    else:
        print(code, file=output_file)
        print('The kernel consists of {} FLOP'.format(flop), file=output_file)

    # if verbose print a little bit more infos
    if args.verbose > 0:
        for arg in vars(args):
            print('{:<40}{:>40}'.format(arg, str(getattr(args, arg))),
                  file=output_file)


def run_bench(args, output_file=sys.stdout):
    """This method creates a benchmark code starting from the C code passed as
    input
    """
    # process kernel

    # machine information
    # Read machine description
    machine = MachineModel(args.machine.name)  # , args=args)
    #compiler, compiler_args = machine.get_compiler()

    code = str(args.code_file.read())
    code = clean_code(code)

    flop = count_ops(code)

    kernel = KernelBench(code, filename=args.code_file.name,
                         machine=machine, block_factor=args.block,
                         flop=flop, initwithrand=args.initwithrand )

    # taken from kerncraft
    # # works only for up to 3 dimensions
    array = []
    required_consts = []
    for v in kernel.variables.values():
        array = []
        if v[1] is not None:
            for c in v[1]:
                if type(c) is not sympy.Integer:
                    array.append(c)
        required_consts.append(array)

    # required_consts = [v[1] for v in kernel.variables.values() if v[1] is not None]
    required_consts = set([i for l in required_consts for i in l])

    if len(required_consts) > 0:
        define_dict = {}
        if args.define:
            # build defines permutations
            for name, values in args.define:
                if name not in define_dict:
                    define_dict[name] = [[name, v] for v in values]
                    continue
                for v in values:
                    if v not in define_dict[name]:
                        define_dict[name].append([name, v])
        else:
            my_constant_size = 100
            for name in list(required_consts):
                if name not in define_dict:
                    define_dict[name] = [[name, my_constant_size]]
        define_product = list(itertools.product(*list(define_dict.values())))

        # Check that all consts have been defined
        # if set(required_consts).difference(set([symbol_pos_int(k) for k in define_dict.keys()])):
        #     raise ValueError("Not all constants have been defined. Required are: {}".format(
        #         required_consts))
    else:
        define_product = [{}]

    # if define_product:
    #     size = define_product[0][0][1]
    # else:
    #     size = 100

    for define in define_product:
        # Reset state of kernel
        kernel.clear_state()

        # Add constants from define arguments
        for k, v in define:
            kernel.set_constant(k, v)

    # get compilable C code
    c_code, kernel = kernel.as_code(from_cli=not args.nocli)

   # Save storage to file or print to STDOUT
    if args.store:
        # build the name of the output file according to dimensions and
        # diameter
        tempname = args.code_file.name + '.tmp'

        with open(tempname, 'w') as out:
            out.write(c_code)
        shutil.move(tempname, args.code_file.name.split('.')[0] + "_compilable.c")
        with open(os.path.join(os.path.dirname(tempname), 'kernel.c'), 'w') as out:
            out.write(kernel)
    else:
        print(c_code, file=output_file)
        print(kernel, file=output_file)


def main():
    """This method is the main, it creates a paerser, uses it and runs the
    business logic
    """
    # Create and populate parser
    parser = create_parser()

    # Parse given arguments
    args = parser.parse_args()

    # BUSINESS LOGIC
    args.func(args, parser)

if __name__ == '__main__':
    main()
