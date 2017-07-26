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

# Version check
if sys.version_info[0] == 2 and sys.version_info < (2, 7) or \
        sys.version_info[0] == 3 and sys.version_info < (3, 4):
    print("Must use python 2.7 or 3.4 and greater.", file=sys.stderr)
    sys.exit(1)


def class_for_name(module_name, class_name):
    """
    this method returns a class given its name and the module it belongs to
    """
    # load the module, will raise ImportError if module cannot be loaded
    mod = importlib.import_module(module_name)
    # get the class, will raise AttributeError if class cannot be found
    myclass = getattr(mod, class_name)
    return myclass


def create_parser():
    """This method creates a parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-D', '--dimensions', metavar=('DIMENSIONS'), type=int,
                        default=2, required=True, help='Define the number of '
                        'dimensions to create in the final C code. Value must '
                        'be integer.')

    parser.add_argument('-r', '--radius', metavar=('RADIUS'), type=int,
                        default=2, required=True, help='Define the radius of '
                        'the stencil in each dimension. Value must be integer.')

    symmetry_group = parser.add_mutually_exclusive_group()
    symmetry_group.add_argument('-s', '--symmetric', action='store_true',
                                help='Define if the coefficient is symmetric '
                                'along the two sides of an axis.')
    symmetry_group.add_argument('-a', '--asymmetric', action='store_true',
                                help='Define if the coefficient is symmetric '
                                'along the two sides of an axis.')

    isotropy_group = parser.add_mutually_exclusive_group()
    isotropy_group.add_argument('-i', '--isotropic', action='store_true',
                                help='Define if the coefficients are isotropic '
                                '(does not depend on the direction).')
    isotropy_group.add_argument('-A', '--anisotropic', action='store_true',
                                help='Define if the coefficients are isotropic '
                                '(does not depend on the direction).')

    parser.add_argument('-k', '--kind', metavar=('KIND'),
                        choices=['star', 'box'], type=str, default='star',
                        help='Kind of stencil to generate. Value must be star '
                        'or box')

    parser.add_argument('-C', '--coefficient', metavar=('COEFF'), type=str,
                        default='constant', choices=['constant', 'variable'],
                        help='Define if the stencil has a fixed coeffient or a'
                        ' matrix of coefficients. Value must be scalar or '
                        'matrix')

    parser.add_argument('-g', '--inputgrids', metavar=('INPUTGRIDS'), type=int,
                        default=1, required=True, help='Define the number of '
                        'input grids to create in the final C code. Value must '
                        'be integer.')

    parser.add_argument('-t', '--datatype', metavar=('DATATYPE'), type=str,
                        choices=['float', 'double'], default='double',
                        help='Define the datatype of the grids used in the '
                        'stencil. Value must be double or float')

    # for s in stencils.__all__:
    #     ag = parser.add_argument_group('arguments for '+s+' stencil', getattr(stencils, s).name)
    #     getattr(stencils, s).configure_arggroup(ag)

    return parser


def check_arguments(args, parser):
    """This method checks that some of the arguments given to the parser
    respect our convention (are what they are supposed to be)
    """
    if args.coefficient not in ['constant', 'variable']:
        parser.error('--coefficient can only be "scalar" or "matrix"')
    if args.datatype not in ['float', 'double']:
        parser.error('--coefficient can only be "float" or "double"')


def run(args):
    """This method creates an object of type Stencil and calls the appropriate
    methods to generate the C code
    """
    # Create a new Stencil
    # first we need to retrive the name of the stencil class out of the "kind"
    # passed via command line
    mykind = (args.kind).title()
    if(args.coefficient) == 'constant':
        mykind = mykind + 'Constant'
    elif(args.coefficient) == 'variable':
        mykind = mykind + 'Variable'

    stencil_class = class_for_name('stencils', mykind)

    if args.asymmetric:
        symmetricity = False
        symmetric = 'asymmetric'
    else:
        symmetricity = True
        symmetric = 'symmetric'

    if args.anisotropic:
        isotropy = False
        isotropic = 'anisotropic'
    else:
        isotropy = True
        isotropic = 'isotropic'

    stencil = stencil_class(dimensions=args.dimensions, radius=args.radius,
                            symmetricity=symmetricity, isotropy=isotropy,
                            datatype=args.datatype, inputgrids=args.inputgrids)

    # build the name of the output file according to dimensions and diameter
    if args.kind == 'box':
        points = (stencil.radius * 2 + 1)**stencil.dimensions
        output_file = '{}d-{}pt-{}.c'.format(stencil.dimensions, points, 'box')
    
    else:
        # diameter = 2 * radius
        # points = diamiter * dimensions + 1 (centerpoint)
        points = 2 * stencil.radius * stencil.dimensions + 1
    
        output_file = '{}d-{}pt-{}-{}-{}_{}.c'.format(stencil.dimensions,
                                                      points,
                                                      isotropic, symmetric,
                                                      args.coefficient,
                                                      'coefficients')

    # create the declaration part of the final C code
    declaration = stencil.declaration()
    # create the loop part of the final C code
    loop = stencil.loop()
    code = declaration + '\n'.join(loop)

    # Save storage to file
    with open(output_file, 'w') as outfile:
        outfile.write(code)


def main():
    """This method is the main, it creates a paerser, uses it and runs the
    business logic
    """
    # Create and populate parser
    parser = create_parser()

    # Parse given arguments
    args = parser.parse_args()

    # Checking arguments
    check_arguments(args, parser)

    # BUSINESS LOGIC IS FOLLOWING
    run(args)

if __name__ == '__main__':
    main()
