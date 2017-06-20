#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

# Version check
import sys
if sys.version_info[0] == 2 and sys.version_info < (2, 7) or \
        sys.version_info[0] == 3 and sys.version_info < (3, 4):
    print("Must use python 2.7 or 3.4 and greater.", file=sys.stderr)
    sys.exit(1)

import argparse
import os.path
import pickle
import shutil
import math
import re
import itertools
import operator
from functools import reduce

import sympy
import six
from six.moves import range

import stencils




def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-D', '--dimensions', metavar=('DIMENSIONS'), type=int, default=2, required=True,
                        help='Define the number of dimensions to create in the final C code. Values must be integer.')

    parser.add_argument('-S', '--simmetricity', nargs='*', metavar=('SIMMETRICITY'), type=tuple, default="111 111", required=True,
                        help='Define the size of the previously specified dimensions. Values must be integer.')

    parser.add_argument('-k', '--kind', metavar=('VALUE'), choices=stencils.__all__, type=str, default='star',
                        help='Kind of stencil to generate. Value must be star or box')

    parser.add_argument('-C', '--coefficient',  metavar=('COEFF'), type=str, default='scalar', choices=['scalar', 'matrix'], 
                        help='Define if the stencil has a fixed coeffient or a matrix of coefficients. Value must be scalar or matrix')

    parser.add_argument('-t', '--datatype', metavar=('DATATYPE'), type=str, choices=['float', 'double'], default='double',
                        help='Define the datatype of the grids used in the stencil. Value must be double or float')

    # for s in stencils.__all__:
    #     ag = parser.add_argument_group('arguments for '+s+' stencil', getattr(stencils, s).name)
    #     getattr(stencils, s).configure_arggroup(ag)

    return parser

def check_arguments(args, parser):
    if args.coefficient not in ['scalar', 'matrix']:
        parser.error('--coefficient can only be "scalar" or "matrix"')
    if args.datatype not in ['float', 'double']:
        parser.error('--coefficient can only be "float" or "double"')
    if len(args.simmetricity) != args.dimensions:
        parser.error('--simmetricty cannot accept a number of simmetricity levels different from dimensions')
    args.simmetricity = [tuple(x for x in y if x != ',') for y in args.simmetricity]

    

def run(parser, args):

    # Create a new Stencil
    stencil = stencils.Star(dimensions=args.dimensions, simmetricity=args.simmetricity, kind=args.kind, coeff=args.coefficient , datatype =args.datatype)

    # get the maximum diameter
    a = []
    for i in range(0,len(stencil.simmetricity)):
        a.append(len(stencil.simmetricity[i]))
    myrange = max(a)

    # build the name of the output file according to dimensions and diameter
    output_file = '{}d-{}pt.c'.format(stencil.dimensions, myrange)

    # create the declaration part of the final C code
    declaration = stencil.declaration()
    code = declaration
    # create the loop part of the final C code
    loop = stencil.loop()
    code = code + '\n'.join(loop)

    # Save storage to file
    with open(output_file, 'w') as f:
        f.write(code)


def main():
    # Create and populate parser
    parser = create_parser()

    # Parse given arguments
    args = parser.parse_args()

    # Checking arguments
    check_arguments(args, parser)

    # BUSINESS LOGIC IS FOLLOWING
    run(parser, args)


if __name__ == '__main__':
    main()
