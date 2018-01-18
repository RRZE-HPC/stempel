#!/usr/bin/env python3
"""
Simple tool to fix the output of likwid_bench_auto.
:author: Danilo Guerrera (University of Basel)
"""

import os
from itertools import islice
import sys
import argparse

def create_parser():
    """This method creates a parser
    """
    example_gen = 'Example usage: fix_yaml -i input.yml -o output.yml'
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=true, help='input filename')
    parser.add_argument('-o', '--output', help='output filename')
    parser.set_defaults(func=fix)

    return parser

def fix(args, parser, output_file=sys.stdout):
    """This method parses the input file and when finds an entry containing
        the PrefixedUnit object reformats it properly
    """
    inputfile = args.input

    if not args.output:
        outputfile = inputfile + '_mod'
    with open(args.input) as f1:
        with open(args.output, 'wb') as f2:
            lines = f1.readlines()
            lit = iter(enumerate(lines))
            for i, line in lit:
                if not '!!python/object:prefixedunit.PrefixedUnit' in line:
                    f2.write(line)
                elif 'size per group' in line:
                    linegroup = lines[i].split(':')[0]
                    prefix = lines[i+1].split(':')[1].strip()
                    if prefix == "''":
                        prefix = ''
                    unit = lines[i+2].split(':')[1].strip()
                    value = lines[i+3].split(':')[1].strip()
                    f2.write(linegroup + ': ' + value + ' ' + prefix + unit + '\n')
                    next(islice(lit, 3,3), None)
                elif 'size per core' in line:
                    firstline = lines[i].split(':')[0]
                    prefix = lines[i+1].split(':')[1].strip()
                    if prefix == "''":
                        prefix = ''
                    unit = lines[i+2].split(':')[1].strip()
                    value = lines[i+3].split(':')[1].strip()
                    f2.write(firstline + ': ' + value + ' ' + prefix + unit + '\n')
                    next(islice(lit, 3,3), None)
                elif 'byte' in line:
                    linebytes = lines[i].split(':')[0]
                    prefix = lines[i+1].split(':')[1].strip()
                    if prefix == "''":
                        prefix = ''
                    unit = lines[i+2].split(':')[1].strip()
                    value = lines[i+3].split(':')[1].strip()
                    linestreams = lines[i+4].split(':')[0]
                    streams = lines[i+4].split(':')[1].strip()
                    f2.write(linebytes + ': ' + value + ' ' + prefix + unit + '\n')
                    f2.write(linestreams + ': ' + streams + '\n')
                    next(islice(lit, 4,4), None)
                else:
                    firstline = lines[i].split('-')[0]
                    prefix = lines[i+1].split(':')[1].strip()
                    if prefix == "''":
                        prefix = ''
                    unit = lines[i+2].split(':')[1].strip()
                    value = lines[i+3].split(':')[1].strip()
                    f2.write(firstline + '- ' + value + ' ' + prefix + unit + '\n')
                    next(islice(lit, 3,3), None)

def main():
    """This method is the main, it creates a paerser, uses it and runs the
    business logic
    """

    # Create and populate parser
    parser = create_parser()

    # # Parse given arguments
    args = parser.parse_args()

    # # BUSINESS LOGIC
    args.func(args, parser)


if __name__ == '__main__':
    main()
