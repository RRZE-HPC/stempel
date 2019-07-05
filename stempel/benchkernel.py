#!/usr/bin/env python
"""This module focuses on the generation of C code manipulating the AST
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

from copy import deepcopy
import operator
import tempfile
import subprocess
import os
import os.path
import sys
import collections
import random
from functools import reduce

#import logging
from distutils.spawn import find_executable

import sympy

from pycparser import CParser, c_ast, plyparser
from pycparser.c_generator import CGenerator

import kerncraft
from kerncraft.kernel import KernelCode
from kerncraft.machinemodel import MachineModel

if sys.version_info[0] == 2 and sys.version_info >= (2, 7):
    def fullmatch(regex, string, flags=0):
        """Emulate python-3.4 re.fullmatch()."""
        return re.match("(?:" + regex + r")\Z", string, flags=flags)
elif sys.version_info[0] == 3 and sys.version_info >= (3, 4):
    from re import fullmatch

def symbol_pos_int(*args, **kwargs):
    """Create a sympy.Symbol with positive and integer assumptions."""
    kwargs.update({'positive': True,
                   'integer': True})
    return sympy.Symbol(*args, **kwargs)


def prefix_indent(prefix, textblock, later_prefix=' '):
    """
    Prefix and indent all lines in *textblock*.
    *prefix* is a prefix string
    *later_prefix* is used on all but the first line, if it is a single character
                   it will be repeated to match length of *prefix*
    """
    textblock = textblock.split('\n')
    line = prefix + textblock[0] + '\n'
    if len(later_prefix) == 1:
        later_prefix = ' ' * len(prefix)
    line = line + '\n'.join([later_prefix + x for x in textblock[1:]])
    if line[-1] != '\n':
        return line + '\n'
    else:
        return line


def transform_multidim_to_1d_decl(decl):
    """
    Transforms ast of multidimensional declaration to a single dimension declaration.
    In-place operation!
    Returns name and dimensions of array (to be used with transform_multidim_to_1d_ref())
    """
    dims = []
    t = decl.type
    while type(t) is c_ast.ArrayDecl:
        dims.append(t.dim)
        t = t.type

    if dims:
        # Multidimensional array
        decl.type.dim = reduce(lambda l, r: c_ast.BinaryOp('*', l, r), dims)
        decl.type.type = t

    return decl.name, dims


def transform_multidim_to_1d_ref(aref, dimension_dict):
    """
    Transform ast of multidimensional reference to a single dimension reference.
    In-place operation!
    """
    dims = []
    name = aref
    while type(name) is c_ast.ArrayRef:
        dims.append(name.subscript)
        name = name.name

    subscript_list = []
    for i, d in enumerate(dims):
        if i == 0:
            subscript_list.append(d)
        else:
            subscript_list.append(c_ast.BinaryOp('*', d, reduce(
                lambda l, r: c_ast.BinaryOp('*', l, r),
                dimension_dict[name.name][-1:-i - 1:-1])))

    aref.subscript = reduce(
        lambda l, r: c_ast.BinaryOp('+', l, r), subscript_list)
    aref.name = name


def transform_array_decl_to_malloc(decl):
    """Transform ast of "type var_name[N]" to "type* var_name = __mm_malloc(N, 32)" (in-place)."""
    if type(decl.type) is not c_ast.ArrayDecl:
        # Not an array declaration, can be ignored
        return

    type_ = c_ast.PtrDecl([], decl.type.type)
    decl.init = c_ast.FuncCall(
        c_ast.ID('aligned_malloc'),
        c_ast.ExprList([
            c_ast.BinaryOp(
                '*',
                c_ast.UnaryOp(
                    'sizeof',
                    c_ast.Typename(None, [], c_ast.TypeDecl(
                        None, [], decl.type.type.type))),
                decl.type.dim),
            c_ast.Constant('size_t', '32')]))
    decl.type = type_


def find_array_references(ast):
    """Return list of array references in AST."""
    if type(ast) is c_ast.ArrayRef:
        return [ast]
    elif type(ast) is list:
        return list(map(find_array_references, ast))
    elif ast is None:
        return []
    else:
        return reduce(operator.add, [find_array_references(o[1]) for o in ast.children()], [])


def force_iterable(f):
    """Will make any functions return an iterable objects by wrapping its result in a list."""
    def wrapper(*args, **kwargs):
        r = f(*args, **kwargs)
        if hasattr(r, '__iter__'):
            return r
        else:
            return [r]
    return wrapper

class KernelBench(KernelCode):
    """
    Kernel information gathered from code using pycparser
    This version allows compilation and generation of code for benchmarking
    """

    def __init__(self, kernel_code, machine, block_factor=None,
                 flop=None, filename=None, initwithrand=False):
        super(KernelBench, self).__init__(kernel_code=kernel_code, machine=machine, filename=filename)
        
        self.block_factor = block_factor
        self.flop = flop
        self.initwithrand = initwithrand
        # need to refer to local lextab, otherwise the systemwide lextab would
        # be imported

    def as_code(self, type_='likwid', from_cli=True):
        """
        generates compilable source code from AST
        *type* can be iaca or likwid.
        """
        assert self.kernel_ast is not None, "AST does not exist, this could be due to running of " \
                                            "kernel description rather than code."

        ast = deepcopy(self.kernel_ast)

        declarations = [d for d in ast.block_items if type(d) is c_ast.Decl]

        # transform multi-dimensional declarations to one dimensional
        # references
        array_dimensions = dict(
            list(map(transform_multidim_to_1d_decl, declarations)))
        # transform to pointer and malloc notation (stack can be too small)
        list(map(transform_array_decl_to_malloc, declarations))


        # if(argc != ) {
        #     printf("Usage: %s N NB M MB L LB repetitions\n", argv[0]);
        #     return 0;
        # }

        sizes_decls_typenames = []
        # add declarations for constants from the executable command line
        if(from_cli):
            var_list = sorted([k.name for k in self.constants])
            blocking = ''
            num_args = 1
            # add declaration of the block
            if self.block_factor == 1:
                var_list.append('block_factor')
                blocking = 'blocking'
                num_args = num_args + 1
            elif self.block_factor > 1:
                for dim in range(0,3):
                    var_list.append('block_factor_' + var_list[dim])
                    blocking = blocking + 'blocking_' + var_list[dim] + ' '
                    num_args = num_args + 1

            i = 1  # subscript for cli input
            for k in var_list:
                # const long N = atoi(argv[1])
                type_decl = c_ast.TypeDecl(k, [], c_ast.IdentifierType(['long']))
                init = c_ast.FuncCall(
                    c_ast.ID('atoll'),
                    c_ast.ExprList([c_ast.ArrayRef(c_ast.ID('argv'), c_ast.Constant('long', str(i)))]))
                i += 1
                #add the variable to the list of variables
                type_name = c_ast.Typename(None, [], type_decl)
                sizes_decls_typenames.append(type_name)
                decl = c_ast.Decl(k, ['const'], [], [], type_decl, init, None)
                ast.block_items.insert(0, decl)

            #build and add the if statement checking the number of arguments passed on the command line
            mysize = ''.join(['size_' + c + ' ' for c in sorted([k.name for k in self.constants])])
            num_args = num_args + len(self.constants)

            argerror_string = 'Wrong number of arguments. Usage:\\n%s {}{}\\n'.format(mysize, blocking)
            argerror_cond = c_ast.BinaryOp('!=',
                    c_ast.ID('argc'), c_ast.Constant('int', num_args))

            argerror = c_ast.FuncCall(c_ast.ID('printf'),
                c_ast.ExprList([
                    c_ast.Constant('string', '"{}"'.format(argerror_string)), c_ast.ID('argv[0]')]))
            arg_check_if = c_ast.If(cond=argerror_cond, iftrue=c_ast.Compound([argerror, c_ast.FuncCall(
                    c_ast.ID('return'), c_ast.ExprList([c_ast.Constant('int', '0')]))]),
                iffalse=None)
            ast.block_items.insert(0, arg_check_if)

        else:
            # add declarations for constants from value passed to stempel
            for name, value in list(self.constants.items()):
                type_decl = c_ast.TypeDecl(str(name), [], c_ast.IdentifierType(['size_t']))
                decl = c_ast.Decl(
                    name, ['const'], [], [],
                        #type_decl, c_ast.Constant('int', str(value)), None)
                        type_decl, c_ast.Constant('long', str(name.name+'_MAX')), None)
                type_name = c_ast.Typename(None, [], type_decl)
                sizes_decls_typenames.append(type_name)
                ast.block_items.insert(0, decl)


        if type_ == 'likwid':
            # Call likwid_markerInit()
            ast.block_items.insert(
                0, c_ast.Constant('string', 'INSERTMACROINIT'))
            # Call likwid_markerThreadInit()
            # Call likwid_markerClose()
            #ast.block_items.append(c_ast.FuncCall(c_ast.ID('likwid_markerClose'), None))
            ast.block_items.append(c_ast.Constant(
                'string', 'INSERTMACROCLOSE'))

        constants = [d.name for d in declarations if d.name.startswith('c')]
        nconstants = len(constants)


        # inject array initialization
        for d in declarations:

            i = ast.block_items.index(d)

            # Build ast to inject
            if array_dimensions[d.name]:
                # this is an array, we need a for loop to initialize it
                # for(init; cond; next) stmt

                # Init: int i = 0;
                counter_name = 'i'
                while counter_name in array_dimensions:
                    counter_name = chr(ord(counter_name) + 1)

                init = c_ast.DeclList([
                    c_ast.Decl(
                        counter_name, [], [], [], c_ast.TypeDecl(
                            counter_name, [], c_ast.IdentifierType(['long'])),
                        c_ast.Constant('long', '0'), None)], None)

                # Cond: i < ... (... is length of array)
                cond = c_ast.BinaryOp(
                    '<',
                    c_ast.ID(counter_name),
                    reduce(lambda l, r: c_ast.BinaryOp('*', l, r), array_dimensions[d.name]))

                # Next: i++
                next_ = c_ast.UnaryOp('++', c_ast.ID(counter_name))


                if d.name == 'W':
                    # get the number of dimensions by fetching the size of W
                    # w_dims can be 2 (1D array + 1 for constants), 3, or 4
                    w_dims = len(array_dimensions.get(d.name))
                    assert (w_dims < 5 and w_dims >
                            1), "STEMPEL can treat stencils up to 3D"
                    if w_dims == 2:
                        if isinstance(d.init.args.exprs[0].right.right, c_ast.Constant): #(sizeof(double)) * (N * 3))
                            factor = float(d.init.args.exprs[0].right.right.value)
                        else:#W[2][N]
                            factor = float(d.init.args.exprs[0].right.left.value)

                    elif w_dims == 3:
                        if isinstance(d.init.args.exprs[0].right.right, c_ast.Constant): #(sizeof(double)) * (M * N) * 3
                            factor = float(d.init.args.exprs[0].right.right.value)
                        elif isinstance(d.init.args.exprs[0].right.left.right, c_ast.Constant):#(sizeof(double)) * (M * 3) * N
                            factor = float(d.init.args.exprs[0].right.left.right.value)
                        else:#(sizeof(double)) * (3 * M) * N
                            factor = float(d.init.args.exprs[0].right.left.left.value)
                    else:  # w_dims == 4
                        if isinstance(d.init.args.exprs[0].right.right, c_ast.Constant): #(sizeof(double)) * (((M * N) * P) * 3)
                            factor = float(d.init.args.exprs[0].right.right.value)
                        elif isinstance(d.init.args.exprs[0].right.left.right, c_ast.Constant):#(sizeof(double)) * (((M * N) * 3) * P)
                            factor = float(d.init.args.exprs[0].right.left.right.value)
                        elif isinstance(d.init.args.exprs[0].right.left.left.right, c_ast.Constant):#(sizeof(double)) * (((M * 3) * N) * P)
                            factor = float(d.init.args.exprs[0].right.left.left.right.value)
                        else:#(sizeof(double)) * (((3 * M) * N) * P)
                            factor = float(d.init.args.exprs[0].right.left.left.left.value)

                    if self.initwithrand is True:
                        # we add 2 in order to get a factor that stabilizes the results
                        # empirical value
                        factor = factor + 2.0
                        rand_max = c_ast.BinaryOp(
                            '/', c_ast.BinaryOp('/', c_ast.FuncCall(
                                c_ast.ID('rand'), c_ast.ExprList([])),
                                                c_ast.Cast(
                                                    c_ast.IdentifierType(['double']),
                                                    c_ast.ID('RAND_MAX'))),
                            c_ast.Constant('float', factor))
                    else:
                        rand_max = c_ast.Constant('float', random.uniform(-23.42,+23.42))
                    #, rand_max))
                    #,
                    # c_ast.Constant('float', factor))
                else:
                    #array a or b
                    if self.initwithrand is True:
                        rand_max = c_ast.BinaryOp(
                                '/', c_ast.FuncCall(
                                    c_ast.ID('rand'),
                                    c_ast.ExprList([])),
                                c_ast.Cast(
                                    c_ast.IdentifierType(['double']),
                                    c_ast.ID('RAND_MAX')))
                    else:
                        rand_max = c_ast.Constant('float', random.uniform(-23.42,+23.42))

                stmt = c_ast.Assignment(
                    '=',
                    c_ast.ArrayRef(c_ast.ID(d.name), c_ast.ID(counter_name)),
                    rand_max)
                pragma_numa = c_ast.Pragma('omp parallel for schedule(runtime)')
                ast.block_items.insert(
                    i + 1, pragma_numa)
                # ast.block_items.insert(
                #     i + 1, c_ast.Pragma('omp parallel for schedule(runtime)')
                ast.block_items.insert(
                    i + 2, c_ast.For(init, cond, next_, stmt))
            else:
                # this is a scalar, so a simple Assignment is enough
                #calculate the factor
                if self.initwithrand is True:
                    factor = 2.0 + float(nconstants)
                    ast.block_items.insert(
                        i + 1, c_ast.Assignment('=', c_ast.ID(d.name),
                            c_ast.BinaryOp(
                            '/', c_ast.BinaryOp(
                                '/', c_ast.FuncCall(
                                    c_ast.ID('rand'),
                                    c_ast.ExprList([])),
                                c_ast.Cast(
                                    c_ast.IdentifierType(['double']),
                                    c_ast.ID('RAND_MAX'))),
                            c_ast.Constant('float', factor))))
                else:
                    ast.block_items.insert(
                        i + 1,
                        c_ast.Assignment(
                            '=', c_ast.ID(d.name),
                            c_ast.Constant('float', random.uniform(-23.42,+23.42))))

        # transform multi-dimensional array references to one dimensional
        # references
        list(map(lambda aref: transform_multidim_to_1d_ref(aref, array_dimensions),
                 find_array_references(ast)))

        dummylist=[]
        # Make sure nothing gets removed by inserting dummy calls
        for d in declarations:
            if array_dimensions[d.name]:
                dummylist.append(c_ast.FuncCall(
                    c_ast.ID('dummy'),
                    c_ast.ExprList([c_ast.ID(d.name)])))
            else:
                dummylist.append(c_ast.FuncCall(
                            c_ast.ID('dummy'),
                            c_ast.ExprList([c_ast.UnaryOp('&', c_ast.ID(d.name))])))

        dummies = c_ast.If(
            cond=c_ast.ID('var_false'),
            iftrue=c_ast.Compound(dummylist),
            iffalse=None)
        ast.block_items.insert(-2,dummies)

        # if we do not want the version accepting inputs from command line,
        # we need to declare the blocking factor
        if not from_cli:
            # add declaration of the block
            if self.block_factor == 1:
                type_decl = c_ast.TypeDecl(
                    'block_factor', [], c_ast.IdentifierType(['size_t']))
                decl = c_ast.Decl(
                    'block_factor', ['const'], [], [],
                    type_decl, c_ast.Constant('size_t', str(self.block_factor)), None)
                ast.block_items.insert(-3, decl)

                # add it to the list of declarations, so it gets passed to the
                # kernel_loop
                declarations.append(decl)
            elif self.block_factor > 1:
                # add declaration of the block
                for dim in range(0,3):
                    name = 'block_factor_' + var_list[dim]
                    type_decl = c_ast.TypeDecl(
                        name, [], c_ast.IdentifierType(['size_t']))
                    decl = c_ast.Decl(
                        name, ['const'], [], [],
                        type_decl, c_ast.Constant('size_t', str(self.block_factor)), None)
                    ast.block_items.insert(-3, decl)

                    # add it to the list of declarations, so it gets passed to the
                    # kernel_loop
                    declarations.append(decl)

        # Wrap everything in a loop
        # size_t repeat = atoi(argv[2])
        type_decl = c_ast.TypeDecl(
            'repeat', [], c_ast.IdentifierType(['size_t']))
        # init = c_ast.FuncCall(
        #     c_ast.ID('atoi'),
        #     c_ast.ExprList([c_ast.ArrayRef(
        #         c_ast.ID('argv'), c_ast.Constant('int', str(len(self.constants)+2)))]))
        # ast.block_items.insert(-3, c_ast.Decl(
        #     'repeat', ['const'], [], [],
        #     type_decl, init, None))
        ast.block_items.insert(-3, c_ast.Decl(
            'repeat', ['const'], [], [],
            type_decl, c_ast.Constant('size_t', '1'), None))

        # timing variables declaration and initialisation
        type_decl = c_ast.TypeDecl(
            'runtime', [], c_ast.IdentifierType(['double']))
        ast.block_items.insert(-3, c_ast.Decl(
            'runtime', ['const'], [], [],
            type_decl, c_ast.Constant('double', '0.0'), None))

        decl = c_ast.Decl('wct_start', [], [], [], c_ast.TypeDecl(
            'wct_start', [], c_ast.IdentifierType(['double'])
        ), None, None)
        ast.block_items.insert(-3, decl)
        decl = c_ast.Decl('wct_end', [], [], [], c_ast.TypeDecl(
            'wct_end', [], c_ast.IdentifierType(['double'])
        ), None, None)
        ast.block_items.insert(-3, decl)
        decl = c_ast.Decl('cput_start', [], [], [], c_ast.TypeDecl(
            'cput_start', [], c_ast.IdentifierType(['double'])
        ), None, None)
        ast.block_items.insert(-3, decl)
        decl = c_ast.Decl('cput_end', [], [], [], c_ast.TypeDecl(
            'cput_end', [], c_ast.IdentifierType(['double'])
        ), None, None)
        ast.block_items.insert(-3, decl)

        # call the timing function at the beginning
        start_timing = c_ast.FuncCall(c_ast.ID('timing'),
                                      c_ast.ExprList([c_ast.UnaryOp('&', c_ast.ID('wct_start')),
                                                      c_ast.UnaryOp('&', c_ast.ID('cput_start'))]))

        # take out the for loop that will be written in a function on top
        forloop = ast.block_items.pop(-2)

        # creating a list of pointer to all the variables of type pointer
        pointers_list = [c_ast.Typename(None, [], c_ast.PtrDecl(
            [], c_ast.TypeDecl(d.name, [], d.type.type))) for d in declarations if type(d.type) is c_ast.PtrDecl]
        first_array_name = pointers_list[0].type.type.declname
        # get the number of dimensions by fetching the size of the first
        # array
        mydims = len(array_dimensions.get(first_array_name))

        # for(n = 0; n < repeat; n++) {...}
        index_name = 'n'
        init = c_ast.DeclList([
            c_ast.Decl(
                index_name, [], [], [], c_ast.TypeDecl(
                    index_name, [], c_ast.IdentifierType(['size_t'])),
                c_ast.Constant('size_t', '0'),
                None)], None)
        cond = c_ast.BinaryOp('<', c_ast.ID(
            index_name), c_ast.ID('repeat'))
        next_ = c_ast.UnaryOp('++', c_ast.ID(index_name))

        expr_list = [c_ast.ID(d.name) for d in declarations] + [c_ast.ID(s) for s in sorted([k.name for k in self.constants])]
        if self.block_factor == 1:
            expr_list = expr_list + [c_ast.ID('block_factor')]
        elif self.block_factor > 1:
            expr_list = expr_list + [c_ast.ID('block_factor_'+var_list[dim]) for dim in range(0,3)]

        stmt = c_ast.FuncCall(c_ast.ID('kernel_loop'),
                              c_ast.ExprList(expr_list))

        swap_tmp = c_ast.Assignment('=', c_ast.ID('tmp'),
                                    c_ast.ID(pointers_list[0].type.type.declname))
        swap_grid = c_ast.Assignment('=', c_ast.ID(pointers_list[0].type.type.declname),
                                     c_ast.ID(pointers_list[1].type.type.declname))
        last_swap = c_ast.Assignment('=', c_ast.ID(pointers_list[1].type.type.declname),
                                     c_ast.ID('tmp'))
        stmt = c_ast.Compound([stmt, swap_tmp, swap_grid, last_swap, dummies] )
        myfor = c_ast.For(init, cond, next_, stmt)

        # call the timing function at the beginning
        end_timing = c_ast.FuncCall(c_ast.ID('timing'),
                                    c_ast.ExprList([c_ast.UnaryOp('&', c_ast.ID('wct_end')),
                                                    c_ast.UnaryOp('&', c_ast.ID('cput_end'))]))

        update_runtime = c_ast.Assignment('=', c_ast.ID('runtime'),
                                          c_ast.BinaryOp('-', c_ast.ID('wct_end'), c_ast.ID('wct_start')))

        update_iter = c_ast.Assignment('*=', c_ast.ID('repeat'),
                                       c_ast.Constant('size_t', '2'))

        # while(runtime<2. || repeat<=2) {...}
        cond = c_ast.BinaryOp( '||',
                c_ast.BinaryOp('<', c_ast.ID( 'runtime'), c_ast.Constant('double', '2.0')),
                c_ast.BinaryOp('<=', c_ast.ID( 'repeat'), c_ast.Constant('size_t', '2')));
        stmt = c_ast.Compound(
            [start_timing, myfor, end_timing, update_runtime, update_iter])

        ast.block_items.insert(-1, c_ast.While(cond, stmt))

        # the variable repeat must be divided by 2 since in the last loop
        # was doubled before exiting
        ast.block_items.insert(-1, c_ast.Assignment('/=',
                                                    c_ast.ID('repeat'), c_ast.Constant('size_t', '2')))

        if type_ == 'likwid':
            # Instrument the outer for-loop with likwid
            # ast.block_items.insert(-2, c_ast.FuncCall(
            #     c_ast.ID('likwid_markerStartRegion'),
            #     c_ast.ExprList([c_ast.Constant('string', '"Sweep"')])))
            ast.block_items.insert(-1,
                                   c_ast.Constant('string', 'INSERTMACROSTART'))

        run_index_name = 'n'
        run_init = c_ast.DeclList([
            c_ast.Decl(
                run_index_name, [], [], [], c_ast.TypeDecl(
                    run_index_name, [], c_ast.IdentifierType(['size_t'])),
                c_ast.Constant('size_t', '0'),
                None)], None)
        run_cond = c_ast.BinaryOp('<', c_ast.ID(run_index_name), c_ast.ID('repeat'))
        run_next = c_ast.UnaryOp('++', c_ast.ID(run_index_name))

        run_expr_list = [c_ast.ID(d.name) for d in declarations] + [c_ast.ID(s) for s in sorted([k.name for k in self.constants])]
        if self.block_factor == 1:
            run_expr_list = run_expr_list + [c_ast.ID('block_factor')]
        elif self.block_factor > 1:
            run_expr_list = run_expr_list + [c_ast.ID('block_factor_'+var_list[dim]) for dim in range(0,3)]

        run_stmt = c_ast.FuncCall(c_ast.ID('kernel_loop'),
                                 c_ast.ExprList(run_expr_list))

        # creating a list of pointer to all the variables of type pointer
        run_pointers_list = [c_ast.Typename(None, [], c_ast.PtrDecl(
            [], c_ast.TypeDecl(d.name, [], d.type.type))) for d in declarations if type(d.type) is c_ast.PtrDecl]

        run_swap_tmp = c_ast.Assignment('=', c_ast.ID('tmp'),
                                    c_ast.ID(run_pointers_list[0].type.type.declname))
        run_swap_grid = c_ast.Assignment('=', c_ast.ID(run_pointers_list[0].type.type.declname),
                                     c_ast.ID(run_pointers_list[1].type.type.declname))
        run_last_swap = c_ast.Assignment('=', c_ast.ID(run_pointers_list[1].type.type.declname),
                                     c_ast.ID('tmp'))
        run_stmt = c_ast.Compound([run_stmt, run_swap_tmp, run_swap_grid, run_last_swap, dummies] )
        run_myfor = c_ast.For(run_init, run_cond, run_next, run_stmt)
        ast.block_items.insert(-1, run_myfor)

        if type_ == 'likwid':
            # close the region "Sweep" of likwid
            ast.block_items.insert(-1,
                                   c_ast.Constant('string', 'INSERTMACROSTOP'))

        # calculate the size of the grid, taking the letters representing
        # its dimensions from the array of constants
        size = '(' + ' * '.join(s for s in sorted([k.name for k in self.constants])) + ')'

        decl = c_ast.Decl('tmp', [], [], [], c_ast.PtrDecl(
            [], c_ast.TypeDecl('tmp', [],
                               pointers_list[0].type.type.type.type)),
                          None, None)
        ast.block_items.insert(-7, decl)

        # creating a list of standard types for all the non-pointer
        # variables
        variables_list = [c_ast.Typename(None, [], c_ast.TypeDecl(
            d.name, [], d.type.type)) for d in declarations if type(d.type) is c_ast.TypeDecl]
        variables_list = variables_list + sizes_decls_typenames

        norm_loop = deepcopy(forloop)
        # generate the LUP expression according to the number of dimensions
        # it is necessary to do so since we do not know a priori how many nested for we have
        # additionally builds the norm of the array
        if mydims == 1:

            if isinstance(forloop.stmt, c_ast.Compound):
                norm_cond_lvalue = norm_loop.stmt.block_items[0].lvalue
                norm_cond_rvalue = norm_loop.stmt.block_items[0].rvalue
            else:
                norm_cond_lvalue = norm_loop.stmt.lvalue
                norm_cond_rvalue = norm_loop.stmt.rvalue

            #extract the end condition of the for and calculate the total number of points on which we iterated
            myleft_lv1 = forloop.cond.right.left
            myright_lv1 = int(forloop.cond.right.right.value)
            myinit_lv1 = int(forloop.init.decls[0].init.value)
            myright_lv1 += myinit_lv1
            mysize_lv1 = c_ast.BinaryOp('-', myleft_lv1, c_ast.Constant('int', myright_lv1))

            lup_expression_long = c_ast.Cast(
                c_ast.IdentifierType(['long long']), mysize_lv1)  # c_ast.ExprList([forloop.cond.right])
            lup_expression = c_ast.Cast(
                c_ast.IdentifierType(['double']), mysize_lv1)

            # # norm
            # point = norm_cond_lvalue
            # # set the name of the grid to the first (the order changed
            # # after the swap)
            # point.name = c_ast.ID(pointers_list[0].type.type.declname)
            # norm_cond_lvalue = c_ast.ID('total')
            # newop = c_ast.BinaryOp('+',
            #                        c_ast.ID('total'),
            #                        c_ast.BinaryOp('*', point, point))
            # norm_cond_rvalue = newop
            # calculate difference between a and b
            point = norm_cond_lvalue
            # set the name of the grid to the first (the order changed
            # after the swap)
            point1 = deepcopy(point)
            point.name = c_ast.ID(pointers_list[0].type.type.declname)
            norm_cond_lvalue = c_ast.ID('total')
            newop = c_ast.BinaryOp('+',
                                   c_ast.ID('total'),
                                   c_ast.BinaryOp('-', point, point1))
            norm_cond_rvalue = newop


            if isinstance(forloop.stmt, c_ast.Compound):
                norm_loop.stmt.block_items[0].lvalue = norm_cond_lvalue
                norm_loop.stmt.block_items[0].rvalue = norm_cond_rvalue
            else:
                norm_loop.stmt.lvalue = norm_cond_lvalue
                norm_loop.stmt.rvalue = norm_cond_rvalue

        elif mydims == 2:

            # mycode = CGenerator().visit(forloop.stmt.block_items[0].cond.right)
            # print(mycode)
            # exit(1)

            if isinstance(forloop.stmt, c_ast.Compound):
                right_cond = forloop.stmt.block_items[0].cond.right
                lv2_init = forloop.stmt.block_items[0].init
                norm_cond_lvalue = norm_loop.stmt.block_items[0].stmt.block_items[0].lvalue
                norm_cond_rvalue = norm_loop.stmt.block_items[0].stmt.block_items[0].rvalue

            else: #no compound
                right_cond = forloop.stmt.cond.right
                lv2_init = forloop.stmt.init
                norm_cond_lvalue = norm_loop.stmt.stmt.lvalue
                norm_cond_rvalue = norm_loop.stmt.stmt.lvalue

            #extract the end condition of the for and calculate the total number of points on which we iterated
            myleft_lv1 = forloop.cond.right.left
            myright_lv1 = int(forloop.cond.right.right.value)
            myinit_lv1 = int(forloop.init.decls[0].init.value)
            myright_lv1 += myinit_lv1
            mysize_lv1 = c_ast.BinaryOp('-', myleft_lv1, c_ast.Constant('int', myright_lv1))
            #extract the end condition of the for and calculate the total number of points on which we iterated
            myleft_lv2 = right_cond.left
            myright_lv2 = int(right_cond.right.value)
            myinit_lv2 = int(lv2_init.decls[0].init.value)
            myright_lv2 += myinit_lv2
            mysize_lv2 = c_ast.BinaryOp('-', myleft_lv2, c_ast.Constant('int', myright_lv2))

            lup_expression_long = c_ast.BinaryOp(
                '*', c_ast.Cast(
                    c_ast.IdentifierType(
                        ['long long']), mysize_lv1),
                c_ast.Cast(c_ast.IdentifierType(
                    ['long long']),
                mysize_lv2))
            lup_expression = c_ast.BinaryOp(
                '*', c_ast.Cast(
                    c_ast.IdentifierType(
                        ['double']), mysize_lv1),
                c_ast.Cast(c_ast.IdentifierType(
                    ['double']),
                mysize_lv2))

            # # norm
            # point = norm_cond_lvalue
            # # set the name of the grid to the first (the order changed
            # # after the swap)
            # point.name = c_ast.ID(pointers_list[0].type.type.declname)
            # norm_cond_lvalue = c_ast.ID('total')
            # newop = c_ast.BinaryOp('+',
            #                        c_ast.ID('total'),
            #                        c_ast.BinaryOp('*', point, point))
            # norm_cond_rvalue = newop
            # calculate difference between a and b
            point = norm_cond_lvalue
            # set the name of the grid to the first (the order changed
            # after the swap)
            point1 = deepcopy(point)
            point.name = c_ast.ID(pointers_list[0].type.type.declname)
            norm_cond_lvalue = c_ast.ID('total')
            newop = c_ast.BinaryOp('+',
                                   c_ast.ID('total'),
                                   c_ast.BinaryOp('-', point, point1))
            norm_cond_rvalue = newop

            #rebuild the norm loop
            if isinstance(norm_loop.stmt, c_ast.Compound):
                norm_loop.stmt.block_items[0].stmt.block_items[0].lvalue = norm_cond_lvalue
                norm_loop.stmt.block_items[0].stmt.block_items[0].rvalue = norm_cond_rvalue

            else: #no compound
                norm_loop.stmt.stmt.lvalue = norm_cond_lvalue
                norm_loop.stmt.stmt.lvalue = norm_cond_rvalue

        elif mydims == 3:

            if isinstance(forloop.stmt, c_ast.Compound):
                right_cond = forloop.stmt.block_items[0].cond.right
                lv2_init = forloop.stmt.block_items[0].init
                norm_cond_lvalue = norm_loop.stmt.block_items[0].stmt.block_items[0].stmt.block_items[0].lvalue
                norm_cond_rvalue = norm_loop.stmt.block_items[0].stmt.block_items[0].stmt.block_items[0].rvalue
                stmt_rcond = forloop.stmt.block_items[0].stmt.block_items[0].cond.right
                lv3_init = forloop.stmt.block_items[0].stmt.block_items[0].init
            else: #no compound
                right_cond = forloop.stmt.cond.right
                lv2_init = forloop.stmt.init
                norm_cond_lvalue = norm_loop.stmt.stmt.stmt.lvalue
                norm_cond_rvalue = norm_loop.stmt.stmt.stmt.lvalue
                stmt_rcond = forloop.stmt.stmt.cond.right
                lv3_init = forloop.stmt.stmt.init

            #extract the end condition of the for and calculate the total number of points on which we iterated
            myleft_lv1 = forloop.cond.right.left
            myright_lv1 = int(forloop.cond.right.right.value)
            myinit_lv1 = int(forloop.init.decls[0].init.value)
            myright_lv1 += myinit_lv1
            mysize_lv1 = c_ast.BinaryOp('-', myleft_lv1, c_ast.Constant('int', myright_lv1))
            #extract the end condition of the for and calculate the total number of points on which we iterated
            myleft_lv2 = right_cond.left
            myright_lv2 = int(right_cond.right.value)
            myinit_lv2 = int(lv2_init.decls[0].init.value)
            myright_lv2 += myinit_lv2
            mysize_lv2 = c_ast.BinaryOp('-', myleft_lv2, c_ast.Constant('int', myright_lv2))
            #extract the end condition of the for and calculate the total number of points on which we iterated
            myleft_lv3 = stmt_rcond.left
            myright_lv3 = int(stmt_rcond.right.value)
            myinit_lv3 = int(lv3_init.decls[0].init.value)
            myright_lv3 += myinit_lv3
            mysize_lv3 = c_ast.BinaryOp('-', myleft_lv3, c_ast.Constant('int', myright_lv3))

            lup_expression_long = c_ast.BinaryOp(
                '*', c_ast.BinaryOp(
                    '*', c_ast.Cast(c_ast.IdentifierType(
                        ['long long']), mysize_lv1),
                    c_ast.Cast(c_ast.IdentifierType(
                        ['long long']), mysize_lv2)),
                c_ast.Cast(c_ast.IdentifierType(
                    ['long long']), mysize_lv3))
            lup_expression = c_ast.BinaryOp(
                '*', c_ast.BinaryOp(
                    '*', c_ast.Cast(c_ast.IdentifierType(
                        ['double']), mysize_lv1),
                    c_ast.Cast(c_ast.IdentifierType(
                        ['double']), mysize_lv2)),
                c_ast.Cast(c_ast.IdentifierType(
                    ['double']), mysize_lv3))


            # # calculate norm
            # point = norm_cond_lvalue
            # # set the name of the grid to the first (the order changed
            # # after the swap)
            # point.name = c_ast.ID(pointers_list[0].type.type.declname)
            # norm_cond_lvalue = c_ast.ID('total')
            # newop = c_ast.BinaryOp('+',
            #                        c_ast.ID('total'),
            #                        c_ast.BinaryOp('*', point, point))
            # calculate difference between a and b
            point = norm_cond_lvalue
            # set the name of the grid to the first (the order changed
            # after the swap)
            point1 = deepcopy(point)
            point.name = c_ast.ID(pointers_list[0].type.type.declname)
            norm_cond_lvalue = c_ast.ID('total')
            newop = c_ast.BinaryOp('+',
                                   c_ast.ID('total'),
                                   c_ast.BinaryOp('-', point, point1))
            norm_cond_rvalue = newop


            if isinstance(norm_loop.stmt, c_ast.Compound):
                norm_loop.stmt.block_items[0].stmt.block_items[0].stmt.block_items[0].lvalue = norm_cond_lvalue
                norm_loop.stmt.block_items[0].stmt.block_items[0].stmt.block_items[0].rvalue = norm_cond_rvalue
            else: #no compound
                norm_loop.stmt.stmt.stmt.lvalue = norm_cond_lvalue
                norm_loop.stmt.stmt.stmt.lvalue = norm_cond_rvalue

        # we build mlup. should be like:
        # (double)iter*(size_x-ghost)*(size_y-ghost)*(size_z-ghost)/runtime/1000000.
        lup_expression_long = c_ast.BinaryOp('*',
                                        c_ast.Cast(
                                            c_ast.IdentifierType(['long long']), c_ast.ID('repeat')),
                                        lup_expression_long)
        lup_expression = c_ast.BinaryOp('*',
                                        c_ast.Cast(
                                            c_ast.IdentifierType(['double']), c_ast.ID('repeat')),
                                        lup_expression)
        # cast it to double since the first variables are ints
        #LUP_expr_cast =  c_ast.Cast(c_ast.IdentifierType(['double']), lup_expression)
        # we put all together to get mlup
        #MLUP = c_ast.BinaryOp('/', LUP_expr_cast, c_ast.BinaryOp('*', c_ast.ID('runtime'), c_ast.Constant('double', '1000000.')))
        lup = c_ast.BinaryOp('/', lup_expression, c_ast.ID('runtime'))
        glup = c_ast.BinaryOp('/', lup, c_ast.Constant('double', '1000000000.'))

        # insert the printf of the stats
        mystring = "iterations: %d\\n"
        ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                                                  c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)),
                                                                  c_ast.ID('repeat')])))
        mystring = "Total iterations: %lld LUP\\n"
        ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                                                  c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)),
                                                                  lup_expression_long])))

        if self.flop:
            #flops calculated before building the benchmark code
            flop_per_lup = "FLOP: {}\\n".format(self.flop)
            ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(flop_per_lup))])))

            #total work FLOP*LUPs
            mystring = "Total work: %lld FLOP\\n"
            total_work_long = c_ast.BinaryOp('*', lup_expression_long, c_ast.Constant('int', self.flop))
            total_work = c_ast.BinaryOp('*', lup_expression, c_ast.Constant('int', self.flop))
            ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)), total_work_long])))

            #GLUPs = lup_expr / runtime / 1e9
            mystring = "performance: %lf GLUP/s\\n"
            ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)), glup])))

            #GFLOPs = lup_expr / runtime / 1e9
            mystring = "Performance in GFLOP/s: %lf\\n"
            flops = c_ast.BinaryOp('/', total_work, c_ast.ID('runtime'))
            gflops = c_ast.BinaryOp('/', flops, c_ast.Constant('double', '1000000000.'))
            ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)), gflops])))


        else:
            mystring = "Performance in GLUP/s: %lf\\n"
            ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                                                  c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)),
                                                                  glup])))
        mystring = "size: %d\\n"
        ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                                                  c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)),
                                                                  c_ast.ID(size)])))
        mystring = "time: %lf\\n"
        ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                                                  c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)),
                                                                  c_ast.ID('runtime')])))


        # insert the loop computing the total squared
        decl = c_ast.Decl('total', [], [], [], c_ast.TypeDecl(
            'total', [], c_ast.IdentifierType(['double'])
        ), c_ast.Constant('double', '0.0'), None)
        ast.block_items.insert(-1, decl)
        ast.block_items.insert(-1, norm_loop)

        #calculate and print the norm of a
        # insert the printf of the norm
        # mystring = "norm(a): %lf\\n"

        # sqrt_total = c_ast.FuncCall(c_ast.ID('sqrt'),
        #                             c_ast.ExprList([c_ast.ID('total')]))
        # ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
        #                                           c_ast.ExprList([
        #                                               c_ast.Constant('string', '"{}"'.format(mystring)), sqrt_total])))

        # calculate and print a(i,j) - b(i,j)
        mystring = "diff(a-b): %lf\\n"
        ast.block_items.insert(-1, c_ast.FuncCall(c_ast.ID('printf'),
                                                  c_ast.ExprList([c_ast.Constant('string', '"{}"'.format(mystring)),
                                                                  c_ast.ID('total')])))

        # embed compound into main FuncDecl
        decl = c_ast.Decl('main', [], [], [], c_ast.FuncDecl(
            c_ast.ParamList([
                c_ast.Typename(
                    None, [], c_ast.TypeDecl(
                        'argc', [], c_ast.IdentifierType(['int']))),
                c_ast.Typename(None, [], c_ast.PtrDecl(
                    [], c_ast.PtrDecl([], c_ast.TypeDecl(
                        'argv', [], c_ast.IdentifierType(
                            ['char'])))))]),
            c_ast.TypeDecl('main', [], c_ast.IdentifierType(['int']))),
                          None, None)

        ast = c_ast.FuncDef(decl, None, ast)


        myvariables = []
        for i in range(0, mydims):
            myvariables.append(chr(ord('i') + i))

        pragma_int = None
        if self.block_factor:
            pragma_int = c_ast.Pragma('omp for schedule(runtime)')
        else:
            pragma_int = c_ast.Pragma('omp parallel for schedule(runtime)')

        # declaring the function of the kernel with the parameters list built
        # before
        #declare and use in the main file
        decl_sig_kernel = c_ast.Decl('kernel_loop', [], [], [], c_ast.FuncDecl(
            c_ast.ParamList(pointers_list + variables_list),
            c_ast.TypeDecl('kernel_loop', [], c_ast.IdentifierType(['extern void']))),
                          None, None)

        #declare to use in kernel.c
        decl_kernel = c_ast.Decl('kernel_loop', [], [], [], c_ast.FuncDecl(
            c_ast.ParamList(pointers_list + variables_list),
            c_ast.TypeDecl('kernel_loop', [], c_ast.IdentifierType(['void']))),
                          None, None)

        mycompound = None
        if self.block_factor and mydims > 1:

            if isinstance(forloop.stmt, c_ast.Compound):
                    myblockstmt = forloop.stmt.block_items[0]
            else:
                myblockstmt = forloop.stmt

            if mydims == 2:  # blocking on the inner-most loop
                beginning = myvariables[0] + 'b'
                end = myvariables[0] + 'end'

                pragma = c_ast.Pragma(
                    'omp parallel')# private({}, {})'.format(beginning, end))

                #mycode = CGenerator().visit(forloop.stmt.block_items[0].init.decls[0])

                init = c_ast.DeclList([
                    c_ast.Decl(
                        beginning, [], [], [], c_ast.TypeDecl(
                            beginning, [], c_ast.IdentifierType(['size_t'])),
                        myblockstmt.init.decls[0].init, None)], None)
                # for(jb = 1; jb < N-1; jb+=block_factor) {...}reduce(lambda l,
                # r: c_ast.BinaryOp('*', l, r), array_dimensions[d.name]))
                cond = c_ast.BinaryOp('<', c_ast.ID(
                    beginning), myblockstmt.cond.right)
                next_ = c_ast.BinaryOp(
                    '+=', c_ast.ID(beginning), c_ast.ID('block_factor'))

                decl = c_ast.Decl(end, [], [], [], c_ast.TypeDecl(
                    end, [], c_ast.IdentifierType(['size_t'])), c_ast.FuncCall(
                    c_ast.ID('min'), c_ast.ExprList([
                        c_ast.BinaryOp(
                                '+', c_ast.ID(beginning), c_ast.ID('block_factor')),
                        myblockstmt.cond.right])), None)

                myblockstmt.init.decls[0].init = c_ast.ID(beginning)
                myblockstmt.cond.right = c_ast.ID(end)

                mycompound = c_ast.Compound(
                    [decl, pragma_int, forloop])

                newfor = c_ast.For(init, cond, next_, mycompound)

                mycompound = c_ast.Compound([pragma, newfor])

            elif mydims == 3:  # blocking on the middle loop
                if self.block_factor == 1:
                    beginning = myvariables[1] + 'b'
                    end = myvariables[1] + 'end'
                    pragma = c_ast.Pragma(
                        'omp parallel')# private({}, {})'.format(beginning, end))

                    init = c_ast.DeclList([
                        c_ast.Decl(
                            beginning, [], [], [], c_ast.TypeDecl(
                                beginning, [], c_ast.IdentifierType(['size_t'])),
                            myblockstmt.init.decls[0].init,
                            None)], None)
                    # for(jb = 1; jb < N-1; jb+=block_factor) {...}reduce(lambda l,
                    # r: c_ast.BinaryOp('*', l, r), array_dimensions[d.name]))
                    cond = c_ast.BinaryOp('<', c_ast.ID(
                        beginning), myblockstmt.cond.right)
                    next_ = c_ast.BinaryOp(
                        '+=', c_ast.ID(beginning), c_ast.ID('block_factor'))

                    decl = c_ast.Decl(end, [], [], [], c_ast.TypeDecl(
                        end, [], c_ast.IdentifierType(['size_t'])), c_ast.FuncCall(
                        c_ast.ID('min'), c_ast.ExprList([
                            c_ast.BinaryOp(
                                    '+', c_ast.ID(beginning), c_ast.ID('block_factor')),
                            myblockstmt.cond.right])), None)

                    myblockstmt.init.decls[0].init = c_ast.ID(beginning)
                    myblockstmt.cond.right = c_ast.ID(end)

                    mycompound = c_ast.Compound(
                        [decl, pragma_int, forloop])

                    newfor = c_ast.For(init, cond, next_, mycompound)

                    mycompound = c_ast.Compound([pragma, newfor])
                elif self.block_factor > 1:
                    blk_forloop = forloop
                    blk = [blk_forloop.stmt.block_items[0].stmt.block_items[0],
                           blk_forloop.stmt.block_items[0],
                           blk_forloop]

                    from_ = [blk[0].init.decls[0].init,blk[1].init.decls[0].init,blk[2].init.decls[0].init]
                    to_   = [blk[0].cond.right,blk[1].cond.right,blk[2].cond.right]

                    for dim in range(0,mydims):
                        beginning = myvariables[dim] + 'b'
                        end = myvariables[dim] + 'end'

                        blk[dim].init.decls[0].init = c_ast.ID(beginning)
                        blk[dim].cond.right = c_ast.ID(end)

                    mycompound = c_ast.Compound([pragma_int,blk_forloop])

                    for dim in range(0,mydims):
                        beginning = myvariables[dim] + 'b'
                        end = myvariables[dim] + 'end'

                        init = c_ast.DeclList([
                            c_ast.Decl(
                                beginning, [], [], [], c_ast.TypeDecl(
                                    beginning, [], c_ast.IdentifierType(['size_t'])),
                                from_[dim],
                                None)], None)
                        # for(jb = 1; jb < N-1; jb+=block_factor) {...}reduce(lambda l,
                        # r: c_ast.BinaryOp('*', l, r), array_dimensions[d.name]))
                        cond = c_ast.BinaryOp('<', c_ast.ID(
                            beginning), to_[dim])
                        next_ = c_ast.BinaryOp(
                            '+=', c_ast.ID(beginning), c_ast.ID('block_factor_' + var_list[2-dim]))

                        decl = c_ast.Decl(end, [], [], [], c_ast.TypeDecl(
                            end, [], c_ast.IdentifierType(['size_t'])), c_ast.FuncCall(
                            c_ast.ID('min'), c_ast.ExprList([
                                c_ast.BinaryOp(
                                        '+', c_ast.ID(beginning), c_ast.ID('block_factor_' + var_list[2-dim])),
                                to_[dim]])), None)

                        mycompound = c_ast.For(init, cond, next_, c_ast.Compound([decl, mycompound]))

                    pragma = c_ast.Pragma('omp parallel')# private({}, {})'.format(beginning, end))
                    mycompound = c_ast.Compound([pragma,mycompound])
        else:
            mycompound = c_ast.Compound([pragma_int, forloop])

        # logging.warning(type(forloop))
        ast_kernel = c_ast.FuncDef(decl_kernel, None, mycompound)

        # ast = c_ast.FuncDef(decl_sig_kernel, None, ast)

        ast = c_ast.FileAST([ast])

        #add the declaration of the kernel function to the AST
        ast.ext.insert(0, decl_sig_kernel)

        # add dummy function declaration
        decl = c_ast.Decl('dummy', [], [], [], c_ast.FuncDecl(
            c_ast.ParamList([c_ast.Typename(None, [], c_ast.PtrDecl(
                [], c_ast.TypeDecl(None, [], c_ast.IdentifierType(['double']))))]),
            c_ast.TypeDecl('dummy', [], c_ast.IdentifierType(['void']))),
                          None, None)
        ast.ext.insert(0, decl)

        # add alloc function declaration
        decl = c_ast.Decl('aligned_malloc', [], [], [], c_ast.FuncDecl(
            c_ast.ParamList([
                c_ast.Typename(None, [], c_ast.TypeDecl(
                    None, [], c_ast.IdentifierType(['size_t']))),
                c_ast.Typename(None, [], c_ast.TypeDecl(
                    None, [], c_ast.IdentifierType(['size_t'])))
            ]),
            c_ast.TypeDecl('aligned_malloc', [], c_ast.IdentifierType(['void*']))),
                          None, None)
        ast.ext.insert(0, decl)

        # add external var_false declaration
        decl = c_ast.Decl('var_false', [], ['extern'], [], c_ast.TypeDecl(
            'var_false', [], c_ast.IdentifierType(['int'])
        ), None, None)
        ast.ext.insert(1, decl)

        # convert to code string
        code = CGenerator().visit(ast)

        # add empty line on top
        code = '\n' + code

        if not from_cli:
            # add defines of the variables storing the size of the dimensions
            # for name, value in list(self.constants.items()):
            #     line = '#define {} {}L\n'.format(name.name+'_MAX', value)
            #     code = line + code
            #code = '\n' + code
            code = '#include "dim_input.h"\n' + code


        ifdefperf = '#ifdef LIKWID_PERFMON\n'
        endif = '#endif\n'

        code = ifdefperf + '#include <likwid.h>\n' + endif + code

        # add "#include"s for dummy, var_false and stdlib (for malloc)
        code = '#include "kernel.c"\n' + code
        code = '#include "kerncraft.h"\n' + code
        code = '#include "timing.h"\n' + code

        code = '#include <stdlib.h>\n' + code

        # substitute the string added with the macro, since there is no way to
        # add MACROs with pycparser. It is a workaround
        # TODO change the code creation in a way to use MACROs
        pragraomp = '  #pragma omp parallel\n  {}\n    ' + '{}' + '\n  {}'

        likwid_init = 'LIKWID_MARKER_INIT;'
        likwid_register = 'LIKWID_MARKER_REGISTER("Sweep");'
        likwid_register = pragraomp.format('{', likwid_register, '}')
        macroinit = '\n  ' + ifdefperf + '  ' + likwid_init + '\n'
        macroinit += likwid_register + '\n  ' + endif
        code = code.replace('INSERTMACROINIT;', macroinit)

        start_sweep = 'LIKWID_MARKER_START("Sweep");'
        pragma_start_sweep = pragraomp.format('{', start_sweep, '}')
        macrostart = '\n  ' + ifdefperf + pragma_start_sweep + '\n'
        code = code.replace('INSERTMACROSTART;', macrostart)

        stop_sweep = 'LIKWID_MARKER_STOP("Sweep");'
        marker_get = '\n    int nevents = 50, count = 50;\n'
        marker_get += '    double events[50];\n'
        marker_get += '    LIKWID_MARKER_GET("Sweep", &nevents, events, &runtime, &count );'
        pragma_stop_sweep = pragraomp.format('{', stop_sweep + marker_get, '}')
        macrostop = pragma_stop_sweep[2:] + '\n  ' + endif
        code = code.replace('INSERTMACROSTOP;', macrostop)

        likwid_close = 'LIKWID_MARKER_CLOSE;'
        macroclose = '\n  ' + ifdefperf + '  ' + likwid_close + '\n  ' + endif
        code = code.replace('INSERTMACROCLOSE;', macroclose)

        # remove trailing ";": it must be fixed in the place where it is
        # accidentally added, i.e. when adding the function kernel_loop.
        # just a workaround
        code = code.rstrip()
        if code.endswith(';'):
            code = code = code[0:-1]

        kernel = CGenerator().visit(ast_kernel)
        if self.block_factor:
            kernel = '#endif\n' + kernel
            kernel = '#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )\n' + kernel
            kernel = '#ifndef min\n' + kernel

        # return mycode
        return code, kernel

