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
from functools import reduce

#import logging
from distutils.spawn import find_executable

import sympy

from pycparser import CParser, c_ast, plyparser
from pycparser.c_generator import CGenerator

import kerncraft
from kerncraft.kernel import Kernel
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
            c_ast.Constant('int', '32')]))
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

class KernelBench(Kernel):
    """
    Kernel information gathered from code using pycparser
    This version allows compilation and generation of code for benchmarking
    """

    def __init__(self, kernel_code, machine, block_factor=None, flop=None, filename=None):
        super(KernelBench, self).__init__(machine=machine)

        # Initialize state
        self.asm_block = None

        self.kernel_code = kernel_code
        self._filename = filename
        self.block_factor = block_factor
        self.flop = flop
        # need to refer to local lextab, otherwise the systemwide lextab would
        # be imported
        parser = CParser()
        try:
            self.kernel_ast = parser.parse(
                self._as_function(),
                filename=filename).ext[0].body
        except plyparser.ParseError as e:
            print('Error parsing kernel code:', e)
            sys.exit(1)

        self._process_code()

        self.check()

    def print_kernel_code(self, output_file=sys.stdout):
        """Print source code of kernel."""
        print(self.kernel_code, file=output_file)

    def _as_function(self, func_name='test', filename=None):
        if filename is None:
            filename = ''
        else:
            filename = '"{}"'.format(filename)
        return '#line 0 \nvoid {}() {{\n#line 1 {}\n{}\n#line 999 \n}}'.format(
            func_name, filename, self.kernel_code)

    def clear_state(self):
        """Clears mutable internal states"""
        super(KernelBench, self).clear_state()
        self.asm_block = None

    def _process_code(self):
        assert type(self.kernel_ast) is c_ast.Compound, "Kernel has to be a compound statement"
        assert all([type(s) in [c_ast.Decl, c_ast.Pragma]
                    for s in self.kernel_ast.block_items[:-1]]), \
            'all statements before the for loop need to be declarations or pragmas'
        assert type(self.kernel_ast.block_items[-1]) is c_ast.For, \
            'last statement in kernel code must be a loop'

        for item in self.kernel_ast.block_items[:-1]:
            if type(item) is c_ast.Pragma:
                continue
            array = type(item.type) is c_ast.ArrayDecl

            if array:
                dims = []
                t = item.type
                while type(t) is c_ast.ArrayDecl:
                    dims.append(self.conv_ast_to_sym(t.dim))
                    t = t.type

                assert len(t.type.names) == 1, "only single types are supported"
                self.set_variable(item.name, t.type.names[0], tuple(dims))

            else:
                assert len(item.type.type.names) == 1, "only single types are supported"
                self.set_variable(item.name, item.type.type.names[0], None)

        floop = self.kernel_ast.block_items[-1]
        self._p_for(floop)

    def conv_ast_to_sym(self, math_ast):
        """
        Convert mathematical expressions to a sympy representation.
        May only contain paranthesis, addition, subtraction and multiplication from AST.
        """
        if type(math_ast) is c_ast.ID:
            return symbol_pos_int(math_ast.name)
        elif type(math_ast) is c_ast.Constant:
            return sympy.Integer(math_ast.value)
        else:  # elif type(dim) is c_ast.BinaryOp:
            op = {
                '*': operator.mul,
                '+': operator.add,
                '-': operator.sub
            }

            return op[math_ast.op](
                self.conv_ast_to_sym(math_ast.left),
                self.conv_ast_to_sym(math_ast.right))

    def _get_offsets(self, aref, dim=0):
        """
        Return a tuple of offsets of an ArrayRef object in all dimensions.
        The index order is right to left (c-code order).
        e.g. c[i+1][j-2] -> (-2, +1)
        If aref is actually a c_ast.ID, None will be returned.
        """
        if isinstance(aref, c_ast.ID):
            return None

        # Check for restrictions
        assert type(aref.name) in [c_ast.ArrayRef, c_ast.ID], \
            "array references must only be used with variables or other array references"
        assert type(aref.subscript) in [c_ast.ID, c_ast.Constant, c_ast.BinaryOp], \
            'array subscript must only contain variables or binary operations'

        # Convert subscript to sympy and append
        idxs = [self.conv_ast_to_sym(aref.subscript)]

        # Check for more indices (multi-dimensional access)
        if type(aref.name) is c_ast.ArrayRef:
            idxs += self._get_offsets(aref.name, dim=dim+1)

        # Reverse to preserver order (the subscripts in the AST are traversed backwards)
        if dim == 0:
            idxs.reverse()

        return tuple(idxs)

    @classmethod
    def _get_basename(cls, aref):
        """
        Return base name of ArrayRef object.
        e.g. c[i+1][j-2] -> 'c'
        """
        if isinstance(aref.name, c_ast.ArrayRef):
            return cls._get_basename(aref.name)
        elif isinstance(aref.name, str):
            return aref.name
        else:
            return aref.name.name

    def _p_for(self, floop):
        # Check for restrictions
        assert type(floop) is c_ast.For, "May only be a for loop"
        assert hasattr(floop, 'init') and hasattr(floop, 'cond') and hasattr(floop, 'next'), \
            "Loop must have initial, condition and next statements."
        assert type(floop.init) is c_ast.DeclList, \
            "Initialization of loops need to be declarations."
        assert len(floop.init.decls) == 1, "Only single declaration is allowed in init. of loop."
        assert floop.cond.op in '<', "only lt (<) is allowed as loop condition"
        assert type(floop.cond.left) is c_ast.ID, 'left of cond. operand has to be a variable'
        assert type(floop.cond.right) in [c_ast.Constant, c_ast.ID, c_ast.BinaryOp], \
            'right of cond. operand has to be a constant, a variable or a binary operation'
        assert type(floop.next) in [c_ast.UnaryOp, c_ast.Assignment], \
            'next statement has to be a unary or assignment operation'
        assert floop.next.op in ['++', 'p++', '+='], 'only ++ and += next operations are allowed'
        assert type(floop.stmt) in [c_ast.Compound, c_ast.Assignment, c_ast.For], \
            'the inner loop may contain only assignments or compounds of assignments'

        if type(floop.cond.right) is c_ast.ID:
            const_name = floop.cond.right.name
            iter_max = symbol_pos_int(const_name)
        elif type(floop.cond.right) is c_ast.Constant:
            iter_max = sympy.Integer(floop.cond.right.value)
        else:  # type(floop.cond.right) is c_ast.BinaryOp
            bop = floop.cond.right
            assert bop.op in '+-*', ('only addition (+), substraction (-) and multiplications (*) '
                                     'are accepted operators')
            iter_max = self.conv_ast_to_sym(bop)

        iter_min = self.conv_ast_to_sym(floop.init.decls[0].init)

        if type(floop.next) is c_ast.Assignment:
            assert type(floop.next.lvalue) is c_ast.ID, \
                'next operation may only act on loop counter'
            assert type(floop.next.rvalue) is c_ast.Constant, 'only constant increments are allowed'
            assert floop.next.lvalue.name == floop.cond.left.name == floop.init.decls[0].name, \
                'initial, condition and next statement of for loop must act on same loop ' \
                'counter variable'
            step_size = int(floop.next.rvalue.value)
        else:
            assert type(floop.next.expr) is c_ast.ID, 'next operation may only act on loop counter'
            assert floop.next.expr.name == floop.cond.left.name == floop.init.decls[0].name, \
                'initial, condition and next statement of for loop must act on same loop ' \
                'counter variable'
            assert isinstance(floop.next, c_ast.UnaryOp), 'only assignment or unary operations ' \
                'are allowed for next statement of loop.'
            assert floop.next.op in ['++', 'p++', '--', 'p--'], 'Unary operation can only be ++ ' \
                'or -- in next statement'
            if floop.next.op in ['++', 'p++']:
                step_size = sympy.Integer('1')
            else:  # floop.next.op in ['--', 'p--']:
                step_size = sympy.Integer('-1')

        # Document for loop stack
        self._loop_stack.append(
            # (index name, min, max, step size)
            (floop.init.decls[0].name, iter_min, iter_max, step_size)
        )

        # Traverse tree
        if type(floop.stmt) is c_ast.For:
            self._p_for(floop.stmt)
        elif type(floop.stmt) is c_ast.Assignment:
            self._p_assignment(floop.stmt)
        # Handle For if it is the last statement, only preceeded by Pragmas
        elif type(floop.stmt.block_items[-1]) is c_ast.For and \
                all([type(s) == c_ast.Pragma for s in floop.stmt.block_items[:-1]]):
            self._p_for(floop.stmt.block_items[-1])
        else:  # type(floop.stmt) is c_ast.Compound
            # Handle Assignments
            for assgn in floop.stmt.block_items:
                self._p_assignment(assgn)

    def _p_assignment(self, stmt):
        # Check for restrictions
        assert type(stmt) is c_ast.Assignment, \
            "Only assignment statements are allowed in loops."
        assert type(stmt.lvalue) in [c_ast.ArrayRef, c_ast.ID], \
            "Only assignment to array element or varialbe is allowed."

        write_and_read = False
        if stmt.op != '=':
            write_and_read = True
            op = stmt.op.strip('=')
            self._flops[op] = self._flops.get(op, 0)+1

        # Document data destination
        # self.destinations[dest name] = [dest offset, ...])
        self.destinations.setdefault(self._get_basename(stmt.lvalue), set())
        self.destinations[self._get_basename(stmt.lvalue)].add(
            self._get_offsets(stmt.lvalue))

        if write_and_read:
            # this means that +=, -= or something of that sort was used
            self.sources.setdefault(self._get_basename(stmt.lvalue), set())
            self.sources[self._get_basename(stmt.lvalue)].add(
                self._get_offsets(stmt.lvalue))

        # Traverse tree
        self._p_sources(stmt.rvalue)

    def _p_sources(self, stmt):
        sources = []
        assert type(stmt) in \
            [c_ast.ArrayRef, c_ast.Constant, c_ast.ID, c_ast.BinaryOp, c_ast.UnaryOp], \
            'only references to arrays, constants and variables as well as binary operations ' + \
            'are supported'
        assert type(stmt) is not c_ast.UnaryOp or stmt.op in ['-', '--', '++', 'p++', 'p--'], \
            'unary operations are only allowed with -, -- and ++'

        if type(stmt) in [c_ast.ArrayRef, c_ast.ID]:
            # Document data source
            bname = self._get_basename(stmt)
            self.sources.setdefault(bname, set())
            self.sources[bname].add(self._get_offsets(stmt))
        elif type(stmt) is c_ast.BinaryOp:
            # Traverse tree
            self._p_sources(stmt.left)
            self._p_sources(stmt.right)

            self._flops[stmt.op] = self._flops.get(stmt.op, 0)+1
        elif type(stmt) is c_ast.UnaryOp:
            self._p_sources(stmt.expr)
            self._flops[stmt.op] = self._flops.get(stmt.op[-1], 0)+1

        return sources

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
            var_list = [k.name for k in self.constants]
            blocking = ''
            num_args = 1
            # add declaration of the block
            if self.block_factor:
                var_list.append('block_factor')
                blocking = 'blocking'
                num_args = num_args + 1

            i = 1  # subscript for cli input
            for k in var_list:
                # cont int N = atoi(argv[1])
                type_decl = c_ast.TypeDecl(k, [], c_ast.IdentifierType(['int']))
                init = c_ast.FuncCall(
                    c_ast.ID('atoi'),
                    c_ast.ExprList([c_ast.ArrayRef(c_ast.ID('argv'), c_ast.Constant('int', str(i)))]))
                i += 1
                #add the variable to the list of variables
                type_name = c_ast.Typename(None, [], type_decl)
                sizes_decls_typenames.append(type_name)
                decl = c_ast.Decl(k, ['const'], [], [], type_decl, init, None)
                ast.block_items.insert(0, decl)

            #build and add the if statement checking the number of arguments passed on the command line
            mysize = 'size ' * len(self.constants)
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
                type_decl = c_ast.TypeDecl(str(name), [], c_ast.IdentifierType(['int']))
                decl = c_ast.Decl(
                    name, ['const'], [], [],
                        #type_decl, c_ast.Constant('int', str(value)), None)
                        type_decl, c_ast.Constant('int', str(name.name+'_MAX')), None)
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
                            counter_name, [], c_ast.IdentifierType(['int'])),
                        c_ast.Constant('int', '0'), None)], None)

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
                    #, rand_max))
                    #,
                    # c_ast.Constant('float', factor))
                else:
                    #array a or b
                    rand_max = c_ast.BinaryOp(
                            '/', c_ast.FuncCall(
                                c_ast.ID('rand'),
                                c_ast.ExprList([])),
                            c_ast.Cast(
                                c_ast.IdentifierType(['double']),
                                c_ast.ID('RAND_MAX')))

                stmt = c_ast.Assignment(
                    '=',
                    c_ast.ArrayRef(c_ast.ID(d.name), c_ast.ID(counter_name)),
                    rand_max)
                ast.block_items.insert(
                    i + 1, c_ast.For(init, cond, next_, stmt))

                # inject dummy access to arrays, so compiler does not over-optimize code
                # with if around it, so code will actually run
                # ast.block_items.insert(
                #     i + 2, c_ast.If(
                #         cond=c_ast.ID('var_false'),
                #         iftrue=c_ast.Compound([
                #             c_ast.FuncCall(
                #                 c_ast.ID('dummy'),
                #                 c_ast.ExprList([c_ast.ID(d.name)]))]),
                #         iffalse=None))
            else:
                # this is a scalar, so a simple Assignment is enough

                #calculate the factor
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

                # inject dummy access to scalar, so compiler does not over-optimize code
                # TODO put if around it, so code will actually run
                # ast.block_items.insert(
                #     i + 2, c_ast.If(
                #         cond=c_ast.ID('var_false'),
                #         iftrue=c_ast.Compound([
                #             c_ast.FuncCall(
                #                 c_ast.ID('dummy'),
                #                 c_ast.ExprList([c_ast.UnaryOp('&', c_ast.ID(d.name))]))]),
                #         iffalse=None))

        # transform multi-dimensional array references to one dimensional
        # references
        list(map(lambda aref: transform_multidim_to_1d_ref(aref, array_dimensions),
                 find_array_references(ast)))

        dummies = []
        # Make sure nothing gets removed by inserting dummy calls
        for d in declarations:
            if array_dimensions[d.name]:
                dummies.append(c_ast.If(
                    cond=c_ast.ID('var_false'),
                    iftrue=c_ast.Compound([
                        c_ast.FuncCall(
                            c_ast.ID('dummy'),
                            c_ast.ExprList([c_ast.ID(d.name)]))]),
                    iffalse=None))
            else:
                dummies.append(c_ast.If(
                    cond=c_ast.ID('var_false'),
                    iftrue=c_ast.Compound([
                        c_ast.FuncCall(
                            c_ast.ID('dummy'),
                            c_ast.ExprList([c_ast.UnaryOp('&', c_ast.ID(d.name))]))]),
                    iffalse=None))

        # if we do not want the version accepting inputs from command line,
        # we need to declare the blocking factor
        if not from_cli:
            # add declaration of the block
            if self.block_factor:
                type_decl = c_ast.TypeDecl(
                    'block_factor', [], c_ast.IdentifierType(['int']))
                decl = c_ast.Decl(
                    'block_factor', ['const'], [], [],
                    type_decl, c_ast.Constant('int', str(self.block_factor)), None)
                ast.block_items.insert(-3, decl)

                # add it to the list of declarations, so it gets passed to the
                # kernel_loop
                declarations.append(decl)

        # Wrap everything in a loop
        # int repeat = atoi(argv[2])
        type_decl = c_ast.TypeDecl(
            'repeat', [], c_ast.IdentifierType(['int']))
        # init = c_ast.FuncCall(
        #     c_ast.ID('atoi'),
        #     c_ast.ExprList([c_ast.ArrayRef(
        #         c_ast.ID('argv'), c_ast.Constant('int', str(len(self.constants)+2)))]))
        # ast.block_items.insert(-3, c_ast.Decl(
        #     'repeat', ['const'], [], [],
        #     type_decl, init, None))
        ast.block_items.insert(-3, c_ast.Decl(
            'repeat', ['const'], [], [],
            type_decl, c_ast.Constant('int', '1'), None))

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
                    index_name, [], c_ast.IdentifierType(['int'])),
                c_ast.Constant('int', '0'),
                None)], None)
        cond = c_ast.BinaryOp('<', c_ast.ID(
            index_name), c_ast.ID('repeat'))
        next_ = c_ast.UnaryOp('++', c_ast.ID(index_name))
        #stmt = c_ast.Compound([ast.block_items.pop(-2)]+dummies)

        expr_list = [c_ast.ID(d.name) for d in declarations] + [c_ast.ID(s.name) for s in self.constants]
        if self.block_factor:
            expr_list = expr_list + [c_ast.ID('block_factor')]

        stmt = c_ast.FuncCall(c_ast.ID('kernel_loop'),
                              c_ast.ExprList(expr_list))

        swap_tmp = c_ast.Assignment('=', c_ast.ID('tmp'),
                                    c_ast.ID(pointers_list[0].type.type.declname))
        swap_grid = c_ast.Assignment('=', c_ast.ID(pointers_list[0].type.type.declname),
                                     c_ast.ID(pointers_list[1].type.type.declname))
        last_swap = c_ast.Assignment('=', c_ast.ID(pointers_list[1].type.type.declname),
                                     c_ast.ID('tmp'))
        stmt = c_ast.Compound([stmt, swap_tmp, swap_grid, last_swap])
        myfor = c_ast.For(init, cond, next_, stmt)

        # call the timing function at the beginning
        end_timing = c_ast.FuncCall(c_ast.ID('timing'),
                                    c_ast.ExprList([c_ast.UnaryOp('&', c_ast.ID('wct_end')),
                                                    c_ast.UnaryOp('&', c_ast.ID('cput_end'))]))

        update_runtime = c_ast.Assignment('=', c_ast.ID('runtime'),
                                          c_ast.BinaryOp('-', c_ast.ID('wct_end'), c_ast.ID('wct_start')))

        update_iter = c_ast.Assignment('*=', c_ast.ID('repeat'),
                                       c_ast.Constant('int', '2'))

        # while(runtime<2. || repeat<=2) {...}
        cond = c_ast.BinaryOp( '||',
                c_ast.BinaryOp('<', c_ast.ID( 'runtime'), c_ast.Constant('double', '2.0')),
                c_ast.BinaryOp('<=', c_ast.ID( 'repeat'), c_ast.Constant('int', '2')));
        stmt = c_ast.Compound(
            [start_timing, myfor, end_timing, update_runtime, update_iter])

        ast.block_items.insert(-1, c_ast.While(cond, stmt))

        # the variable repeat must be divided by 2 since in the last loop
        # was doubled before exiting
        ast.block_items.insert(-1, c_ast.Assignment('/=',
                                                    c_ast.ID('repeat'), c_ast.Constant('int', '2')))

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
                    run_index_name, [], c_ast.IdentifierType(['int'])),
                c_ast.Constant('int', '0'),
                None)], None)
        run_cond = c_ast.BinaryOp('<', c_ast.ID(run_index_name), c_ast.ID('repeat'))
        run_next_ = c_ast.UnaryOp('++', c_ast.ID(run_index_name))
        #run_stmt = c_ast.Compound([ast.block_items.pop(-2)]+dummies)

        run_expr_list = [c_ast.ID(d.name) for d in declarations] + [c_ast.ID(s.name) for s in self.constants]
        if self.block_factor:
            run_expr_list = expr_list + [c_ast.ID('block_factor')]

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
        run_stmt = c_ast.Compound([run_stmt, run_swap_tmp, run_swap_grid, run_last_swap])
        run_myfor = c_ast.For(run_init, run_cond, run_next_, run_stmt)
        ast.block_items.insert(-1, run_myfor)

        if type_ == 'likwid':
            # close the region "Sweep" of likwid
            ast.block_items.insert(-1,
                                   c_ast.Constant('string', 'INSERTMACROSTOP'))

        # calculate the size of the grid, taking the letters representing
        # its dimensions from the array of constants
        size = '(' + ' * '.join(k.name for k in self.constants) + ')'

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


        # else:
        #     ast.block_items += dummies

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
                            beginning, [], c_ast.IdentifierType(['int'])),
                        myblockstmt.init.decls[0].init, None)], None)
                # for(jb = 1; jb < N-1; jb+=block_factor) {...}reduce(lambda l,
                # r: c_ast.BinaryOp('*', l, r), array_dimensions[d.name]))
                cond = c_ast.BinaryOp('<', c_ast.ID(
                    beginning), myblockstmt.cond.right)
                next_ = c_ast.BinaryOp(
                    '+=', c_ast.ID(beginning), c_ast.ID('block_factor'))
                #stmt = c_ast.Compound([ast.block_items.pop(-2)]+dummies)

                decl = c_ast.Decl(end, [], [], [], c_ast.TypeDecl(
                    end, [], c_ast.IdentifierType(['int'])), c_ast.FuncCall(
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
                beginning = myvariables[1] + 'b'
                end = myvariables[1] + 'end'
                pragma = c_ast.Pragma(
                    'omp parallel')# private({}, {})'.format(beginning, end))

                init = c_ast.DeclList([
                    c_ast.Decl(
                        beginning, [], [], [], c_ast.TypeDecl(
                            beginning, [], c_ast.IdentifierType(['int'])),
                        myblockstmt.init.decls[0].init,
                        None)], None)
                # for(jb = 1; jb < N-1; jb+=block_factor) {...}reduce(lambda l,
                # r: c_ast.BinaryOp('*', l, r), array_dimensions[d.name]))
                cond = c_ast.BinaryOp('<', c_ast.ID(
                    beginning), myblockstmt.cond.right)
                next_ = c_ast.BinaryOp(
                    '+=', c_ast.ID(beginning), c_ast.ID('block_factor'))
                #stmt = c_ast.Compound([ast.block_items.pop(-2)]+dummies)

                decl = c_ast.Decl(end, [], [], [], c_ast.TypeDecl(
                    end, [], c_ast.IdentifierType(['int'])), c_ast.FuncCall(
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

        code = '#include <math.h>\n\n' + code
        code = '#include <stdlib.h>\n' + code

        # substitute the string added with the macro, since there is no way to
        # add MACROs with pycparser. It is a workaround
        # TODO change the code creation in a way to use MACROs
        pragraomp = '  #pragma omp parallel\n  {}\n    ' + '{}' + '\n  {}'

        likwid_init = 'LIKWID_MARKER_INIT;'
        likwid_register = 'LIKWID_MARKER_REGISTER("Sweep");'
        macroinit = '\n  ' + ifdefperf + '  ' + likwid_init + '\n  '
        macroinit += likwid_register + '\n  ' + endif
        code = code.replace('INSERTMACROINIT;', macroinit)

        start_sweep = 'LIKWID_MARKER_START("Sweep");'
        pragma_start_sweep = pragraomp.format('{', start_sweep, '}')
        macrostart = '\n  ' + ifdefperf + pragma_start_sweep + '\n  ' + endif
        code = code.replace('INSERTMACROSTART;', macrostart)

        stop_sweep = 'LIKWID_MARKER_STOP("Sweep");'
        pragma_stop_sweep = pragraomp.format('{', stop_sweep, '}')
        macrostop = '\n  ' + ifdefperf + pragma_stop_sweep + '\n  ' + endif
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
        kernel = '#endif\n' + kernel
        kernel = '#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )\n' + kernel
        kernel = '#ifndef min\n' + kernel

        # return mycode
        return code, kernel

