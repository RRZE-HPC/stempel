from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import

import subprocess
from functools import reduce
import operator
import sys
from distutils.spawn import find_executable
from pprint import pprint
import re
import string
import six


class Stencil(object):

    name = "star"

    @classmethod
    def configure_arggroup(cls, parser):
        pass

    def __init__(self, dimensions=2, simmetricity=((-0.5,1,2), (1,1,1)), kind='star', coeff=None , datatype = 'double', args=None, parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to 2 (2 dimensional stencil)
        *sizes* is the size of the stencil in each of the dimensions. It defaults to 50 in x and 50 in y
        *simmetricity* is a tuple representing the grid points from which the stencil depends. values of the stencil. For each dimension
        we have a tuple representing the coefficients of the neighbours in that side of the subdimension (left or right). A tuple consists
        always of an odd number of values: if a value equals to 0, it means that that neighbour does not play a role in the stencil computation
        *kind* represents the kind of stencil we deal with. It can be for instance a star or box stencil. It defaults to star
        *coeff* represents the coefficient of the stencil, in case they ae not a constant, thus it is not possible to specify them in the simmetricity input.
        *datatype* represents the type of the data to be store in the grids. By default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions
        #save the letter of the dimensions in a variable (if 2 dimensions they are "M" and "N")
        for i in range(0, self.dimensions):
            self.dims.append(string.ascii_uppercase[12+i])
        self.simmetricity = simmetricity
        self.kind = kind
        self.coeff = coeff
        self.datatype = datatype
        
        #to be changed in future to allow stencil on more than 2 grids
        self.output = string.ascii_lowercase[1]
        self.input = string.ascii_lowercase[0]

        #if there is a matrix of coefficient we name it after the third letter of the alphabet
        if self.coeff=='matrix':
            self.coeff_matrix =  string.ascii_lowercase[2]

        #get the indices saved in a variable (if 2 dimension they are "j" and "i")
        for i in reversed(range(0, self.dimensions)):
            self.loop_variables.append(string.ascii_lowercase[8+i])
  
        self._args = args
        self._parser = parser


        if args:
            # handle CLI info
            pass

    def declaration(self):
        dims = ['M', 'N', 'P']
        #initialization of the input matrix
        init_line1 = self.datatype + ' a'
        #initialization of the output matrix
        init_line2 = self.datatype + ' b'
        #add the dimensions to the matrices
        for i in range(0, dimensions):
            init_line1 = init_line1 + '[' + dims[i] + ']'
            init_line2 = init_line2 + '[' + dims[i] + ']'

        init_line1 = init_line1 + ';'
        init_line2 = init_line2 + ';'

        #declare a variable for the initialization of the matrix of coefficients
        init_line3=None
        #if a matrix of coefficients exists, then we declare it
        if self.coeff=='matrix':
            split = init_line2.split()
            init_line3 =  split[0] + ' c' + split[1][1:]
        else if self.coeff=='scalar':
            init_line3 = self.datatype + ' s;'

        declaration = init_line1 + '\n' + init_line2 + '\n'

        if init_line3:
            declaration = declaration + init_line3 + '\n'

    def left(centerpoint='a[k][j][i]', dimension=1, myrange=1):
        letter = loop_variables[dimension-1]
        newpoint = centerpoint.replace(letter, letter+'-{}'.format(abs(myrange)))
        return newpoint

    def right(centerpoint='a[k][j][i]', dimension=1, myrange=1):
        letter = loop_variables[dimension-1]
        newpoint = centerpoint.replace(letter, letter+'+{}'.format(abs(myrange)))
        return newpoint

    def loop(self):
        loop_variables = ['k', 'j', 'i']
        dims = ['M', 'N', 'P']
        loop_lines=[]
        #build the lines of the foor loop, according to the dimensions we have
        for i in range (0, self.dimensions):
            loop_lines.insert(i, 'for(int {}={}; {} < {}-{}; {}++)'.format(loop_variables[i], len(self.simmetricity[i])/2, loop_variables[i], dims[i], len(self.simmetricity[i])/2, loop_variables[i]))

        centerpoint = self.input
        lefthand = self.output
        
        for i in range(0, self.dimensions):
            lefthand = lefthand + '[' + loop_variabless[i] + ']'
            centerpoint = centerpoint + '[' + loop_variables[i] + ']'

        for i in self.simmetricity:
            for k in i:
                stencil = stencil + ' {}'.format(left(centerpoint), len(self.simmetricity), k)

        if self.coeff=='scalar':
            scalar = ' * {}'.format(self.coeff)
        else:
            scalar = ''


        righthand = '({}){};'.format(stencil, scalar)

        print(righthand)







    def perfctr(self, cmd, group='MEM', cpu='S0:0', code_markers=True, pin=True):
        '''
        runs *cmd* with likwid-perfctr and returns result as dict
        
        *group* may be a performance group known to likwid-perfctr or an event string.
        Only works with single core!
        '''

        # Making sure iaca.sh is available:
        if find_executable('likwid-perfctr') is None:
            print("likwid-perfctr was not found. Make sure likwid is installed and found in PATH.",
                  file=sys.stderr)
            sys.exit(1)

        # FIXME currently only single core measurements support!
        perf_cmd = ['likwid-perfctr', '-f', '-O', '-g', group]

        if pin:
            perf_cmd += ['-C', cpu]
        else:
            perf_cmd += ['-c', cpu]

        if code_markers:
            perf_cmd.append('-m')

        perf_cmd += cmd
        if self._args.verbose > 1:
            print(' '.join(perf_cmd))
        try:
            output = subprocess.check_output(perf_cmd).decode('utf-8').split('\n')
        except subprocess.CalledProcessError as e:
            print("Executing benchmark failed: {!s}".format(e), file=sys.stderr)
            sys.exit(1)

        results = {}
        ignore = True
        for l in output:
            l = l.split(',')
            try:
                # Metrics
                results[l[0]] = float(l[1])
            except:
                pass
            try:
                # Event counters
                counter_value = int(l[2])
                if re.fullmatch(r'[A-Z_]+', l[0]) and re.fullmatch(r'[A-Z0-9]+', l[1]):
                    results.setdefault(l[0], {})
                    results[l[0]][l[1]] = counter_value
            except (IndexError, ValueError):
                pass

        return results

    def analyze(self):
        bench = self.kernel.build(self.machine['compiler'],
                                  cflags=self.machine['compiler flags'],
                                  verbose=self._args.verbose > 1)

        # Build arguments to pass to command:
        args = [bench] + [six.text_type(s) for s in list(self.kernel.constants.values())]

        # Determan base runtime with 100 iterations
        runtime = 0.0
        time_per_repetition = 0.2/10.0

        while runtime < 0.15:
            # Interpolate to a 0.2s run
            if time_per_repetition != 0.0:
                repetitions = 0.2//time_per_repetition
            else:
                repetitions *= 10

            result = self.perfctr(args+[six.text_type(repetitions)])
            runtime = result['Runtime (RDTSC) [s]']
            time_per_repetition = runtime/float(repetitions)

        self.results = {'raw output': result}

        self.results['Runtime (per repetition) [s]'] = time_per_repetition
        # TODO make more generic to support other (and multiple) constantnames
        # TODO support SP (devide by 4 instead of 8.0)
        iterations_per_repetition = reduce(
            operator.mul,
            [self.kernel.subs_consts(max_-min_)/self.kernel.subs_consts(step)
             for idx, min_, max_, step in self.kernel._loop_stack],
            1)
        self.results['Iterations per repetition'] = iterations_per_repetition
        iterations_per_cacheline = float(self.machine['cacheline size'])/8.0
        cys_per_repetition = time_per_repetition*float(self.machine['clock'])
        self.results['Runtime (per cacheline update) [cy/CL]'] = \
            (cys_per_repetition/iterations_per_repetition)*iterations_per_cacheline
        self.results['MEM volume (per repetition) [B]'] = \
            result['Memory data volume [GBytes]']*1e9/repetitions
        self.results['Performance [MFLOP/s]'] = \
            sum(self.kernel._flops.values())/(time_per_repetition/iterations_per_repetition)/1e6
        if 'Memory bandwidth [MBytes/s]' in result:
            self.results['MEM BW [MByte/s]'] = result['Memory bandwidth [MBytes/s]']
        else:
            self.results['MEM BW [MByte/s]'] = result['Memory BW [MBytes/s]']
        self.results['Performance [MLUP/s]'] = (iterations_per_repetition/time_per_repetition)/1e6
        self.results['Performance [MIt/s]'] = (iterations_per_repetition/time_per_repetition)/1e6

    def report(self, output_file=sys.stdout):
        if self._args.verbose > 0:
            print('Runtime (per repetition): {:.2g} s'.format(
                      self.results['Runtime (per repetition) [s]']),
                  file=output_file)
        if self._args.verbose > 0:
            print('Iterations per repetition: {!s}'.format(
                     self.results['Iterations per repetition']),
                  file=output_file)
        print('Runtime (per cacheline update): {:.2f} cy/CL'.format(
                  self.results['Runtime (per cacheline update) [cy/CL]']),
              file=output_file)
        print('MEM volume (per repetition): {:.2f} Byte'.format(
                  self.results['MEM volume (per repetition) [B]']),
              file=output_file)
        print('Performance: {:.2f} MFLOP/s'.format(self.results['Performance [MFLOP/s]']),
              file=output_file)
        print('Performance: {:.2f} MLUP/s'.format(self.results['Performance [MLUP/s]']),
              file=output_file)
        print('Performance: {:.2f} It/s'.format(self.results['Performance [MIt/s]']),
              file=output_file)
        if self._args.verbose > 0:
            print('MEM bandwidth: {:.2f} MByte/s'.format(self.results['MEM BW [MByte/s]']),
                  file=output_file)
        print('', file=output_file)