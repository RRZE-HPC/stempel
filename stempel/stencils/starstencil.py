from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import

import sys
import string


def signum(number=1):
    if type(number)==unicode or type(number)==str:
        number = int(number)
    if number < 0:
        sig = '-'
    else:
        sig = '+'
    return sig

def value(number = 1):
    if type(number)==unicode or type(number)==str:
        number = int(number)
    if abs(number) != 1:
        mystring = str(abs(number))+'*'
    else:
        mystring = ''
    return mystring

def left(centerpoint='a[j][i]', dimension=2, loop_variables=['j', 'i'], myrange=1):
    letter = loop_variables[dimension-1]
    newpoint = centerpoint.replace(letter, letter+'-{}'.format(abs(myrange)))
    return newpoint

def right(centerpoint='a[j][i]', dimension=2, loop_variables=['j', 'i'], myrange=1):
    letter = loop_variables[dimension-1]
    newpoint = centerpoint.replace(letter, letter+'+{}'.format(abs(myrange)))
    return newpoint
               
class Star(object):

    name = "star"

    @classmethod
    def configure_arggroup(cls, parser):
        pass

    def __init__(self, dimensions=2, simmetricity=((1,1,1), (1,1,1)), coeff=None , datatype = 'double', args=None, parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to 2 (2 dimensional stencil)
        *sizes* is the size of the stencil in each of the dimensions. It defaults to 50 in x and 50 in y
        *simmetricity* is a tuple representing the grid points from which the stencil depends. values of the stencil. For each dimension
        we have a tuple representing the coefficients of the neighbours in that side of the subdimension (left or right). A tuple consists
        always of an odd number of values: if a value equals to 0, it means that that neighbour does not play a role in the stencil computation
        *coeff* represents the coefficient of the stencil, in case they are not a constant, thus it is not possible to specify them in the simmetricity input.
        *datatype* represents the type of the data to be store in the grids. By default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions

        self.dims = []
        #save the letter of the dimensions in a variable (if 2 dimensions they are "M" and "N")
        for i in range(0, self.dimensions):
            self.dims.append(string.ascii_uppercase[12+i])

        self.simmetricity = simmetricity
        self.coeff = coeff
        self.datatype = datatype
        
        #to be changed in future to allow stencil on more than 2 grids
        self.output = string.ascii_lowercase[1]
        self.input = string.ascii_lowercase[0]

        #if there is a matrix of coefficient we name it after the third letter of the alphabet
        if self.coeff=='matrix':
            self.coefficient =  string.ascii_lowercase[2]
        else:
            self.coefficient = 's'

        #get the indices saved in a variable (if 2 dimension they are "j" and "i")
        self.loop_variables = []
        for i in reversed(range(0, self.dimensions)):
            self.loop_variables.append(string.ascii_lowercase[8+i])
  
        self._args = args
        self._parser = parser

        if args:
            # handle CLI info
            pass

    def declaration(self):
        '''
        builds the declaration of thevariables for the C code and returns declaration as unicode
        
        '''
        #initialization of the input matrix
        init_line1 = self.datatype + ' a'
        #initialization of the output matrix
        init_line2 = self.datatype + ' b'
        #add the dimensions to the matrices
        for i in range(0, self.dimensions):
            init_line1 = init_line1 + '[' + self.dims[i] + ']'
            init_line2 = init_line2 + '[' + self.dims[i] + ']'

        init_line1 = init_line1 + ';'
        init_line2 = init_line2 + ';'

        #declare a variable for the initialization of the matrix of coefficients
        init_line3=None
        #if a matrix of coefficients exists, then we declare it
        if self.coeff=='matrix':
            split = init_line2.split()
            init_line3 =  split[0] + ' c' + split[1][1:]
        elif self.coeff=='scalar':
            init_line3 = self.datatype + ' ' + self.coefficient + ';'

        declaration = init_line1 + '\n' + init_line2 + '\n'

        if init_line3:
            declaration = declaration + init_line3 + '\n'

        return declaration



    def loop(self):
        '''
        builds the loop for the C code and returns loop_lines as list
        
        '''
        loop_lines=[]

        #build the lines of the foor loop, according to the dimensions we have
        for i in range (0, self.dimensions):
            line = 'for(int {}={}; {} < {}-{}; {}++)'.format(self.loop_variables[i], len(self.simmetricity[i])//2, self.loop_variables[i], self.dims[i], len(self.simmetricity[i])//2, self.loop_variables[i]) + '{'
            loop_lines.insert(i, line)

        centerpoint = self.input
        lefthand = self.output
        coefficient = self.coefficient
        # build the centerpoint and the lefthand of the equation
        for i in range(0, self.dimensions):
            lefthand = lefthand + '[' + self.loop_variables[i] + ']'
            centerpoint = centerpoint + '[' + self.loop_variables[i] + ']'

        # declare an empty stencil line (string)
        stencil = ''
        #count if the coefficient of the centerpoint in each dimension is not 0
        count = 0
        # store how many of the coefficient of the centerpoint (in each axis) are negative
        sig = 0
        #store the product of all the coefficient of the centerpoint in each axis
        centercoeff = 1
        # for each of the dimensions of the simmetricity (==dimensions)
        for i in range(len(self.simmetricity)):
            # take each value (simmetricity[0] == 1 1 1)
            for k in range(len(self.simmetricity[i])):
                # each element which is in the left part of the array, it means that is a value on the left of the centerpoint in the grid and is not 0 (otherwise we would not consider it)
                if k < len(self.simmetricity[i])//2 and int(self.simmetricity[i][k]) != 0:
                    # build the value as a gridpoint
                    stencil = stencil + '{}{}{} '.format(signum(self.simmetricity[i][k]), value(self.simmetricity[i][k]), left(centerpoint, i, self.loop_variables, k-len(self.simmetricity[i])//2))
                # each element which is in the right part of the array, it means that is a value on the right of the centerpoint in the grid and is not 0 (otherwise we would not consider it)
                elif k > len(self.simmetricity[i])//2 and int(self.simmetricity[i][k]) != 0:
                    stencil = stencil + '{}{}{} '.format(signum(self.simmetricity[i][k]), value(self.simmetricity[i][k]), right(centerpoint, i, self.loop_variables, k-len(self.simmetricity[i])//2))

                assert bool(len(self.simmetricity[i]) & 1), "each tuple of simmetricity must contain an odd number of values (always with centerpoint)"
                if k == len(self.simmetricity[i])//2 and int(self.simmetricity[i][k]) != 0:
                    count = count + 1
                    if signum(self.simmetricity[i][k]) == '-':
                        sig = sig + 1
                    centercoeff = centercoeff * int(self.simmetricity[i][k])

        if count == len(self.simmetricity):
            # if sig is odd
            if  bool(sig & 1):
                centersig = '-'
            else:
                centersig = '+'

            stencil = stencil + '{}{}{} '.format(centersig, value(centercoeff), centerpoint)
        

        
        scalar = ' * {}'.format(self.coefficient)

        #remove trailing +
        if stencil[0] == "+":
            stencil = stencil[1:]
        

        righthand = '({}){};'.format(stencil, scalar)

        closing = '}\n' * self.dimensions

        computation = lefthand + ' = ' + righthand + '\n' + closing

        loop_lines.append(computation)

        return loop_lines




class StarMat(object):

    name = "starMat"

    @classmethod
    def configure_arggroup(cls, parser):
        pass

    def __init__(self, dimensions=2, simmetricity=((1,1,1), (1,1,1)), coeff=None , datatype = 'double', args=None, parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to 2 (2 dimensional stencil)
        *sizes* is the size of the stencil in each of the dimensions. It defaults to 50 in x and 50 in y
        *simmetricity* is a tuple representing the grid points from which the stencil depends. values of the stencil. For each dimension
        we have a tuple representing the coefficients of the neighbours in that side of the subdimension (left or right). A tuple consists
        always of an odd number of values: if a value equals to 0, it means that that neighbour does not play a role in the stencil computation
        *coeff* represents the coefficient of the stencil, in case they are not a constant, thus it is not possible to specify them in the simmetricity input.
        *datatype* represents the type of the data to be store in the grids. By default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions

        self.dims = []
        #save the letter of the dimensions in a variable (if 2 dimensions they are "M" and "N")
        for i in range(0, self.dimensions):
            self.dims.append(string.ascii_uppercase[12+i])

        self.simmetricity = simmetricity
        self.coeff = coeff
        self.datatype = datatype
        
        #to be changed in future to allow stencil on more than 2 grids
        self.output = string.ascii_lowercase[1]
        self.input = string.ascii_lowercase[0]

        #if there is a matrix of coefficient we name it after the third letter of the alphabet
        if self.coeff=='matrix':
            self.coefficient =  string.ascii_lowercase[2]
        else:
            self.coefficient = 's'

        #get the indices saved in a variable (if 2 dimension they are "j" and "i")
        self.loop_variables = []
        for i in reversed(range(0, self.dimensions)):
            self.loop_variables.append(string.ascii_lowercase[8+i])
  
        self._args = args
        self._parser = parser

        if args:
            # handle CLI info
            pass

    def declaration(self):
        '''
        builds the declaration of thevariables for the C code and returns declaration as unicode
        
        '''
        #initialization of the input matrix
        init_line1 = self.datatype + ' a'
        #initialization of the output matrix
        init_line2 = self.datatype + ' b'
        #add the dimensions to the matrices
        for i in range(0, self.dimensions):
            init_line1 = init_line1 + '[' + self.dims[i] + ']'
            init_line2 = init_line2 + '[' + self.dims[i] + ']'

        init_line1 = init_line1 + ';'
        init_line2 = init_line2 + ';'

        #declare a variable for the initialization of the matrix of coefficients
        init_line3=None
        #if a matrix of coefficients exists, then we declare it
        if self.coeff=='matrix':
            split = init_line2.split()
            init_line3 =  split[0] + ' c' + split[1][1:]
        elif self.coeff=='scalar':
            init_line3 = self.datatype + ' ' + self.coefficient + ';'

        declaration = init_line1 + '\n' + init_line2 + '\n'

        if init_line3:
            declaration = declaration + init_line3 + '\n'

        return declaration


    def loop(self):
        '''
        builds the loop for the C code and returns loop_lines as list
        
        '''
        loop_lines=[]

        #build the lines of the foor loop, according to the dimensions we have
        for i in range (0, self.dimensions):
            line = 'for(int {}={}; {} < {}-{}; {}++)'.format(self.loop_variables[i], len(self.simmetricity[i])//2, self.loop_variables[i], self.dims[i], len(self.simmetricity[i])//2, self.loop_variables[i]) + '{'
            loop_lines.insert(i, line)

        centerpoint = self.input
        lefthand = self.output
        coefficient = self.coefficient
        # build the centerpoint and the lefthand of the equation
        for i in range(0, self.dimensions):
            lefthand = lefthand + '[' + self.loop_variables[i] + ']'
            centerpoint = centerpoint + '[' + self.loop_variables[i] + ']'
            coefficient = coefficient + '[' + self.loop_variables[i] + ']'

        # declare an empty stencil line (string)
        stencil = ''
        #count if the coefficient of the centerpoint in each dimension is not 0
        count = 0
        # store how many of the coefficient of the centerpoint (in each axis) are negative
        sig = 0
        #store the product of all the coefficient of the centerpoint in each axis
        centercoeff = 1
        # for each of the dimensions of the simmetricity (==dimensions)
        for i in range(len(self.simmetricity)):
            # take each value (simmetricity[0] == 1 1 1)
            for k in range(len(self.simmetricity[i])):
                # each element which is in the left part of the array, it means that is a value on the left of the centerpoint in the grid and is not 0 (otherwise we would not consider it)
                if k < len(self.simmetricity[i])//2 and int(self.simmetricity[i][k]) != 0:
                    # build the value as a gridpoint
                    stencil = stencil + '{}{}{}{}{} '.format(signum(self.simmetricity[i][k]), value(self.simmetricity[i][k]), left(coefficient, i, self.loop_variables, k-len(self.simmetricity[i])//2), ' * ', left(centerpoint, i, self.loop_variables, k-len(self.simmetricity[i])//2))
                # each element which is in the right part of the array, it means that is a value on the right of the centerpoint in the grid and is not 0 (otherwise we would not consider it)
                elif k > len(self.simmetricity[i])//2 and int(self.simmetricity[i][k]) != 0:
                    stencil = stencil + '{}{}{}{}{} '.format(signum(self.simmetricity[i][k]), value(self.simmetricity[i][k]), right(coefficient, i, self.loop_variables, k-len(self.simmetricity[i])//2), ' * ', right(centerpoint, i, self.loop_variables, k-len(self.simmetricity[i])//2))

                assert bool(len(self.simmetricity[i]) & 1), "each tuple of simmetricity must contain an odd number of values (always with centerpoint)"
                if k == len(self.simmetricity[i])//2 and int(self.simmetricity[i][k]) != 0:
                    count = count + 1
                    if signum(self.simmetricity[i][k]) == '-':
                        sig = sig + 1
                    centercoeff = centercoeff * int(self.simmetricity[i][k])

        if count == len(self.simmetricity):
            # if sig is odd
            if  bool(sig & 1):
                centersig = '-'
            else:
                centersig = '+'

            stencil = stencil + '{}{}{}{}{} '.format(centersig, value(centercoeff), coefficient, ' * ', centerpoint)
        

        #remove trailing +
        if stencil[0] == "+":
            stencil = stencil[1:]
        

        righthand = '({});'.format(stencil)

        closing = '}\n' * self.dimensions

        computation = lefthand + ' = ' + righthand + '\n' + closing

        loop_lines.append(computation)

        return loop_lines

