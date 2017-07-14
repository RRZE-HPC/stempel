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

            
class StarConstant(object):
    """class for the stencil with constant coefficients
    It can be:
        simmetric->isotropic
        simmetric->anisotropic
        asimmetric->anisotropic
    """

    name = "star"

    @classmethod
    def configure_arggroup(cls, parser):
        pass

    def __init__(self, dimensions=2, radius=1, simmetricity=True, isotropy=True , datatype = 'double', inputgrids=1, args=None, parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to 2 (2 dimensional stencil)
        *radius* represents the radius of the stencil on each side of each dimension
        *simmetricity* is a boolean representing the simmetricity of the stencil with respect to the coefficients
        *isotropy* is a boolean representing the isotropy of the stencil (no dependency on the direction)
        *coeff* represents the coefficients of the stencil: can be either constant or variable.
        *datatype* represents the type of the data to be store in the grids. By default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions

        self.dims = []
        #save the letter of the dimensions in a variable (if 2 dimensions they are "M" and "N")
        for i in range(0, self.dimensions):
            self.dims.append(string.ascii_uppercase[12+i])

        self.radius = radius
        self.simmetricity = simmetricity
        self.isotropy = isotropy
        self.datatype = datatype
        
        self.inputgrids = inputgrids
        #to be changed in future to allow stencil on more than 2 grids
        self.inputs = [string.ascii_lowercase[i] for i in range(inputgrids)]
        self.output = string.ascii_lowercase[inputgrids]
        

        # #if there is a matrix of coefficients we name it after the next available letter of the alphabet
        # if self.coeff=='variable':
        #     self.coefficient =  string.ascii_uppercase[inputgrids+1]
        # #else we name it with a constant w
        # else:
        #     self.coefficient = string.ascii_lowercase[inputgrids+1]
        if self.isotropy and self.simmetricity:
            num_coefficients = radius+1
        elif self.simmetricity and not self.isotropy:
            num_coefficients = (radius * self.dimensions) + 1
        else:#not simmetricity and not isotropy)
            num_coefficients = (2 * radius * self.dimensions) + 1

        # myformat = '{}' * num_coefficients
        myletter = string.ascii_lowercase[inputgrids+1]
        self.coefficients = [ myletter+str(i) for i in range(num_coefficients)]


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
        #initialization of the input matrices (array of input matrices)
        init_inputs = []
        for i in self.inputs:
            init_inputs.append(self.datatype + ' ' + i)
            #init_inputs = [u'double a']

        #initialization of the output matrix
        init_output = self.datatype + ' ' + self.output
        
        #add the dimensions to the matrices
        for lineno in range(len(init_inputs)):
            for i in range(0, self.dimensions):
                line = init_inputs[lineno] + '[' + self.dims[i] + ']'
                init_inputs.pop(lineno)
                init_inputs.insert(lineno, line)
                init_output = init_output + '[' + self.dims[i] + ']'
        

        for lineno in range(len(init_inputs)):
                line = init_inputs[lineno] + ';'
                init_inputs.pop(lineno)
                init_inputs.insert(lineno, line)
        init_output = init_output + ';'

        #declare a variable for the initialization of the coefficients. It is a string
        init_coefficients=''
        #we declare all the coefficients line by line (can be modified to make only a 1 line)
        for i in self.coefficients:
            init_coefficients = init_coefficients + self.datatype + ' ' + i + ';\n'

        declaration = ''
        for i in init_inputs:
            declaration = declaration + i + '\n'

        declaration = declaration + init_output + '\n' + init_coefficients + '\n'

        return declaration



    def loop(self):
        '''
        builds the loop for the C code and returns loop_lines as list
        
        '''
        loop_lines=[]

        #build the lines of the foor loop, according to the dimensions we have
        for i in range (0, self.dimensions):
            line = 'for(int {}={}; {} < {}-{}; {}++)'.format(self.loop_variables[i], self.radius, self.loop_variables[i], self.dims[i], self.radius, self.loop_variables[i]) + '{'
            loop_lines.insert(i, line)

        centerpoint = self.inputs[0]
        lefthand = self.output
        coefficients = self.coefficients
        radius = self.radius
        dimensions = self.dimensions

        # build the centerpoint and the lefthand of the equation
        for i in range(0, dimensions):
            lefthand = lefthand + '[' + self.loop_variables[i] + ']'
            centerpoint = centerpoint + '[' + self.loop_variables[i] + ']'

        # declare an empty stencil line (string)
        stencil = ''

        # isotropic and simmetric
        if self.simmetricity and self.isotropy:
            assert (len(self.coefficients) == (self.radius + 1)), "In case of an isotropic and simmetric stencil with constant coefficient, the number of the coefficient must be equal to (radius + 1)"
            stencil = self.coefficients[0] + '*' + centerpoint + '\n'
            for i in range(self.radius):
                for j in range(self.dimensions):
                    stencil = stencil + '+ {} * ({} + {})'.format(self.coefficients[i+1], left(centerpoint, j, self.loop_variables, i+1), right(centerpoint, j, self.loop_variables, i+1)) + '\n'
        
        # asimmetric
        elif not self.simmetricity:
            assert (len(self.coefficients) == (2 * self.radius * self.dimensions + 1)), "In case of an asimmetric stencil with constant coefficient, the number of the coefficient must be equal to (2 * radius * dimensions + 1)"
            stencil = self.coefficients[0] + '*' + centerpoint + '\n'
            for i in range(self.radius):
                for j in range(self.dimensions):
                    stencil = stencil + '+ {} * {} + {} * {}'.format(coefficients[j + (i+1)*(j+1)], left(centerpoint, j, self.loop_variables, i+1), coefficients[j + (i+1)*(j+1) + 1], right(centerpoint, j, self.loop_variables, i+1)) + '\n'
        
        # anisotropic and simmetric
        elif not self.isotropy and self.simmetricity:
            assert (len(self.coefficients) == (self.radius * self.dimensions + 1)), "In case of anisotropic and simmetric stencil with constant coefficient, the number of the coefficient must be equal to (radius * dimensions + 1)"
            stencil = self.coefficients[0] + '*' + centerpoint + '\n'
            for i in range(self.radius):
                for j in range(self.dimensions):
                    stencil = stencil + '+ {} * ({} + {})'.format(self.coefficients[(i+1)*(j+1)], left(centerpoint, j, self.loop_variables, i+1), right(centerpoint, j, self.loop_variables, i+1)) + '\n'
        

        righthand = '({});'.format(stencil)

        closing = '}\n' * self.dimensions

        computation = lefthand + ' = ' + righthand + '\n' + closing

        loop_lines.append(computation)

        return loop_lines



 

class StarVariable(object):
    """class for the stencil with variable coefficients
    It can be:
        simmetric->isotropic
        simmetric->anisotropic
        asimmetric->anisotropic
    """
    name = "starVariable"

    @classmethod
    def configure_arggroup(cls, parser):
        pass
    
    def __init__(self, dimensions=2, radius=1, simmetricity=True, isotropy=True , datatype = 'double', inputgrids=1, args=None, parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to 2 (2 dimensional stencil)
        *radius* represents the radius of the stencil on each side of each dimension
        *simmetricity* is a boolean representing the simmetricity of the stencil with respect to the coefficients
        *isotropy* is a boolean representing the isotropy of the stencil (no dependency on the direction)
        *coeff* represents the coefficients of the stencil: can be either constant or variable.
        *datatype* represents the type of the data to be store in the grids. By default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions

        self.dims = []
        #save the letter of the dimensions in a variable (if 2 dimensions they are "M" and "N")
        for i in range(0, self.dimensions):
            self.dims.append(string.ascii_uppercase[12+i])

        self.radius = radius
        self.simmetricity = simmetricity
        self.isotropy = isotropy
        self.datatype = datatype
        
        self.inputgrids = inputgrids
        #to be changed in future to allow stencil on more than 2 grids
        self.inputs = [string.ascii_lowercase[i] for i in range(inputgrids)]
        self.output = string.ascii_lowercase[inputgrids]

        if self.isotropy and self.simmetricity:
            num_coefficients = radius+1
        elif self.simmetricity and not self.isotropy:
            num_coefficients = (radius * self.dimensions) + 1
        else:#not simmetricity and not isotropy)
            num_coefficients = (2 * radius * self.dimensions) + 1

        # myformat = '{}' * num_coefficients
        myletter = string.ascii_uppercase[inputgrids+1]
        self.coefficients = [ myletter+str(i) for i in range(num_coefficients)]

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
        #initialization of the input matrices (array of input matrices)
        init_inputs = []
        for i in self.inputs:
            init_inputs.append(self.datatype + ' ' + i)
            #init_inputs = [u'double a']

        #initialization of the output matrix
        init_output = self.datatype + ' ' + self.output

        init_coefficients = []
        for i in self.coefficients:
            init_coefficients.append(self.datatype + ' ' + i)
        
        #add the dimensions to the matrices
        for lineno in range(len(init_inputs)):
            for i in range(0, self.dimensions):
                line = init_inputs[lineno] + '[' + self.dims[i] + ']'
                init_inputs.pop(lineno)
                init_inputs.insert(lineno, line)
                init_output = init_output + '[' + self.dims[i] + ']'

        for lineno in range(len(init_coefficients)):
            for i in range(0, self.dimensions):
                line = init_coefficients[lineno] + '[' + self.dims[i] + ']'
                init_coefficients.pop(lineno)
                init_coefficients.insert(lineno, line)

        # add ";" to close the line
        for lineno in range(len(init_inputs)):
                line = init_inputs.pop(lineno) + ';'
                init_inputs.insert(lineno, line)
        init_output = init_output + ';'

        for lineno in range(len(init_coefficients)):
            line = init_coefficients.pop(lineno) + ';'
            init_coefficients.insert(lineno, line)
        
        declaration = ''
        for i in init_inputs:
            declaration = declaration + i + '\n'

        declaration = declaration + init_output + '\n'

        for i in init_coefficients:
            declaration = declaration + i + '\n\n'

        return declaration



    def loop(self):
        '''
        builds the loop for the C code and returns loop_lines as list
        
        '''
        loop_lines=[]

        #build the lines of the foor loop, according to the dimensions we have
        for i in range (0, self.dimensions):
            line = 'for(int {}={}; {} < {}-{}; {}++)'.format(self.loop_variables[i], self.radius, self.loop_variables[i], self.dims[i], self.radius, self.loop_variables[i]) + '{'
            loop_lines.insert(i, line)

        centerpoint = self.inputs[0]
        lefthand = self.output
        coefficients = self.coefficients
        radius = self.radius
        dimensions = self.dimensions

        # build the centerpoint, the lefthand and the coefficients of the equation
        for i in range(0, dimensions):
            lefthand = lefthand + '[' + self.loop_variables[i] + ']'
            centerpoint = centerpoint + '[' + self.loop_variables[i] + ']'
            for c in range(0, len(coefficients)):
                coeff = coefficients.pop(c)
                coeff = coeff + '[' + self.loop_variables[i] + ']'
                coefficients.insert(c, coeff)

        # declare an empty stencil line (string)
        stencil = ''

        # isotropic and simmetric
        if (self.simmetricity and self.isotropy):
            print("simmetric and isotropic")
            assert (len(self.coefficients) == (self.radius + 1)), "In case of an isotropic and simmetric stencil with constant coefficient, the number of the coefficient must be equal to (radius + 1)"
            stencil = self.coefficients[0] + '*' + centerpoint + '\n'
            for i in range(self.radius):
                for j in range(self.dimensions):
                    stencil = stencil + '+ {} * ({} + {})'.format(self.coefficients[i+1], left(centerpoint, j, self.loop_variables, i+1), right(centerpoint, j, self.loop_variables, i+1)) + '\n'
        
        # asimmetric
        elif not self.simmetricity:
            print("asimmetric")
            assert (len(self.coefficients) == (2 * self.radius * self.dimensions + 1)), "In case of an asimmetric stencil with constant coefficient, the number of the coefficient must be equal to (2 * radius * dimensions + 1)"
            stencil = self.coefficients[0] + '*' + centerpoint + '\n'
            for i in range(self.radius):
                for j in range(self.dimensions):
                    stencil = stencil + '+ {} * {} + {} * {}'.format(coefficients[j + (i+1)*(j+1)], left(centerpoint, j, self.loop_variables, i+1), coefficients[j + (i+1)*(j+1) + 1], right(centerpoint, j, self.loop_variables, i+1)) + '\n'
        
        # anisotropic and simmetric
        elif (not self.isotropy and self.simmetricity):
            print("simmetric and anisotropic")
            assert (len(self.coefficients) == (self.radius * self.dimensions + 1)), "In case of anisotropic and simmetric stencil with constant coefficient, the number of the coefficient must be equal to (radius * dimensions + 1)"
            stencil = self.coefficients[0] + '*' + centerpoint + '\n'
            for i in range(self.radius):
                for j in range(self.dimensions):
                    stencil = stencil + '+ {} * ({} + {})'.format(self.coefficients[(i+1)*(j+1)], left(centerpoint, j, self.loop_variables, i+1), right(centerpoint, j, self.loop_variables, i+1)) + '\n'
        

        righthand = '({});'.format(stencil)

        closing = '}\n' * self.dimensions

        computation = lefthand + ' = ' + righthand + '\n' + closing

        loop_lines.append(computation)

        return loop_lines

