#!/usr/bin/env python
"""This module defines the Stencil class of type Box with constant or variable
    coefficients and provides the functions to generate C code out of the
    stencil specification. It is used by stempel.py
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import

import string


def signum(number=1):
    """This function takes in input a number (either as unicode, string or int
    and returns its signum as a string
    """
    if isinstance(number, unicode) or isinstance(number, str):
        number = int(number)
    if number < 0:
        sig = '-'
    else:
        sig = '+'
    return sig

def value(number=1):
    """This function takes in input a number (either as unicode, string or int
    and returns its value as a string
    """
    if isinstance(number, unicode) or isinstance(number, str):
        number = int(number)
    if abs(number) != 1:
        mystring = str(abs(number))+'*'
    else:
        mystring = ''
    return mystring

def left(centerpoint='a[j][i]', dimension=2, loop_variables=None,
         myrange=1):
    """This function takes in input a point (the centerpoint of the stencil),
    the dimension (first, second or third) on which to act and the loop
    variables. It returns the point on its left, in the chosen dimension,
    as a string
    """
    letter = loop_variables[dimension-1]
    newpoint = centerpoint.replace(letter, letter+'-{}'.format(abs(myrange)))
    return newpoint

def distance_from_center(point='a[j][i]', loop_variables=None):
    """This function takes in input a point and the loop variables.
    It returns the distance of the point from the centerpoint
    """
    chars = ['a', 'b', '[', ']', '+', '-', ' ']
    chars += loop_variables

    for c in chars:
        point = point.replace(c, '')

    distance = 0
    for i in point:
        distance += int(i)

    return distance

def points_at_distance(points=None, loop_variables=None, distance=1):
    """This function takes in input a list of points, the loop variables and a
    distance. It returns a list containing all the points at the given distance
    """
    good_points = []
    for point in points:
        if distance_from_center(point, loop_variables) == distance:
            good_points.append(point)
        else:
            pass
    #print(good_points)
    return [x for x in good_points if x]



def boxpoint(centerpoint='a[j][i]', point=0, dimensions=2, radius=2,
             loop_variables=None):
    """This function takes in input a point (the centerpoint of the stencil),
    the dimensions of the grid, the radius of the stencil and the loop
    variables. It returns all the points of the box with chosen radius
    as a string
    """
    newpoint = centerpoint
    for i in range(dimensions):
        if i == 0:
            myvalue = point % (radius * 2 + 1)
        elif i == 1:
            myvalue = (point // (radius * 2 + 1)) % (radius * 2 + 1)
        elif i == 2:
            #elements_per_plane = (radius * 2 + 1)**2
            myvalue = point // ((radius * 2 + 1)**2)

        if myvalue < radius:
            sig = '-'
        elif myvalue > radius:
            sig = '+'
        else:
            sig = ''

        number = abs(int(myvalue - radius))
        if number == 0:
            strnumber = ''
        else:
            strnumber = str(number)


        newpoint = newpoint.replace(loop_variables[i], loop_variables[i]
                                    + '{}{}'.format(sig, strnumber))
    # print(str(point) + ' ' + newpoint)
    if newpoint == centerpoint:
        return None
    else:
        return newpoint
    #     point
    # if dimensions == 3:
    #     #to find out in which plan the point is:
    #     elements_per_plane = (radius * 2 + 1)**2
    #     z = point // elements_per_plane
    #     if z < radius:
    #         signum = '-'
    #     elif z > radius:
    #         signum = '+'
    #     else:
    #         signum = ''

    # x = point % (radius * 2 + 1)

    # y = (point // (radius * 2 + 1)) % (radius * 2 + 1)

    # if x < radius:
    #     signum = '-'
    # elif x > radius:
    #     signum = '+'
    # else:
    #     signum = ''

    # if y < radius:
    #     signum = '-'
    # elif y > radius:
    #     signum = '+'
    # else:
    #     signum = ''

    # xnumber = abs(int(point - radius))
    # if xnumber == 0:
    #     xvalue = ''
    # else:
    #     xvalue = str(xnumber)

    # ynumber = abs(int(point - radius))
    # if ynumber == 0:
    #     yvalue = ''
    # else:
    #     yvalue = str(ynumber)

    # newpoint = centerpoint
    # for letter in loop_variables:
    #     newpoint.replace(letter, letter+'{}{}'.format(signum, value))

def right(centerpoint='a[j][i]', dimension=2, loop_variables=None,
          myrange=1):
    """This function takes in input a point (the centerpoint of the stencil),
    the dimension (first, second or third) on which to act and the loop
    variables. It returns the point on its right, in the chosen dimension,
    as a string
    """
    letter = loop_variables[dimension-1]
    newpoint = centerpoint.replace(letter, letter+'+{}'.format(abs(myrange)))
    return newpoint


class BoxConstant(object):
    """class for the stencil with constant coefficients
    It can be:
        symmetric->isotropic
        symmetric->anisotropic
        asymmetric->anisotropic
    """

    name = "box"

    @classmethod
    def configure_arggroup(cls, parser):
        pass

    def __init__(self, dimensions=2, radius=1, symmetricity='symmetryx',
                 isotropy=True, datatype='double', inputgrids=1, args=None,
                 parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to
            2 (2 dimensional stencil)
        *radius* represents the radius of the stencil on each side of each
            dimension
        *symmetricity* is a boolean representing the symmetricity of the stencil
            with respect to the coefficients
        *isotropy* is a boolean representing the isotropy of the stencil (no
            dependency on the direction)
        *coeff* represents the coefficients of the stencil: can be either
            constant or variable.
        *datatype* represents the type of the data to be store in the grids. By
            default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions

        self.dims = []
        #save the letter of the dimensions in a variable (if 2 dimensions they
        #are "M" and "N")
        myascii = string.ascii_uppercase.replace("O", "")
        for i in range(0, self.dimensions):
            self.dims.append(myascii[12+i])

        self.radius = radius
        self.symmetricity = symmetricity
        self.isotropy = isotropy
        self.datatype = datatype

        self.inputgrids = inputgrids
        #to be changed in future to allow stencil on more than 2 grids
        self.inputs = [string.ascii_lowercase[i] for i in range(inputgrids)]
        self.output = string.ascii_lowercase[inputgrids]


        if self.isotropy and self.symmetricity == 'symmetric':#iso-sym
            self.num_coefficients = (radius * self.dimensions) + 1
        elif self.symmetricity == 'symmetric' and not self.isotropy:#aniso-sym
            self.num_coefficients = (2 * radius * self.dimensions) + 1
        elif self.symmetricity == 'homogeneous':
            self.num_coefficients = 1
        else:#not symmetricity and not isotropy)
            self.num_coefficients = (self.radius * 2 + 1)**self.dimensions

        # myformat = '{}' * num_coefficients
        # myletter = string.ascii_lowercase[inputgrids+1]
        # self.coefficients = [
        # myletter+str(i) for i in range(self.num_coefficients)
        # ]
        myletter = string.ascii_lowercase[inputgrids+1]
        self.coefficients = [
            myletter+str(i) for i in range(self.num_coefficients)
        ]


        #get the indices saved in a variable
        #(if 2 dimension they are "j" and "i")
        self.loop_variables = []
        for i in reversed(range(0, self.dimensions)):
            self.loop_variables.append(string.ascii_lowercase[8+i])

        self._args = args
        self._parser = parser

        if args:
            # handle CLI info
            pass

    def declaration(self):
        """
        builds the declaration of thevariables for the C code and returns
        declaration as unicode
        """
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

        #declare a variable for the initialization of the coefficients.
        #it is a string
        init_coefficients = ''
        #we declare all the coefficients line by line (can be modified to make
        #only a 1 line)
        for i in self.coefficients:
            init_coefficients = init_coefficients + self.datatype + ' ' + i \
            + ';\n'

        declaration = ''
        for i in init_inputs:
            declaration = declaration + i + '\n'

        declaration = declaration + init_output + '\n' + init_coefficients \
                    + '\n'

        return declaration



    def loop(self):
        '''
        builds the loop for the C code and returns loop_lines as list
        '''
        loop_lines = []

        #build the lines of the foor loop, according to the dimensions we have
        for i in range(0, self.dimensions):
            line = 'for(int {}={}; {} < {}-{}; {}++)'.format(
                self.loop_variables[i],
                self.radius,
                self.loop_variables[i],
                self.dims[i], self.radius,
                self.loop_variables[i]
                ) + '{'
            loop_lines.insert(i, line)

        centerpoint = self.inputs[0]
        lefthand = self.output

        # build the centerpoint and the lefthand of the equation
        for i in range(0, self.dimensions):
            lefthand = lefthand + '[' + self.loop_variables[i] + ']'
            centerpoint = centerpoint + '[' + self.loop_variables[i] + ']'

        # declare an empty stencil line (string)
        stencil = ''
        max_distance = self.dimensions * self.radius

        points = []
        for i in range((self.radius * 2 + 1)**self.dimensions):
            mypoint = boxpoint(centerpoint, i, self.dimensions, self.radius,
                               self.loop_variables)
            if mypoint:
                points.append(mypoint)

        ordered_points = []
        for i in range(1, max_distance+1):
            ordered_points.insert(i, points_at_distance(points,
                                                        self.loop_variables, i))

        if self.symmetricity == 'symmetric' and self.isotropy:
            print("symmetric and isotropic")

            stencil = self.coefficients[0] + ' *' + centerpoint + '\n'
            count = 1

            for i in range(max_distance):

                stencil = stencil + '+ {} * ({})\n'.format(
                    self.coefficients[count], ' + '.join(ordered_points[i]))

                count += 1

        elif self.symmetricity == 'asymmetric':
            print("asymmetric")

            stencil = self.coefficients[0] + ' * ' + centerpoint + '\n'
            count = 1

            for point in points:
                stencil = stencil + '+ {} * {}\n'.format(
                    self.coefficients[count], point)
                count += 1


        # anisotropic and symmetric
        elif not self.isotropy and self.symmetricity == 'symmetric':
            print("symmetric and anisotropic")
            print(ordered_points)

            stencil = self.coefficients[0] + ' *' + centerpoint + '\n'
            count = 1

            for i in range(max_distance):

                stencil = stencil + '+ {} * ({})\n'.format(
                    self.coefficients[count], ' + '.join(ordered_points[i]))

                count += 1

        elif self.symmetricity == 'homogeneous':
            print("homogeneous")

            stencil = self.coefficients[0] + ' * (' + centerpoint + '\n'

            for point in points:
                stencil = stencil + '+ {}'.format(point) + '\n'

            stencil = stencil + ')'


        righthand = '{};'.format(stencil)

        closing = '}\n' * self.dimensions

        computation = lefthand + ' = ' + righthand + '\n' + closing

        loop_lines.append(computation)

        return loop_lines



class BoxVariable(object):
    """class for the stencil with variable coefficients
    It can be:
        symmetric->isotropic
        symmetric->anisotropic
        asymmetric->anisotropic
    """
    name = "boxVariable"

    @classmethod
    def configure_arggroup(cls, parser):
        pass

    def __init__(self, dimensions=2, radius=1, symmetricity=True, isotropy=True,
                 datatype='double', inputgrids=1, args=None, parser=None):
        """
        *dimensions* is the number of dimensions of the stencil. It defaults to
            2 (2 dimensional stencil)
        *radius* represents the radius of the stencil on the max dimension
        *symmetricity* is a boolean representing the symmetricity of the
            stencil with respect to the coefficients
        *isotropy* is a boolean representing the isotropy of the stencil
            (no dependency on the direction)
        *coeff* represents the coefficients of the stencil: can be either
            constant or variable.
        *datatype* represents the type of the data to be store in the grids.
            By default double precision.
        *args* (optional) are the parsed arguments from the comand line
        """

        self.dimensions = dimensions

        self.dims = []
        #save the letter of the dimensions in a variable
        #(if 2 dimensions they are "M" and "N")
        for i in range(0, self.dimensions):
            self.dims.append(string.ascii_uppercase[12+i])

        self.radius = radius
        self.symmetricity = symmetricity
        self.isotropy = isotropy
        self.datatype = datatype

        self.inputgrids = inputgrids
        #to be changed in future to allow stencil on more than 2 grids
        self.inputs = [string.ascii_lowercase[i] for i in range(inputgrids)]
        self.output = string.ascii_lowercase[inputgrids]


        if self.isotropy and self.symmetricity == 'symmetric':#iso-sym
            self.num_coefficients = (radius * self.dimensions) + 1
        elif self.symmetricity == 'symmetric' and not self.isotropy:#aniso-sym
            self.num_coefficients = (2 * radius * self.dimensions) + 1
        elif self.symmetricity == 'homogeneous':
            self.num_coefficients = 1
        else:#not symmetricity and not isotropy)
            self.num_coefficients = (self.radius * 2 + 1)**self.dimensions

        self.coefficients = ['W']

        #get the indices saved in a variable
        #(if 2 dimension they are "j" and "i")
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
        builds the declaration of thevariables for the C code and returns
        declaration as unicode
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
        #input/output matrices
        for lineno in range(len(init_inputs)):
            for i in range(0, self.dimensions):
                line = init_inputs[lineno] + '[' + self.dims[i] + ']'
                init_inputs.pop(lineno)
                init_inputs.insert(lineno, line)
                init_output = init_output + '[' + self.dims[i] + ']'

        #coefficients matrix
        for lineno in range(len(init_coefficients)):
            for i in range(0, self.dimensions):
                line = init_coefficients[lineno] + '[' + self.dims[i] + ']'
                init_coefficients.pop(lineno)
                init_coefficients.insert(lineno, line)
            #add extra dimension for the weighting factor
            line = init_coefficients[lineno] + '[' + str(self.num_coefficients)\
            + ']'
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
        loop_lines = []

        #build the lines of the foor loop, according to the dimensions we have
        for i in range(0, self.dimensions):
            line = 'for(int {}={}; {} < {}-{}; {}++)'.format(
                self.loop_variables[i], self.radius, self.loop_variables[i],
                self.dims[i], self.radius, self.loop_variables[i]) + '{'
            loop_lines.insert(i, line)

        centerpoint = self.inputs[0]
        lefthand = self.output
        coefficients = self.coefficients

        # build centerpoint, lefthand and coefficients of the equation
        for i in range(0, self.dimensions):
            lefthand = lefthand + '[' + self.loop_variables[i] + ']'
            centerpoint = centerpoint + '[' + self.loop_variables[i] + ']'
            for coefficient in range(0, len(coefficients)):
                coeff = coefficients.pop(coefficient)
                coeff = coeff + '[' + self.loop_variables[i] + ']'
                coefficients.insert(coefficient, coeff)

        # declare an empty stencil line (string)
        stencil = ''

        max_distance = self.dimensions * self.radius

        points = []
        for i in range((self.radius * 2 + 1)**self.dimensions):
            mypoint = boxpoint(centerpoint, i, self.dimensions, self.radius,
                               self.loop_variables)
            if mypoint:
                points.append(mypoint)

        ordered_points = []
        for i in range(1, max_distance+1):
            ordered_points.insert(i, points_at_distance(points,
                                                        self.loop_variables, i))

        # isotropic and symmetric
        if self.symmetricity and self.isotropy:
            print("symmetric and isotropic")

            stencil = self.coefficients[0] + '[0]' + ' * ' + centerpoint + '\n'
            count = 1

            for i in range(max_distance):

                stencil = stencil + '+ {} * ({})\n'.format(
                    str(self.coefficients[0]) + '[' + str(count) + ']',
                    ' + '.join(ordered_points[i]))

                count += 1

        # asymmetric
        elif not self.symmetricity:
            print("asymmetric")

            stencil = self.coefficients[0] + '[0]' + ' * ' + centerpoint + '\n'
            count = 1

            for point in points:
                stencil = stencil + '+ {} * {}\n'.format(
                    str(self.coefficients[0]) + '[' + str(count) + ']', point)
                count += 1

        # anisotropic and symmetric
        elif not self.isotropy and self.symmetricity:
            print("symmetric and anisotropic")

            stencil = self.coefficients[0] + '[0]' + ' * ' + centerpoint + '\n'
            count = 1

            for i in range(max_distance):

                stencil = stencil + '+ {} * ({})\n'.format(
                    str(self.coefficients[0]) + '[' + str(count) + ']',
                    ' + '.join(ordered_points[i]))

                count += 1

        elif self.symmetricity == 'homogeneous':
            print("homogeneous")

            stencil = self.coefficients[0] + '[0]' + ' * ' + centerpoint + '\n'

            for point in points:
                stencil = stencil + '+ {}'.format(point) + '\n'

            stencil = stencil + ')'


        righthand = '({});'.format(stencil)

        closing = '}\n' * self.dimensions

        computation = lefthand + ' = ' + righthand + '\n' + closing

        loop_lines.append(computation)

        return loop_lines
