#!/usr/bin/env python
"""This module defines some utility functions used in the Stencil classes.
"""
import re

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

def remove_for(loop_code):
    mycode = loop_code.split('\n')
    index = []
    for i in range(len(mycode)):
        if mycode[i].startswith('for'):
            index.append(i)
    line = index[-1]

    mycode = mycode[line+1:]
    mycode = '\n'.join(mycode)

    return mycode


def count_ops(kernel_code):
    """
    Count the number of operations a kernel consists of.
    It accepts a string as input.
    The code must be only the kernel, with the for(s) 
    """
    mycode = remove_for(kernel_code)
    
    mycode = re.sub("[\[].*?[\]]", "", mycode)

    operations =['+', '-', '*', '/']
    count = 0
    for op in operations:
        count += mycode.count(op)

    return count