#!/usr/bin/env python
"""This module defines some utility functions used in the Stencil classes.
"""
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