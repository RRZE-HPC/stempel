'''Collection of stencils
The exported classes must have the following class level attributes:
  * name (str) is the name of the stencil kind (no abreviatation)
  * configure_subparser(parser) classmethod that configures the parser for cli usage
  * declaration()
  * loop()
'''
from .starstencil import StarConstant
from .starstencil import StarVariable
from .boxstencil import Box


__all__ = ['StarConstant', 'StarVariable', 'Box']