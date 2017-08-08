'''Collection of stencils
	The exported classes must have the following class level attributes:
  * name (str) is the name of the stencil kind (no abreviatation)
  * declaration()
  * loop()
'''
from .starstencil import StarConstant
from .starstencil import StarVariable
from .boxstencil import BoxConstant
from .boxstencil import BoxVariable


__all__ = ['StarConstant', 'StarVariable', 'BoxConstant', 'BoxVariable']