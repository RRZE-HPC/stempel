'''
High-level tests for the overall functionallity and things in stempel.py
'''
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import division

import sys
import os
import unittest
import tempfile
import shutil
import pickle
from pprint import pprint
from io import StringIO
from distutils.spawn import find_executable
import platform

import six
import sympy

sys.path.insert(0, '..')
print(os.getcwd())
from stempel import stempel


class TestStempel(unittest.TestCase):
    def setUp(self):
        # Create a temporary directory
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.temp_dir)

    def _find_file(self, name):
        testdir = os.path.dirname(__file__)
        name = os.path.join(testdir, 'test_files', name)
        assert os.path.exists(name)
        return name

    def test_2d5pt_star(self):
        store_file = os.path.join(self.temp_dir, 'test_2d5pt_star')

        parser = stempel.create_parser()
        args = parser.parse_args(['-D', '2', '-r', '1', '-s', '-i', '-k',
        						  'star', '-C', 'constant', '--store',
        						  store_file])
        stempel.check_arguments(args, parser)
        stempel.run(args)


        with open(store_file + '.c', 'rb') as f:
       		results = f.read()

       	with open(self._find_file('2d-5pt.c'), 'rb') as f:
       		test_code = f.read()

       	six.assertCountEqual(self, results, test_code)
        # Check for correct variations of constants
        # print((sympy.var('M'))
        # six.assertCountEqual(self,
        #     [sorted(map(str, r)) for r in results],
        #     [sorted(map(str, r)) for r in [
        #         (sympy.var('M'), sympy.var('N')),
        #         (sympy.var('M'), sympy.var('N'))]])


        # Output of first result:
        # result = results['2d-5pt.c'][[k for k in results['2d-5pt.c']
        #                               if (sympy.var('N'), 1000) in k][0]]

        # six.assertCountEqual(self, result, ['ECMData'])

        # ecmd = result['ECMData']
        # 2 arrays * 1000*50 doubles/array * 8 Bytes/double = 781kB
        # -> fully cached in L3
        # self.assertAlmostEqual(ecmd['L1-L2'], 6, places=1)
        # self.assertAlmostEqual(ecmd['L2-L3'], 6, places=1)
        # self.assertAlmostEqual(ecmd['L3-MEM'], 0.0, places=0)