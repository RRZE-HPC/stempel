stempel
=========

Stencil TEMPlate Engineering Library

This tool allows automatic generation of loop kernels for several kind of stencil patterns.
Stempel provides a framework to generate, perform static analysis of the code (exploiting kerncraft),
generate a full C code and parallelize it via openMP and blocking.

For a detailed documentation see in `<docs/>`_.

Installation
============

On most systems with python pip and setuputils installed, just run:
``pip install --user stempel`` for the latest release.

If you want to build from source:
Clone this repository and run ``python ./setup.py install``.

If you are unfamiliar with python, here is a tutorial on how to install python packages: https://packaging.python.org/installing/ . The use of virtual enviornments is usually a good choice.

Additional requirements are:
 * `kerncraft <https://github.com/RRZE-HPC/kerncraft>`_ (used to generate performance models)
 * `Intel Achitecture Code Analyzer (IACA) <https://software.intel.com/en-us/articles/intel-architecture-code-analyzer>`_, with (working) ``iaca.sh`` in PATH environment variable (used by ECM, ECMCPU and RooflineIACA models)
 * `likwid <https://github.com/RRZE-HPC/likwid>`_ (used in Benchmark model)

Usage
=====

1. Run stempel C code generator

``stempel gen -D 2 -r 1 -s -i``

2. Run stempel benchmark generator

``stempel bench code.c -m machine_file.yaml -b 32``


Credits
=======

Implementation: Danilo Guerrera

kerncraft : Julian Hammer

ECM Model (theory): Georg Hager, Holger Stengel, Jan Treibig

License
=======
AGPLv3
