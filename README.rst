stempel
=========

Stencil TEMPlate Engineering Library

This tool allows automatic generation of loop kernels for several kind of stencil patterns.
Stempel provides a framework to generate, perform static analysis of the code (exploiting kerncraft),
generate a full C code and parallelize it via openMP and blocking.

For a detailed documentation see in `<docs/>`_.

Installation
============

Remember to setup properly your envirnment in advance. For example, using conda Python:

``conda create -n stempel python=3.4``

``source activate stempel``

Clone this repository and run:
``python ./setup.py install``


If you are unfamiliar with python, here is a tutorial on how to install python packages: https://packaging.python.org/installing/ . The use of virtual enviornments is usually a good choice.

Additional requirements are:
 * `kerncraft <https://github.com/RRZE-HPC/kerncraft>`_ (used to generate performance models)
 * `Intel Achitecture Code Analyzer (IACA) <https://software.intel.com/en-us/articles/intel-architecture-code-analyzer>`_, with (working) ``iaca.sh`` in PATH environment variable (used by ECM, ECMCPU and RooflineIACA models)
 * `likwid <https://github.com/RRZE-HPC/likwid>`_ (used in Benchmark model)

Usage
=====

1. Run stempel C code generator

``stempel gen -D 2 -r 1 -i``

or

``stempel gen -D 2 -r 2 -k box -C variable -p``

2. Run stempel benchmark generator

``stempel bench code.c -m machine_file.yaml -b 32``

3. Run a full analysis (a stencil is generated, analysed through kerncraft applying ECM/Data, Roofline and Layer Condition models; a project is set up in PROVA! and an experiment is executed. The outputs are stored to the STEMPEL workspace):

``analysis -w ~/Desktop/stempelwork -p ~/PROVA ~/Desktop/provastempel -k star -m BroadwellEP_E5-2697_CoD.yml -r 2 -d 2 -e 5 -t 2 4 8 10 --method_type OpenMP-4.0-GCC-4.9.3-2.25 -C constant -c isotropic -l /apps/likwid/system/include/ /apps/likwid/system/lib/ --iaca``


Credits
=======

Implementation: Danilo Guerrera

kerncraft : Julian Hammer

ECM Model (theory): Georg Hager, Holger Stengel, Jan Treibig

License
=======
AGPLv3
