.. role:: bash(code)
   :language: bash

h-NUMO Directory Structure
*******************************

A typical copy of NUMA consists of the following:

| :bash:`Makefile`
| :bash:`config.user`
| :bash:`config.p4est`
| :bash:`INSTALL`
| :bash:`docs`
| :bash:`Graphics`
| :bash:`input_files`
| :bash:`lib`
| :bash:`make_depend.pl`
| :bash:`output`
| :bash:`run_scripts`
| :bash:`src`
| :bash:`p4est`
| :bash:`unused`

bin
-------

Upon successfully compiling h-NUMO, the ``bin`` directory will be created along with the executable (``numo3d``).

docs
---------

This directory contains the html files that form the documentation that you are now reading.

Graphics
----------

This directory contains a collection of Matlab m-files that can be used to plot two-dimensional figures (e.g., contours, slices, etc.). Many of these m-files are denoted with a specific case number that allows you to use the m-files without worrying too much about the test case being run.

include
-----------

This directory contains the ``*.mod`` files that are created upon compilation.

input_files
---------------

This directory contains a collection of the input files that may be used to run h-NUMO with specific case numbers. All of the files in this directory should work exactly as they are. Simply copy one of these files to the directory you wish to run from and rename the file to ``numo3d.in``.

lib
-----

This directory contains LAPACK libraries compiled from VERSION 3.11.0: November 2022.

output
--------

This directory contains a list of outputs for various test cases. It is sometimes helpful to see what the outputs of each test case should be.

run_scripts
-------------

This directory contains some simple run scripts for use with h-NUMO. It should give the user an idea of how to generate BASH scripts to run with h-NUMO.

src
------

This is the heart of the code and contains all of the source files actively used in h-NUMO. To add a new file to h-NUMO, you have to include it in the file ``Makefile``.

p4est
------

h-NUMO uses the p4est library :cite:p:`BursteddeWilcoxGhattas11` for the data structures and algorithms for parallel mesh generation, partitioning, and load balancing.

