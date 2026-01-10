.. role:: bash(code)
   :language: bash

h-NUMO Directory Structure
*******************************

A typical copy of h-NUMO consists of the following:

| :bash:`Makefile`
| :bash:`config.user`
| :bash:`config.p4est`
| :bash:`INSTALL`
| :bash:`docs`
| :bash:`Examples`
| :bash:`lib`
| :bash:`make_depend.pl`
| :bash:`src`
| :bash:`p4est`
| :bash:`unused`

bin
-------

Upon successfully compiling h-NUMO, the ``bin`` directory will be created along with the executable (``numo3d``).

docs
---------

This directory contains the html files that form the documentation that you are now reading.

Examples
----------
This directory contains several example test cases that can be run with h-NUMO. Each test case is contained in its own sub-directory, which includes input files, run scripts, a clear_data script to remove generated files, and plotting scripts.

include
-----------

This directory contains the ``*.mod`` files that are created upon compilation.

lib
-----

This directory contains the libraries that are used by h-NUMO.

src
------

This is the heart of the code and contains all of the source files actively used in h-NUMO. To add a new file to h-NUMO, you have to include it in the file ``Makefile``.

p4est
------

h-NUMO uses the p4est library :cite:p:`BursteddeWilcoxGhattas11` for the data structures and algorithms for parallel mesh generation, partitioning, and load balancing.

