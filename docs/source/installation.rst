.. role:: bash(code)
   :language: bash

Installation
================

Libraries
---------------

To compile the h-NUMO model, one needs the following modules

| :bash:`gcc version 8 and above`
| :bash:`openmpi`
| :bash:`make`
| :bash:`lapack`
| :bash:`netcdf`

Installation instructions
--------------------------------

h-NUMO will need the :bash:`mpi90` command in order to compile. It is up to the user to ensure that these packages exist on his computing platform.

1. Clone the repository to a local directory

   * ``git clone git@github.com:ygahounzo/h-NUMO.git``

2. Go to the base directory of the repository

   * cd ``h-NUMO``

2. Install and compile h-NUMO

   Inside the ``h-NUMO`` directory type the following command to compile p4est and h-NUMO

   ``make hnumo``



   .. warning::

      Sometimes simple make does not work, and you need to compile p4est 

      * Compile p4est by typing from ``h-NUMO`` root directory
         ``make p4est/local/lib/libp4est.a``

      * Export path to the p4est library
         To make sure this is exported every time you log-in, put the following line to :bash:`~/.bash_profile` file:

         ``export LD_LIBRARY_PATH=$PATH:<path_to_h-NUMO>/p4est/local/lib``

         where <path_to_h-NUMO> is the location where you have installed h-NUMO

         Just in case type source ~./bash_profile, and it is not a bad idea to include this source command in your run scripts

      * Compile NUMO by going to its directory and typing:
         ``make hnumo`` (check the config.user if you want to compile on a specific platforms)


Once you successfully make h-NUMO, additional directories will result. These are:


         *  :bash:`bin`

         *  :bash:`depend.mk`

         *  :bash:`include`





