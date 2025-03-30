.. role:: bash(code)
   :language: bash

Running h-NUMO
***************

Setting up the Workspace
==========================

Once you have successfully compiled h-NUMO, we recommend the following steps to running h-NUMO.

	1.  In the home/root directory create a directory called ``tests``. Inside ``tests`` create another directory for a specific simulation (for example, let's call it ``bump``.
	2.  Inside of bump copy ``numa3d.in.bump`` from input_files and rename it ``numa3d.in``.
	3.  Inside of bump copy one of the ``run_numo3d`` from ``run_scripts``. Make sure that the other batch submission commands make sense.

Input file: numo3d.in
-------------------

.. literalinclude:: numo3d.in
	:encoding: latin-1
	:linenos:

.. note::
	Check the sections below for more details.

gridnl
-------

.. literalinclude:: gridnl.in
   :linenos:

input
-------

.. literalinclude:: input.in
   :linenos:


h-NUMO Output
=================


    1.  Inside of ``bump`` directory run h-NUMO using your run script.

    	* :bash:`sbatch run_script.sh`

    2.  At the end of the run, you will have a collection of output files. In the input 	file, you can dump the simulation data as a ``vtk`` file for Paraview visualizations or simply as a ``txt`` file for Matlab visualizations. For the Matlab visualization, turn matlab_viz to true; otherwise, it will output ``vtk`` files.
