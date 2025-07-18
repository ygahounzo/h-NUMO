.. h-NUMO documentation master file, created by
   sphinx-quickstart on Wed Apr 10 23:15:02 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

h-NUMO
==================================

**h-NUMO** is an isopycnal ocean circulation model developed  in the Galerkin Numerical Modeling Environment(GNuME) framework (Giraldo, 2016). The GNuME framework uses an arbitrary polynomial basis function expansion and offers a choice of continuous Galerkin (CG) and discontinuous Galerkin (DG) methods Abdi and Giraldo (2016). The framework was previously used to construct the Non-Hydrostatic Unified Model of the Ocean Kopera et al. (2023). The h-NUMO model use arbitrary number of layers. We are thankful to Prof. Robert L. Higdon for sharing his 2D code, which h-NUMO is based on.

.. note::

   This project is under active development.

h-NUMO is

- `Easy to install and configure <installation.html>`_ on a compute node.

- Easy to learn with [online documentation](https://ygahounzo.github.io/h-NUMO/index.html), including a complete description of `the physisc and the numerics <numo_model.html>`_.

- Verified with a suites of `test cases <test.html>`_.

.. toctree::
   :numbered: 3
   :maxdepth: 2
   :caption: Contents

   numo_model
   installation
   directory
   running
   test
   bc_time_method
   papers
   ref



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
