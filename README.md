# h-NUMO
[![tests](https://github.com/ygahounzo/h-NUMO/actions/workflows/ci.yml/badge.svg)](https://github.com/ygahounzo/h-NUMO/actions/workflows/ci.yml)
[![Documentation Status](https://github.com/ygahounzo/h-NUMO/actions/workflows/docs.yml/badge.svg)](https://github.com/ygahounzo/h-NUMO/actions/workflows/docs.yml)

Hydrostatic version of the Non-hydrostatic Model of the Ocean

We implement h-NUMO in the Galerkin Numerical Modeling Environment (GNuME) framework (Giraldo, 2016). The GNuME framework uses an arbitrary polynomial basis function expansion and provides a choice between continuous Galerkin (CG) and discontinuous Galerkin (DG) methods (Abdi and Giraldo, 2016). The framework was previously used to construct the Non-Hydrostatic Unified Model of the Ocean (Kopera et al., 2023), and we refer to the current implementation of multilayer shallow water equations as h-NUMO to signify the hydrostatic aspect of the model. We are thankful to Prof. Robert L. Higdon for sharing his 2D code, which h-NUMO is based on.

h-NUMO is 

- Easy to [install](https://ygahounzo.github.io/h-NUMO/installation.html) either on laptop or HPC systems.

- Simple to configure: all the parameters are specified in a single [input file](https://ygahounzo.github.io/h-NUMO/running.html). The user can also define their own initial conditions.

- Easy to learn with an extensive [online documentation](https://ygahounzo.github.io/h-NUMO/index.html).

- Tested on idealized [test cases](https://ygahounzo.github.io/h-NUMO/test.html).

- Verified: The model was succesfully compared with HYCOM model.
