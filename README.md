# h-NUMO
Hydrostatic version of the Non-hydrostatic Model of the Ocean

We implement h-NUMO in the Galerkin Numerical Modeling Environment(GNuME) (Giraldo, 2016) framework. The GNuME framework uses an arbitrary polynomial basis function expansion and offers a choice of continuous Galerkin (CG) and discontinuous Galerkin (DG) methods Abdi and Giraldo (2016). The framework was previously used to construct the Non-Hydrostatic Unified Model of the Ocean Kopera et al. (2023), and we call the current implementation of multilayer shallow water equations h-NUMO to signify the hydrostatic aspect of the
model.

h-NUMO is 

- [Easy to install](https://ygahounzo.github.io/h-NUMO/installation.html) either on laptop or HPC systems.

- Easy to configure: all the parameters are specified in a single [input file](https://ygahounzo.github.io/h-NUMO/running.html). The user can also define their own initial conditions.

- Easy to learn with an extensive [online documentation](https://ygahounzo.github.io/h-NUMO/index.html).

- Tested on idealized [test cases](https://ygahounzo.github.io/h-NUMO/test.html).

- Verified: The model was succesfully compared with HYCOM model.