---
title: 'h-NUMO: A Hydrostatic Version of the Non-hydrostatic Model of the Ocean'
tags:
  - isopycnal model
  - Ocean modeling
  - Discontinuous Galerkin
  - general circulation models
authors:
 - name: Yao Gahounzo
   orcid:
   affiliation: 1
 - name: Michal Kopera
   affiliation: 2
affiliations:
 - name: Center for Ocean and Atmospheric Studies, Florida State University, FL, USA
   index: 1
 - name: Department of Mathematics, Boise State University, ID, USA
   index: 2
date: 21 June 2025
bibliography: paper.bib
---

# Summary

h-NUMO is a hydrostatic version of the Non-hydrostatic Model of the Ocean (NUMO). It is designed for large-scale ocean applications while maintaining computational efficiency. The model is implemented in the Galerkin Numerical Modeling Environment(GNuME) framework, which connects the work in [@NUMA] and the present model. The GNuME framework employs an arbitrary polynomial basis function expansion and provides a choice between continuous Galerkin (CG) and discontinuous Galerkin (DG) methods [@abdi2016efficient]. The framework was previously used to construct the Non-Hydrostatic Unified Model of the Ocean [@kopera2023non]. h-NUMO is implemented in Fortran 90 and uses MPI for parallel computing, making it ideal for high-performance computing environments.

h-NUMO is an isopycnal model, meaning it approximates the ocean as a stack of layers, each with a constant density specified by the user. Other widely used vertical coordinates in models of the general ocean circulation include level ($z$) and terrain-fitted ($\sigma$) [@griffies2000developments,@chassignet2006generalized]. Some models, including Hybrid Coordinate Ocean Model ([HYCOM](https://github.com/HYCOM)) and Modular Ocean Model version 6 ([MOM6](https://github.com/NOAA-GFDL/MOM6)), use a hybrid coordinate system that employs all of these, with different coordinates in different regions.

# Statement of need

Most ocean models employ finite difference or finite volume methods, but these approaches have limitations in terms of accuracy and flexibility. Considering the various timescales of relevant ocean processes, ocean modeling is intrinsically multiscale, and accurately representing these physical processes presents computational challenges. The numerical methods used in ocean models should include desirable features such as low artificial dissipation, efficient resolution of localized flow features, and the capability to handle complex coastline geometries.

High-order element-based discontinuous Galerkin (DG) methods promise to address these needs [@escobar2012high]. With the DG method, h-NUMO achieves higher accuracy with fewer degrees of freedom compared to finite difference or finite volume methods, making it well-suited for large-scale applications. Implementing high-order finite volume methods is particularly complicated on unstructured meshes due to the need for large reconstruction stencils. The DG method utilizes higher-order basis functions within each cell to address this issue, allowing the model to use unstructured meshes.

h-NUMO is designed for testing the applicability of high-order numerical methods in idealized and large-scale ocean applications. It enables exploration of how increasing the polynomial order of approximation and decreasing the number of elements in the computational domain affect the model's accuracy. This makes h-NUMO ideal for high-resolution and more accurate simulations of both idealized and complex ocean flows.

# h-NUMO features

h-NUMO solves a stack of shallow water equations and uses a time-stepping scheme that splits the governing equations into barotropic (2D) and baroclinic (3D) components. The barotropic equations are solved using strong stability-preserving Runge-Kutta (SSPRK) schemes [@ruuth2006global], while the baroclinic equations are handled with a predictor-corrector method [@higdon2015multiple]. For the barotropic solver, the user can select RK2, SSPRK33, SSPRK34, or SSPRK35 by setting the number of stages to 2, 3, 4, or 5 in the input files. In addition, the user can specify an arbitrary polynomial order, which offers computational efficiency gains since increasing the approximation order is often more effective than refining the mesh.

h-NUMO provides two boundary condition types: free-slip and no-slip. The user is provided with the option to use linear or quadratic bottom drag and specify vertical diffusivity and viscosity. The model is designed to be extensible, providing the flexibility to add new features and capabilities as needed.


# Acknowledgements

This code was based on Prof. Robert L. Higdon's 2D code, for which we are grateful to him for sharing it. We are thankful to Prof. Frank X. Giraldo for reviewing this code. The Office of Naval Research supports this work through Grant N00014-20-1-2038.

# References
