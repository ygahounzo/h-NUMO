---
title: 'h-NUMO: A high-performance discontinuous-Galerkin code for hydrostatic ocean modelling'
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

hNUMO is a high-order discontinuous-Galerkin (DG) isopycnal ocean circulation model that uses barotropic–baroclinic splitting to separate fast and slow motions into coupled subsystems. Extending the work of [@higdon2015multiple] to two horizontal dimensions, hNUMO applies this methodology to two horizontal dimensions using an arbitrary-order, nodal discontinuous Galerkin discretization. 

As an isopycnal model, hNUMO represents the ocean as a stack of constant-density layers. In contrast to traditional $z$-level and terrain-following ($\sigma$) coordinates [@griffies2000developments,@chassignet2006generalized], hybrid-coordinate models such as Hybrid Coordinate Ocean Model ([HYCOM](https://github.com/HYCOM)) and Modular Ocean Model version 6 ([MOM6](https://github.com/NOAA-GFDL/MOM6)) combine multiple vertical coordinates across different regions.

The model's accuracy and performance is assessed through a series of numerical experiments, including a double-gyre circulation test case, with results benchmarked against those from the HYCOM model. While hNUMO incurs a higher computational cost per degree of freedom, it attains comparable resolved kinetic energy in substantially less wall-clock time, achieving an order-of-magnitude speedup and requiring fewer computational resources. Additionally, hNUMO exhibits strong parallel scalability and sustains high efficiency even when run with a fewer number of grid elements per computational core. A detailed description of the model and its numerical discretization as well as model performance is provided in [@GAHOUNZO2026114496].

# Statement of need

Barotropic–baroclinic splitting is commonly used in layered ocean circulation models, which are typically discretized using finite-difference or finite-volume methods. Although these methods are effective, they exhibit limitations in both accuracy and geometric flexibility. Specifically, implementing high-order finite-volume schemes on unstructured meshes is challenging because they require wide reconstruction stencils.

High-order element-based discontinuous Galerkin (DG) methods represent a promising alternative, as they deliver high-order accuracy, reduced numerical dissipation, and efficient resolution of localized flow features [@escobar2012high]. Nevertheless, the application of DG methods to barotropic–baroclinic splitting in ocean circulation models has been limited. The hNUMO implements an arbitrary-order nodal discontinuous Galerkin discretization for both split subsystems. This approach provides a flexible and accurate framework for ocean circulation modeling, facilitating high-order simulations using non-uniform grids while maintaining the efficiency of barotropic–baroclinic coupling.

# h-NUMO features

hNUMO is implemented in the Galerkin Numerical Modeling Environment (GNUME) framework, which supports arbitrary-order polynomial bases and both continuous and discontinuous Galerkin methods [@abdi2016efficient]. The framework was previously used to develop the Non-Hydrostatic Unified Model of the Ocean [@kopera2023non]. The model is written in Fortran 90 and parallelized with MPI for high-performance computing. It is simple and easy to use for idealized and large-scale ocean applications for high-resolution and accurate results.

h-NUMO solves a stack of shallow water equations and uses a time-stepping scheme that splits the governing equations into barotropic (2D) and baroclinic (3D) components. The barotropic equations are solved using strong stability-preserving Runge-Kutta (SSPRK) schemes [@ruuth2006global], while the baroclinic equations are handled with a predictor-corrector method [@higdon2015multiple]. For the barotropic solver, the user can select RK2, SSPRK33, SSPRK34, or SSPRK35 by setting the number of stages to 2, 3, 4, or 5 in the input files. In addition, the user can specify an arbitrary polynomial order, which offers computational efficiency gains since increasing the approximation order is often more effective than refining the mesh.

h-NUMO provides two boundary condition types: free-slip and no-slip. The user is provided with the option to use linear or quadratic bottom drag and specify vertical diffusivity and viscosity. The model is designed to be extensible, providing the flexibility to add new features and capabilities as needed. The model also suport external grid mesh, generated with Gmsh [@geuzaine2009gmsh] and it can also run using unstructured meshes.

# Acknowledgements

This code was based on Prof. Robert L. Higdon's 2D (one horizantal and a vertical dimension) code, for which we are grateful to him for sharing it. We are thankful to Prof. Frank X. Giraldo for reviewing this code. The Office of Naval Research supports this work through Grant N00014-20-1-2038.

# References
