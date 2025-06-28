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

h-NUMO is a hydrostatic version of the Non-hydrostatic Model of the Ocean (NUMO). It is designed for large-scale ocean applications while maintaining computational efficiency. The model is implemented in the Galerkin Numerical Modeling Environment(GNuME) framework, which connects the work in [@NUMA] and the present model. The GNuME framework employs an arbitrary polynomial basis function expansion and provides a choice between continuous Galerkin (CG) and discontinuous Galerkin (DG) methods [@abdi2016efficient]. The framework was previously used to construct the Non-Hydrostatic Unified Model of the Ocean [@kopera2023non].

h-NUMO is an isopycnal model, meaning it approximates the ocean as a stack of layers, each with a constant density specified by the user. Other widely used vertical coordinates in models of the general ocean circulation include level (z) and terrain-fitted (\sigma) [@griffies2000developments,@chassignet2006generalized]. Some models, including the Hybrid Coordinate Ocean Model ([HYCOM](https://github.com/HYCOM)) and Modular Ocean Model version 6 ([MOM6](https://github.com/NOAA-GFDL/MOM6)), use a hybrid coordinate that employs all of these, with different coordinates in different regions.

Most ocean models employ finite difference or finite volume methods, but these methods have limitations in terms of accuracy and flexibility. Considering the various timescales of relevant ocean processes implies that ocean modeling is intrinsically multiscale, and accurately representing these physical processes presents computational challenges. The numerical methods used in ocean models should include desirable features such as low artificial dissipation, efficient resolution of localized flow features, and the capability to handle complex coastline geometries. High-order element-based discontinuous Galerkin (DG) methods promise to address those needs [@escobar2012high].

Additionally, DG methods can achieve higher accuracy with fewer degrees of freedom compared to finite difference or finite volume methods. The h-NUMO leverages these DG features, making it well-suited for large-scale applications. It is implemented in Fortran 90 and uses MPI for parallel computing, making it ideal for high-performance computing environments. The model is also extensible, allowing for the addition of new features and capabilities as needed.

Currently, h-NUMO is designed as a prototype for testing the applicability of high-order numerical methods in idealized and large-scale ocean applications. Explore how an increase in polynomial order approximation and a decrease in the number of elements in the computational domain affect the model's accuracy. Analyze the results of such changes and compare them with those from other ocean models.


# References