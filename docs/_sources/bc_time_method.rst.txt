Boundary conditions and time intergration methods
*****************************************************

Boundary Conditions
---------------------

NUMA reads in boundary conditions via the file numa3d.in through the 
variables x_boundary, y_boundary, and z_boundary. NUMA allows for the following 
types of boundary conditions:

* 0 = do nothing boundary conditions
* 2 = no-slip boundary condition
* 4 = free-slip boundary condition


1. Do Nothing Boundary Condition
	Setting z_boundary=(0,0) tells h-NUMO to ignore these boundaries.

2. No-slip boundary condition
	This boundary condition cancels all direction velocity at a boundary (:math:`\mathbf{u} = 0`). 

3. No-slip or No-flux Boundary Condition
	This boundary condition cancels the normal velocity at a boundary (:math:`\mathbf{n}\cdot\mathbf{u} = 0`). 

Time-Integrators
----------------------------------------------------

h-NUMO is equipped with a suite of time-integrators which include:

* Explicit 5-stage 3rd Order Runge-Kutta Methods (used only for the barotropic equations),
* Two-level or predictor-corrector method for both barotropic and baroclinic equations.
