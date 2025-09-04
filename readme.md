## RBF Incompressible Navier-Stokes

### Detailed Description

This MATLAB code implements a mesh-free method to solve the incompressible Navier-Stokes equations for flow around a cylinder using Polyharmonic Spline Radial Basis Functions (PHS-RBF). The approach uses a staggered grid formulation and a fractional step method for time integration.

### How to Run the Code

1. **Requirements**:
   - MATLAB (the code has been tested on MATLAB R2019b and newer versions)
   - No additional toolboxes are required

2. **Running the Simulation**:
   - Open MATLAB and navigate to the repository directory
   - Run the main script by typing `cylinder` in the MATLAB command window
   - Alternatively, open `cylinder.m` in the MATLAB editor and click the "Run" button

3. **What to Expect**:
   - The script will generate a mesh for the flow domain
   - Set up the differentiation matrices using RBF-FD
   - Run the time-stepping loop for 5000 iterations (by default)
   - Display the final velocity field as scatter plots showing the u and v components

4. **Modifying Simulation Parameters**:
   - Reynolds number: Look for the variable `nu` in `cylinder.m` (currently set to 1/100)
   - Time step: Look for the variable `dt` in `cylinder.m` (currently set to 1e-2)
   - Number of time steps: Look for the variable `Nt` in `cylinder.m` (currently set to 5000)
   - These parameters are defined in the simulation setup section of the code

5. **Computational Requirements**:
   - The simulation is computationally intensive and may take time to complete
   - For a full simulation with 5000 time steps, expect several minutes of computation time depending on your hardware

#### Key Components

##### 1. Domain and Mesh Generation
- The computational domain is a rectangle ($-8 \leq x \leq 24$, $-8 \leq y \leq 8$) with a cylinder of radius 0.5 at the origin
- The mesh is generated using the DistMesh library (included in the distmesh/ folder)
- Local grid refinement is applied near the cylinder and in the wake region
- The code separates nodes into:
  - Interior nodes (`xy_s` for pressure, `xy` for velocity)
  - Boundary nodes (inlet, outlet, top/bottom walls, and cylinder surface)

##### 2. RBF-FD Discretization
- `RBF_PHS_FD_all.m` constructs differentiation matrices using PHS-RBF with polynomial augmentation
- For each node, a local stencil of k=35 nearest neighbors is used
- The RBF interpolant is constructed as:
  - PHS basis functions: $\phi(r) = r^m$ (where m=3 typically)
  - Polynomial terms up to degree d=3 (or higher in some regions)
- The method computes weights for various differential operators:
  - First derivatives ($\partial/\partial x$, $\partial/\partial y$)
  - Laplacian ($\nabla^2$)
  - Second derivatives ($\partial^2/\partial x^2$, $\partial^2/\partial y^2$, $\partial^2/\partial x \partial y$)

##### 3. Staggered Grid Approach
- Pressure values are stored at nodes `xy_s` (P-grid)
- Velocity components are stored at nodes `xy` (V-grid)
- Differentiation operators are constructed to map between grids:
  - `D0_{12}_x`, `D0_{12}_y`: Derivatives from V-grid to P-grid
  - `D0_{21}_x`, `D0_{21}_y`: Derivatives from P-grid to V-grid

##### 4. Boundary Conditions
- Inlet: Uniform flow (U=1, V=0)
- Outlet: Neumann boundary conditions
- Top/bottom walls: No-slip conditions
- Cylinder surface: No-slip conditions
- Special differentiation matrices are constructed for each boundary

##### 5. Time Integration (Fractional Step Method)
The `NS_2d_cylinder_PHS.m` function implements a fractional step method:

1. **Advection Step**:
   - Uses Adams-Bashforth scheme ($\frac{3}{2}$ current - $\frac{1}{2}$ previous)
   - Computes non-linear terms: $H_U = -U \cdot \nabla U$, $H_V = -V \cdot \nabla V$

2. **Boundary Conditions**:
   - Enforces velocity boundary conditions after advection

3. **Viscous Step**:
   - Uses Crank-Nicolson scheme for diffusion terms
   - Implicit treatment: $(I - \frac{dt \cdot \nu \cdot \nabla^2}{2})U = U^* + \frac{dt \cdot \nu \cdot \nabla^2 U}{2}$
   - Solved using precomputed LU decomposition

4. **Pressure Correction**:
   - Computes divergence of intermediate velocity field
   - Solves Poisson equation for pressure: $\nabla^2 p = \frac{\nabla \cdot U^*}{dt}$
   - Uses precomputed inverse operator L_inv_s

5. **Projection Step**:
   - Projects velocity field to be divergence-free: $U = U^* - \nabla p$
   - Enforces final boundary conditions

##### 6. Simulation Parameters
- Reynolds number: $Re = 100$ ($\nu = 0.01$)
- Time step: $dt = 0.01$
- Number of time steps: $5000$

##### 7. Visualization
- The code plots the u and v velocity components at the final time step
- Uses scatter plots colored by velocity values

#### Mathematical Foundation

The incompressible Navier-Stokes equations being solved are:

1. Momentum equation: $\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \nu \nabla^2 \mathbf{u}$
2. Continuity equation: $\nabla \cdot \mathbf{u} = 0$

The fractional step method splits these equations into:
1. $\frac{\partial \mathbf{u}^*}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = \nu \nabla^2 \mathbf{u}$
2. $\mathbf{u}^{n+1} = \mathbf{u}^* - \nabla p$
3. $\nabla \cdot \mathbf{u}^{n+1} = 0$

Which leads to the pressure Poisson equation: $\nabla^2 p = \frac{\nabla \cdot \mathbf{u}^*}{dt}$

#### Notable Implementation Details
- The code uses LU decomposition for efficient solution of linear systems
- Higher-order RBF-FD stencils (more polynomial terms) are used near the cylinder
- The implementation handles the pressure-velocity coupling through the projection method
- The cylinder simulation should capture vortex shedding (von Kármán vortex street) at $Re=100$

### References
- T. Chu, O. T. Schmidt, RBF-FD discretization of the Navier-Stokes equations on scattered but staggered nodes, Journal of Computational Physics 474, 111756, 2023
- T. Chu, O. T. Schmidt, Mesh-free hydrodynamic stability, Journal of Computational Physics, Volume 478, 112268, 2023

### License (DistMesh)

DistMesh is a collection of MATLAB functions for generation and
manipulation of unstructured meshes. DistMesh is Copyright (C) 2004-2012
Per-Olof Persson.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

If you use DistMesh in any program or publication, please acknowledge
its authors by adding a reference to: Per-Olof Persson and Gilbert
Strang, "A Simple Mesh Generator in MATLAB," SIAM Review Vol. 46 (2)
2004.