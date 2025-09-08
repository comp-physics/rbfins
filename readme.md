## RBF Incompressible Navier-Stokes

### Detailed Description

This MATLAB code implements a mesh-free method to solve the incompressible Navier-Stokes equations for flow around obstacles using Polyharmonic Spline Radial Basis Functions (PHS-RBF). The approach uses a staggered grid formulation and a fractional step method for time integration.

**Supported Geometries**:
- **Cylinder**: Circular obstacles (original implementation)
- **Ellipse**: Elliptical obstacles with customizable semi-major and semi-minor axes
- **Rectangle**: Rectangular obstacles with customizable width, height, and center position

### How to Run the Code

1. **Requirements**:
   - MATLAB (the code has been tested on MATLAB R2019b and newer versions)
   - No additional toolboxes are required

2. **Running the Simulation**:
   - Open MATLAB and navigate to the repository directory
   - Run the main script by typing `simulate` in the MATLAB command window
   - Alternatively, open `simulate.m` in the MATLAB editor and click the "Run" button

3. **Selecting Geometry Types**:
   The code supports three geometry types, each with its own configuration:
   
   **Default (Cylinder)**:
   ```matlab
   simulate  % Runs cylinder geometry by default
   ```
   
   **Specific Geometry**:
   ```matlab
   % Set environment variables to specify geometry
   setenv('GEOMETRY_TYPE', 'cylinder');    % or 'ellipse' or 'rectangle'
   simulate
   ```
   
   **Direct Configuration**:
   ```matlab
   cfg = config('cylinder');    % Load cylinder configuration
   cfg = config('ellipse');     % Load ellipse configuration  
   cfg = config('rectangle');   % Load rectangle configuration
   ```

4. **What to Expect**:
   - The script will generate a mesh for the flow domain around the selected obstacle
   - Set up the differentiation matrices using RBF-FD
   - Run the time-stepping loop for the specified number of iterations (default: 5000)
   - Display the final velocity field as scatter plots showing the u and v components

5. **Modifying Simulation Parameters**:
   - **Unified Configuration**: All simulation parameters are handled by the `config(geometry_type)` function
   - **Key Parameters Available**:
     - Reynolds number: `config.simulation.reynolds_number` (default: 100)
     - Time step: `config.simulation.time_step` (default: 1e-2)
     - Number of time steps: `config.simulation.num_time_steps` (default: 5000)
     - Domain size: `config.domain.x_min/x_max/y_min/y_max`
     - **Geometry-Specific Parameters**:
       - **Cylinder**: `config.geometry.obstacle_radius` (default: 0.5)
       - **Ellipse**: `config.geometry.ellipse_a` (semi-major axis), `config.geometry.ellipse_b` (semi-minor axis)
       - **Rectangle**: `config.geometry.rect_width`, `config.geometry.rect_height`, `config.geometry.rect_x_center`, `config.geometry.rect_y_center`
     - RBF parameters: Various stencil sizes and polynomial degrees
   - **To modify parameters**: Edit the values in the `config.m` function for the specific geometry case
   - **Configuration Structure**:
     - `config.domain`: Domain boundaries and dimensions
     - `config.mesh`: Mesh generation parameters (automatically optimized per geometry)
     - `config.rbf`: RBF-FD algorithm parameters
     - `config.simulation`: Time stepping and physics parameters
     - `config.schemes`: Numerical scheme coefficients
     - `config.visualization`: Plot settings and visualization parameters

6. **Computational Requirements**:
   - The simulation is computationally intensive and may take time to complete
   - For a full simulation with 5000 time steps, expect several minutes of computation time depending on your hardware
   - Rectangle geometries require finer meshes and may take longer to converge

### Testing with GitHub Actions

This repository includes comprehensive automated testing using GitHub Actions to ensure the code runs correctly across different environments and geometries.

1. **Test Structure**:
   - **Geometry-Specific Tests**: Each geometry has its own test class (`TestCylinderGeometry.m`, `TestEllipseGeometry.m`, `TestRectangleGeometry.m`)
   - **Golden File Tests**: Regression tests that compare against known-good results (`TestGoldenCylinder.m`, `TestGoldenEllipse.m`, `TestGoldenRectangle.m`)
   - **Base Test Class**: `BaseGeometryTest.m` provides common functionality to avoid code duplication
   - **CI Environment**: Tests run with reduced time steps (20 instead of 5000) for faster execution

2. **CI Workflows**:
   - **MATLAB CI** (`.github/workflows/matlab.yml`): Runs all geometry tests and golden file comparisons
   - **MATLAB Code Quality** (`.github/workflows/ci-quality.yml`): Runs code analyzer (linting) and formatting checks

3. **What the Tests Verify**:
   - **Smoke Tests**: Ensure each geometry simulation runs without errors
   - **Configuration Validation**: Verify that geometry-specific config parameters are present
   - **Golden File Comparison**: Ensure numerical reproducibility by comparing against reference results
   - **Parameter Validation**: Confirm that all required fields are present in configuration

4. **Running Tests Locally**:
   ```matlab
   % Run all tests
   addpath('tests');
   runtests('tests');
   
   % Run tests for specific geometry
   runtests('tests/TestCylinderGeometry.m');
   runtests('tests/TestEllipseGeometry.m');
   runtests('tests/TestRectangleGeometry.m');
   
   % Run golden file tests
   runtests('tests/TestGoldenCylinder.m');
   runtests('tests/TestGoldenEllipse.m');
   runtests('tests/TestGoldenRectangle.m');
   ```

5. **Golden File Testing**:
   - **Purpose**: Ensures numerical reproducibility and detects unintended changes in simulation behavior
   - **Coverage**: Separate golden files for each geometry type
   - **Generation**: Golden files are automatically maintained and updated when needed
   - **Comparison**: Tests compare final velocity fields within specified tolerances (5e-3 relative/absolute)
   - **Deterministic**: Uses fixed random seed (42) for reproducible mesh generation

6. **Test Architecture**:
   - **Inheritance-Based**: All geometry tests inherit from `BaseGeometryTest` to share common logic
   - **Parameterized**: Tests automatically adapt to different geometry types
   - **Environment Variables**: CI tests use `CONFIG_FUNC` and `GEOMETRY_TYPE` environment variables

7. **Adding New Tests**:
   - For new geometries: Create a new test class inheriting from `BaseGeometryTest`
   - Define `GEOMETRY_TYPE` and `EXPECTED_FIELDS` properties
   - Golden files are automatically generated and maintained
   - Follow the existing pattern in `TestCylinderGeometry.m`

8. **Code Quality Tools**:
   - **MATLAB Code Analyzer**: Integrated linting checks for code quality
   - **Unified Code Formatting**: 
     - **Local**: Run `./format.sh` or `matlab -batch "format_code"` to format all code
     - **CI**: Automatically checks if code is properly formatted on pull requests
     - **Idempotent**: Formatter only modifies files when changes are needed
   - **CI Integration**: Both linting and formatting checks run automatically on pull requests

#### Key Components

##### 1. Domain and Mesh Generation
- The computational domain is a rectangle ($-8 \leq x \leq 24$, $-8 \leq y \leq 8$) with an obstacle at the origin
- **Supported Obstacle Types**:
  - **Cylinder**: Circular obstacle with configurable radius (default: 0.5)
  - **Ellipse**: Elliptical obstacle with configurable semi-major and semi-minor axes
  - **Rectangle**: Rectangular obstacle with configurable width, height, and center position
- The mesh is generated using the DistMesh library (included in the `lib/distmesh/` folder)
- **Adaptive Mesh Refinement**: Automatically applied near obstacles and in wake regions
- **Geometry-Specific Optimization**: Mesh parameters are automatically tuned for each geometry type
- The code separates nodes into:
  - Interior nodes (`xy_s` for pressure, `xy` for velocity)
  - Boundary nodes (inlet, outlet, top/bottom walls, and obstacle surface)

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
- **Inlet**: Uniform flow (U=1, V=0)
- **Outlet**: Neumann boundary conditions
- **Top/bottom walls**: No-slip conditions
- **Obstacle surface**: No-slip conditions (applies to all geometry types)
- **Adaptive Implementation**: Special differentiation matrices are constructed for each boundary type
- **Geometry-Agnostic**: Boundary condition implementation works seamlessly across all supported geometries

##### 5. Time Integration (Fractional Step Method)
The `NS_2d_fractional_step_PHS.m` function implements a fractional step method:

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
- **Efficient Linear Solvers**: Uses LU decomposition for efficient solution of linear systems
- **Adaptive Stencils**: Higher-order RBF-FD stencils (more polynomial terms) are used near obstacles
- **Pressure-Velocity Coupling**: Handled through the projection method for incompressible flow
- **Vortex Shedding**: All geometries can capture vortex shedding phenomena at appropriate Reynolds numbers
- **Geometry-Specific Optimizations**: 
  - Rectangle geometries use finer meshes and robust DistMesh parameters for convergence
  - Ellipse geometries balance computational efficiency with accuracy
  - Cylinder geometries use the original optimized parameters
- **Robust Mesh Generation**: Custom `distmesh2d_robust.m` handles challenging geometries with sharp corners

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