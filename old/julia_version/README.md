# Julia Port of RBF-FD Navier-Stokes Solver

This directory contains a Julia port of the MATLAB RBF-FD incompressible Navier-Stokes solver.

## Overview

The original MATLAB code implements a mesh-free Radial Basis Function Finite Difference (RBF-FD) method for solving 2D incompressible Navier-Stokes equations around obstacles (cylinder, ellipse) using:

- **RBF-FD discretization**: Polyharmonic Spline (PHS) basis functions with polynomial augmentation
- **Fractional step method**: Adams-Bashforth advection + Crank-Nicolson diffusion + pressure correction
- **Staggered grids**: Separate velocity (V-grid) and pressure (P-grid) node distributions

## Files

### Core Modules
- `RBFNS.jl` - Main module that includes all components
- `config.jl` - Configuration structure and default parameters
- `nearest_interp.jl` - K-nearest neighbor search using KD-trees
- `RBF_PHS_FD_all.jl` - RBF-FD differentiation matrix construction
- `precompute_velocity_solvers.jl` - Velocity system factorization for Crank-Nicolson
- `build_pressure_system.jl` - Pressure Poisson system construction
- `NS_2d_fractional_step_PHS.jl` - Main time-stepping algorithm

### Simulation and Utilities
- `simulate.jl` - Complete simulation driver (equivalent to simulate.m)
- `geometry.jl` - Mesh generation and boundary extraction
- `plotting.jl` - Visualization utilities using Plots.jl
- `example.jl` - Basic usage examples and component tests
- `test_simple.jl` - Comprehensive test suite

### Dependencies
- `Project.toml` - Julia package dependencies

## Installation

1. Navigate to the `julia_version/` directory
2. Start Julia and activate the environment:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

## Dependencies

- **LinearAlgebra.jl** - Linear algebra operations and sparse matrix support
- **SparseArrays.jl** - Sparse matrix operations (included in Julia standard library)
- **NearestNeighbors.jl** - Efficient k-nearest neighbor search using KD-trees

## Usage

### Quick Start - Run Complete Simulation
```julia
using Pkg
Pkg.activate("julia_version")

# Run the complete simulation
include("simulate.jl")
```

### Component Testing
```julia
# Run basic component tests
include("test_simple.jl")

# Run individual examples
include("example.jl")
```

### Manual Usage
```julia
include("RBFNS.jl")
using .RBFNS

# Create default configuration
cfg = default_config()

# Generate geometry
G = build_geometry(cfg)

# Build RBF-FD operators
xy1 = vcat(G.xy, G.boundary_in, G.boundary_out, G.boundary_y, G.boundary_obs)
neighbors = nearest_interp(xy1, xy1, cfg.stencil_size_main)
Dx, Dy, L, Dxx, Dyy, Dxy = rbf_phs_fd_all(xy1, xy1, neighbors, 
                                           cfg.stencil_size_main, cfg.order_main, cfg.poly_degree_main)

# Visualize results (requires Plots.jl)
using Plots
plot_mesh(G)
```

## Key Differences from MATLAB

1. **1-based vs 0-based indexing**: Julia uses 1-based indexing like MATLAB
2. **Function syntax**: Julia uses `function` keyword and different array syntax
3. **Broadcasting**: Julia uses `.` for element-wise operations (e.g., `x .+ y`)
4. **Sparse matrices**: Similar interface but different underlying implementation
5. **Matrix operations**: Julia's `\` operator for linear solves, `lu()` for factorization
6. **Memory management**: Julia's garbage collector vs MATLAB's automatic memory management

## Performance Notes

- Julia's sparse matrix operations use UMFPACK (same as MATLAB)
- JIT compilation may cause slower first runs but faster subsequent executions
- Memory allocation patterns may differ from MATLAB

## Porting Status

✅ **Complete**: Core algorithms and data structures  
✅ **Complete**: RBF-FD operator construction  
✅ **Complete**: Fractional step time integration  
✅ **Complete**: Configuration management  
✅ **Complete**: Geometry generation (simplified mesh approach)  
✅ **Complete**: Visualization utilities (Plots.jl integration)  
✅ **Complete**: Full simulation driver (`simulate.jl`)  
✅ **Complete**: Test suite and examples  

⚠️ **Simplified**: Boundary condition setup (basic implementation)  
⚠️ **Missing**: Advanced mesh generation (DistMesh equivalent for irregular meshes)  
⚠️ **Missing**: Streamline computation and advanced visualization  

## Next Steps

1. Implement geometry generation (mesh generation around obstacles)
2. Add visualization capabilities using Plots.jl or Makie.jl
3. Create complete simulation driver equivalent to `simulate.m`
4. Add test suite to verify numerical equivalence with MATLAB version
5. Performance optimization and benchmarking

## References

- [1] T. Chu, O. T. Schmidt, "RBF-FD discretization of the Navier-Stokes equations on scattered but staggered nodes", Journal of Computational Physics 474, 111756, 2023
- [2] T. Chu, O. T. Schmidt, "Mesh-free hydrodynamic stability", Submitted to Journal of Computational Physics
