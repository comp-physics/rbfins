using LinearAlgebra
using SparseArrays

"""
    build_pressure_system(L1)

Build pressure system solver from assembled pressure matrix.

This function creates a pressure Poisson system solver from the pre-assembled
system matrix including Laplacian, boundary conditions, and regularization.

INPUTS:
  L1 - Complete pressure system matrix [interior Laplacian; boundary conditions; regularization]

OUTPUTS:
  L_inv_s - Function handle for solving pressure system
"""
function build_pressure_system(L1::SparseMatrixCSC)
    # Precompute LU decomposition for efficient pressure solves
    # Add regularization to handle near-singular pressure matrices
    n = size(L1, 1)
    regularization = 1e-12 * I(n)
    L1_reg = L1 + regularization
    
    F_p = try
        lu(L1_reg)  # Sparse LU factorization (UMFPACK)
    catch e
        if isa(e, SingularException)
            # If still singular, add more regularization
            L1_reg2 = L1 + 1e-10 * I(n)
            lu(L1_reg2)
        else
            rethrow(e)
        end
    end
    
    # Create solver function handle
    L_inv_s = (b -> F_p \ b)  # Pressure solver
    
    return L_inv_s
end

"""
    build_complete_pressure_system(G, xy_s, xy1_s, cfg)

Build complete pressure system including boundary conditions and regularization.

This is a placeholder for the full pressure system construction that would
include building boundary condition operators, assembling the system matrix,
and adding regularization constraints.

INPUTS:
  G      - Geometry structure
  xy_s   - Interior pressure nodes  
  xy1_s  - Complete pressure grid (interior + boundary)
  cfg    - Configuration structure

OUTPUTS:
  P - Structure containing:
      .L_inv_s        - Pressure solver function handle
      .D0_21_x_obs    - Obstacle x-gradient operator
      .D0_21_y_obs    - Obstacle y-gradient operator
      .L_s            - Laplacian operator
"""
function build_complete_pressure_system(G, xy_s, xy1_s, cfg)
    # This would implement the full MATLAB build_pressure_system.m functionality
    # For now, this is a structural placeholder
    error("build_complete_pressure_system not fully implemented - use build_pressure_system with pre-assembled matrix")
end
