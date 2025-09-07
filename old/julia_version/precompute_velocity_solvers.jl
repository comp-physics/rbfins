using LinearAlgebra
using SparseArrays

"""
    precompute_velocity_solvers(L0, Dy_b_0, Dx_b_0, xy_len, boundary_y_len, boundary_out_len, cfg)

Create velocity system solvers with boundary conditions for the Crank-Nicolson scheme.

This function creates the implicit diffusion operators for the Crank-Nicolson
scheme and precomputes their LU decompositions for efficient solving during
time stepping.

INPUTS:
  L0              - Laplacian operator on velocity grid
  Dy_b_0          - Wall boundary condition operator (du/dy = 0)
  Dx_b_0          - Outlet boundary condition operator (du/dx = 0)
  xy_len          - Number of interior velocity nodes
  boundary_y_len  - Number of wall boundary nodes
  boundary_out_len - Number of outlet boundary nodes
  cfg             - Configuration structure

OUTPUTS:
  L_u_inv - Function handle for solving u-velocity system
  L_v_inv - Function handle for solving v-velocity system
"""
function precompute_velocity_solvers(L0::SparseMatrixCSC, Dy_b_0::SparseMatrixCSC, Dx_b_0::SparseMatrixCSC,
                                     xy_len::Int, boundary_y_len::Int, boundary_out_len::Int, cfg)
    
    # Extract simulation parameters
    nu = cfg.viscosity  # Kinematic viscosity (1/Reynolds number)
    dt = cfg.time_step  # Time step size
    crank_nicolson = cfg.crank_nicolson  # Crank-Nicolson coefficient
    
    N = size(L0, 1)  # Total grid size from Laplacian operator
    boundary_size = N - xy_len

    # Setup Laplacian operator for velocity diffusion
    # Zero out Laplacian on boundary nodes (boundary conditions applied separately)
    L = copy(L0)
    if boundary_size > 0
        L[xy_len+1:end, :] .= 0.0
    end

    # Create implicit diffusion operator for Crank-Nicolson scheme
    # I - (dt*nu/2)*∇² for implicit half of viscous terms
    L_I = spdiagm(0 => ones(N)) - (dt * nu * crank_nicolson) * L

    # Create separate operators for u and v velocity components
    L_u = copy(L_I)  # Copy base diffusion operator for u-velocity
    L_v = copy(L_I)  # Copy base diffusion operator for v-velocity

    # Apply boundary conditions for u-velocity
    # Wall boundary: du/dy = 0
    wall_start = xy_len + 1
    wall_end = xy_len + boundary_y_len
    if boundary_y_len > 0
        L_u[wall_start:wall_end, :] = Dy_b_0
    end

    # Outlet boundary: du/dx = 0
    outlet_start = xy_len + boundary_y_len + 1
    outlet_end = xy_len + boundary_y_len + boundary_out_len
    if boundary_out_len > 0
        L_u[outlet_start:outlet_end, :] = Dx_b_0
    end

    # Apply boundary conditions for v-velocity
    # Outlet boundary: dv/dx = 0
    if boundary_out_len > 0
        L_v[outlet_start:outlet_end, :] = Dx_b_0
    end

    # Note: Wall and inlet boundaries have Dirichlet conditions (v=0) applied directly
    # Obstacle boundary has no-slip condition (u=v=0) applied directly

    # Precompute LU decompositions for efficient velocity solves
    F_u = lu(L_u)  # Sparse LU factorization (UMFPACK)
    F_v = lu(L_v)  # Sparse LU factorization (UMFPACK)

    # Create solver function handles
    L_u_inv = (b -> F_u \ b)  # u-velocity solver
    L_v_inv = (b -> F_v \ b)  # v-velocity solver

    return L_u_inv, L_v_inv
end
