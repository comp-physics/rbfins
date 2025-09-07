using LinearAlgebra
using SparseArrays

"""
    ns2d_fractional_step_phs(dt, nu, W1, W2, Dy, Dx, L_inv, L_u_inv, L_v_inv, L0,
                              L_B, L_B_obs, L_W, L_B_y, L_B_S,
                              D0_12_x, D0_12_y, D0_21_x, D0_21_y,
                              Dy_b, Dy_b_1, D0_12_x_obs, D0_12_y_obs, p0, W0, cfg)

Fractional step method for 2D incompressible Navier-Stokes equations.

This function implements the fractional step method for incompressible Navier-Stokes:
1. Advection-diffusion step using Adams-Bashforth + Crank-Nicolson
2. Pressure correction step to enforce incompressibility  
3. Velocity correction using pressure gradient

INPUTS:
  dt      - Time step size
  nu      - Kinematic viscosity (1/Reynolds number)
  W1, W2  - Velocity fields at previous time steps [U1; V1], [U2; V2]
  Dy, Dx  - Spatial derivative operators (y and x directions)
  L_inv   - Pressure Poisson solver (function handle)
  L_u_inv, L_v_inv - Velocity diffusion solvers (function handles)
  L0      - Laplacian operator for velocity diffusion
  L_B, L_B_obs, L_W, L_B_y, L_B_S - Boundary indexing parameters
  D0_12_x, D0_12_y - Divergence operators (V-grid to P-grid)
  D0_21_x, D0_21_y - Gradient operators (P-grid to V-grid)
  Dy_b, Dy_b_1     - Boundary derivative operators
  D0_12_x_obs, D0_12_y_obs - Obstacle boundary operators
  p0      - Previous pressure field
  W0      - Initial velocity field (for reference)
  cfg     - Configuration structure

OUTPUTS:
  W3      - Updated velocity field [U3; V3] at next time step
  p       - Updated pressure field
"""
function ns2d_fractional_step_phs(dt, nu, W1, W2, Dy, Dx, L_inv, L_u_inv, L_v_inv, L0,
                                  L_B, L_B_obs, L_W, L_B_y, L_B_S,
                                  D0_12_x, D0_12_y, D0_21_x, D0_21_y,
                                  Dy_b, Dy_b_1, D0_12_x_obs, D0_12_y_obs, p0, W0, cfg)
    
    # Extract numerical scheme coefficients from configuration
    ADAMS_BASHFORTH_COEFF_CURRENT = cfg.adams_bashforth_current    # 3/2 coefficient for current time step
    ADAMS_BASHFORTH_COEFF_PREVIOUS = cfg.adams_bashforth_previous  # 1/2 coefficient for previous time step
    CRANK_NICOLSON_COEFF = cfg.crank_nicolson                     # 1/2 coefficient for implicit diffusion

    # Extract velocity components from input vectors
    L_W2 = length(W2)  # Total length of velocity vector (2 * number of nodes)

    # Current time step velocities (time level n)
    U = @view W2[1:L_W2÷2]           # u-velocity component at current time
    V = @view W2[L_W2÷2+1:end]       # v-velocity component at current time

    # Previous time step velocities (time level n-1)
    U_1 = @view W1[1:L_W2÷2]         # u-velocity component at previous time
    V_1 = @view W1[L_W2÷2+1:end]     # v-velocity component at previous time

    ## STEP 1: Advection using Adams-Bashforth method
    # Compute nonlinear advection terms: -u*du/dx - v*du/dy, -u*dv/dx - v*dv/dy

    # Advection terms at current time step (time level n)
    H_U = -(U .* (Dx * U) .+ V .* (Dy * U))  # Advection of u-momentum: -(u*du/dx + v*du/dy)
    H_V = -(U .* (Dx * V) .+ V .* (Dy * V))  # Advection of v-momentum: -(u*dv/dx + v*dv/dy)

    # Advection terms at previous time step (time level n-1)
    H_U_1 = -(U_1 .* (Dx * U_1) .+ V_1 .* (Dy * U_1))  # Previous u-momentum advection
    H_V_1 = -(U_1 .* (Dx * V_1) .+ V_1 .* (Dy * V_1))  # Previous v-momentum advection

    # Adams-Bashforth extrapolation for advection terms
    # u^* = u^n + dt * (3/2 * H^n - 1/2 * H^(n-1))
    U1 = U .+ dt .* (ADAMS_BASHFORTH_COEFF_CURRENT .* H_U .- ADAMS_BASHFORTH_COEFF_PREVIOUS .* H_U_1)
    V1 = V .+ dt .* (ADAMS_BASHFORTH_COEFF_CURRENT .* H_V .- ADAMS_BASHFORTH_COEFF_PREVIOUS .* H_V_1)

    # Apply boundary conditions after advection step
    # Obstacle and inlet boundaries: maintain previous values (no-slip, prescribed inlet)
    U1[end-L_B+1:end] .= U[end-L_B+1:end]  # Obstacle + inlet BCs for u

    # Wall boundaries: enforce du/dy = 0 (slip condition for u-velocity)
    U1[L_W+1:L_W+L_B_y] .= -(Dy_b * U1) ./ Dy_b_1  # Wall BC: du/dy = 0

    # v-velocity boundary conditions
    V1[L_W+1:L_W+L_B_y] .= V[L_W+1:L_W+L_B_y]  # Wall BC: v = 0 (no penetration)
    V1[end-L_B+1:end] .= V[end-L_B+1:end]       # Obstacle + inlet BCs for v

    ## STEP 2: Viscous diffusion using Crank-Nicolson method
    # Treat diffusion terms implicitly for stability: u** = u* + dt*nu*∇²*u

    # Add explicit part of viscous terms (from previous time step)
    U2 = U1 .+ dt * nu .* (L0 * U) .* CRANK_NICOLSON_COEFF  # Explicit viscous term for u
    V2 = V1 .+ dt * nu .* (L0 * V) .* CRANK_NICOLSON_COEFF  # Explicit viscous term for v

    # Apply boundary conditions for right-hand side
    U2[end-L_B+1:end] .= U[end-L_B+1:end]  # Obstacle + inlet BCs for u
    V2[end-L_B+1:end] .= V[end-L_B+1:end]  # Obstacle + inlet BCs for v

    # Set boundary node values in RHS vector
    U2[L_W+1:end-L_B] = zeros(L_W2÷2-L_W-L_B)        # Outlet + wall BCs for u
    V2[L_W+L_B_y+1:end-L_B] = zeros(L_W2÷2-L_W-L_B-L_B_y)  # Outlet BCs for v

    U2[L_W+1:L_W+L_B_y] = zeros(L_B_y)  # Wall BC for u
    V2[L_W+1:L_W+L_B_y] = zeros(L_B_y)  # Wall BC for v

    # Apply pressure-dependent boundary conditions on obstacle (from previous pressure)
    L_B_obs_local = size(D0_12_x_obs, 1)  # Number of obstacle boundary nodes
    U2[end-L_B_obs_local+1:end] .= D0_12_x_obs * p0  # Obstacle BC with pressure: du/dn = dp/dx
    V2[end-L_B_obs_local+1:end] .= D0_12_y_obs * p0  # Obstacle BC with pressure: dv/dn = dp/dy

    # Solve implicit viscous step: (I - dt*nu/2*∇²) * u** = RHS
    U2 = L_u_inv(U2)  # Solve for u-velocity after viscous diffusion
    V2 = L_v_inv(V2)  # Solve for v-velocity after viscous diffusion

    ## STEP 3: Pressure correction to enforce incompressibility
    # Solve Poisson equation for pressure: ∇²p = div(u**)/dt

    # Compute divergence of intermediate velocity field u**
    div_full = (D0_12_x * U2 .+ D0_12_y * V2)  # Divergence: du**/dx + dv**/dy
    # Extract only interior divergence (matching MATLAB where F0 has size length(xy_s))
    F0 = div_full[1:length(div_full)-L_B_S]  # Use only interior divergence

    # Add boundary conditions and regularization to pressure system
    F = vcat(F0, zeros(eltype(F0), L_B_S + 1))  # Add zeros for boundary conditions and regularization

    # Solve pressure Poisson equation: ∇²p = div(u**)/dt
    p = L_inv(F)  # Solve using precomputed LU factorization

    # Extract pressure field (remove regularization constraint)
    p = p[1:length(F0)+L_B_S]  # Remove regularization component

    ## STEP 4: Velocity correction using pressure gradient
    # Apply pressure correction: u^(n+1) = u** - dt * grad(p)

    # Subtract pressure gradient from intermediate velocities  
    # MATLAB approach: [(D0_21_x*p); zeros(L_B,1)]
    grad_p_x = D0_21_x * p
    grad_p_y = D0_21_y * p
    
    # Check if we need to pad with zeros to match velocity vector size
    if length(grad_p_x) + L_B == length(U2)
        # MATLAB-style concatenation
        U3 = U2 .- vcat(grad_p_x, zeros(L_B))  # u^(n+1) = u** - dp/dx
        V3 = V2 .- vcat(grad_p_y, zeros(L_B))  # v^(n+1) = v** - dp/dy
    else
        # Fallback: truncate/pad as needed
        U3 = U2 .- grad_p_x[1:length(U2)]  # u^(n+1) = u** - dp/dx
        V3 = V2 .- grad_p_y[1:length(V2)]  # v^(n+1) = v** - dp/dy
    end

    ## Final boundary condition enforcement
    # Apply boundary conditions to final velocity field

    # u-velocity boundary conditions
    U3[end-L_B+1:end] .= U[end-L_B+1:end]        # Obstacle + inlet BCs: u = prescribed
    U3[L_W+1:L_W+L_B_y] .= -(Dy_b * U3) ./ Dy_b_1   # Wall BC: du/dy = 0

    # v-velocity boundary conditions
    V3[L_W+1:L_W+L_B_y] .= V[L_W+1:L_W+L_B_y]   # Wall BC: v = 0 (no penetration)
    V3[end-L_B+1:end] .= V[end-L_B+1:end]        # Obstacle + inlet BCs: v = prescribed

    # Assemble final velocity vector for next time step
    W3 = vcat(U3, V3)  # Combined velocity vector [u^(n+1); v^(n+1)]

    return W3, p
end
