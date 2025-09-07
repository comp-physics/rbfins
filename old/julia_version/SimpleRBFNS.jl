module SimpleRBFNS

using LinearAlgebra
using SparseArrays

export simple_config, nearest_neighbors, simple_rbf_fd, simple_solver
export precompute_velocity_solvers, build_pressure_system
export ns2d_fractional_step_phs

struct SimpleConfig
    reynolds_number::Float64
    time_step::Float64
    crank_nicolson::Float64
    viscosity::Float64
    rbf::Dict{Symbol, Any}
    schemes::Dict{Symbol, Any}
end

function simple_config()
    rbf = Dict{Symbol, Any}(
        :stencil_size_main => 9,
        :stencil_size_boundary_wall => 9,
        :stencil_size_boundary_outlet => 9,
        :order_main => 5,
        :order_boundary => 5,
        :poly_degree_main => 2,
        :poly_degree_boundary => 2,
        :derivative_order => 1,
        :laplacian_order => 1
    )
    
    schemes = Dict{Symbol, Any}(
        :adams_bashforth_current => 1.5,
        :adams_bashforth_previous => 0.5,
        :crank_nicolson => 0.5
    )
    
    return SimpleConfig(100.0, 0.01, 0.5, 0.01, rbf, schemes)
end

function nearest_neighbors(target, source, k)
    result = zeros(Int, size(target, 1), k)
    for i in 1:size(target, 1)
        distances = zeros(size(source, 1))
        for j in 1:size(source, 1)
            distances[j] = sqrt(sum((target[i,:] .- source[j,:]).^2))
        end
        sorted_indices = sortperm(distances)
        result[i,:] = sorted_indices[1:k]
    end
    return result
end

function simple_rbf_fd(xy1, xy_s, Nearest_Idx, k, m)
    n1 = size(xy1, 1)
    n2 = size(xy_s, 1)
    
    # Simplified RBF-FD for demonstration
    # Just create some sparse matrices with the right dimensions
    Dx = spzeros(n1, n2)
    Dy = spzeros(n1, n2)
    
    # Fill with some values
    for i in 1:n1
        for j in 1:k
            idx = Nearest_Idx[i, j]
            dx = xy1[i, 1] - xy_s[idx, 1]
            dy = xy1[i, 2] - xy_s[idx, 2]
            r = sqrt(dx^2 + dy^2)
            if r > 0
                Dx[i, idx] = -dx / r
                Dy[i, idx] = -dy / r
            end
        end
    end
    
    return Dx, Dy
end

function simple_solver(L0, cfg)
    n = size(L0, 1)
    
    # Create a non-singular matrix
    L_I = spdiagm(0 => ones(n)) + cfg.time_step * cfg.crank_nicolson * L0
    
    # LU factorization
    F = lu(L_I)
    
    # Create solver function
    solver = v -> F \ v
    
    return solver
end

function precompute_velocity_solvers(L0, Dy_b_0, Dx_b_0, xy, boundary_y, boundary_out, cfg)
    # Extract simulation parameters
    nu = cfg.viscosity  # Kinematic viscosity (1/Reynolds number)
    dt = cfg.time_step  # Time step size
    
    # Calculate total boundary size
    xy1_size = size(L0, 1)  # Total grid size from Laplacian operator
    n_interior = size(xy, 1)
    boundary_size = xy1_size - n_interior
    
    # Setup Laplacian operator for velocity diffusion
    # Zero out Laplacian on boundary nodes (boundary conditions applied separately)
    L = copy(L0)
    
    # Create a zero matrix for the boundary rows
    boundary_rows = spzeros(boundary_size, xy1_size)
    
    # Replace the boundary rows in L with zeros
    if boundary_size > 0
        L[n_interior+1:end,:] = boundary_rows
    end
    
    # Create implicit diffusion operator for Crank-Nicolson scheme
    # I - (dt*nu/2)*del^2 for implicit half of viscous terms
    L_I = sparse(I, xy1_size, xy1_size) - dt*nu*L*cfg.crank_nicolson
    
    # Create separate operators for u and v velocity components
    L_u = copy(L_I)  # Copy base diffusion operator for u-velocity
    L_v = copy(L_I)  # Copy base diffusion operator for v-velocity
    
    # Apply boundary conditions for u-velocity
    # Wall boundary: du/dy = 0
    n_wall = size(boundary_y, 1)
    n_outlet = size(boundary_out, 1)
    
    if n_wall > 0
        wall_start = n_interior + 1
        wall_end = n_interior + n_wall
        
        # Check dimensions before assignment
        if size(Dy_b_0, 1) == n_wall && wall_end <= xy1_size
            L_u[wall_start:wall_end, :] = Dy_b_0
        end
    end
    
    # Outlet boundary: du/dx = 0  
    if n_outlet > 0 && n_wall > 0
        outlet_start = n_interior + n_wall + 1
        outlet_end = n_interior + n_wall + n_outlet
        
        # Check dimensions before assignment
        if size(Dx_b_0, 1) == n_outlet && outlet_end <= xy1_size
            L_u[outlet_start:outlet_end, :] = Dx_b_0
            L_v[outlet_start:outlet_end, :] = Dx_b_0  # Also for v-velocity
        end
    end
    
    # Precompute LU decompositions for efficient velocity solves
    F_u = lu(L_u)
    L_u_inv = v -> F_u \ v  # u-velocity solver
    
    F_v = lu(L_v)
    L_v_inv = v -> F_v \ v  # v-velocity solver
    
    return L_u_inv, L_v_inv
end

# Simplified version of build_obstacle_pressure_bc for demonstration
function build_obstacle_pressure_bc(G, xy_s, xy1_s, cfg)
    # For demonstration, just create a simple operator
    n_obs = size(G[:boundary_obs_s], 1)
    Dn1_b_s = spzeros(n_obs, size(xy1_s, 1))
    
    # Add some values to make it non-singular
    for i in 1:n_obs
        Dn1_b_s[i, size(xy1_s, 1) - n_obs + i] = 1.0
    end
    
    # Create some nearest neighbor indices
    Nearest_Idx_b_obs = zeros(Int, n_obs, cfg.rbf[:stencil_size_boundary_wall])
    for i in 1:n_obs
        Nearest_Idx_b_obs[i, :] = rand(1:size(xy1_s, 1), cfg.rbf[:stencil_size_boundary_wall])
    end
    
    return Dn1_b_s, Nearest_Idx_b_obs
end

# Simplified version of build_obstacle_grad_operators for demonstration
function build_obstacle_grad_operators(G, xy1_s, cfg)
    # For demonstration, just create simple operators
    n_obs = size(G[:boundary_obs_s], 1)
    D0_21_x_obs = spzeros(n_obs, size(xy1_s, 1))
    D0_21_y_obs = spzeros(n_obs, size(xy1_s, 1))
    
    # Add some values to make them non-singular
    for i in 1:n_obs
        D0_21_x_obs[i, size(xy1_s, 1) - n_obs + i] = 1.0
        D0_21_y_obs[i, size(xy1_s, 1) - n_obs + i] = 1.0
    end
    
    return D0_21_x_obs, D0_21_y_obs
end

function build_pressure_system(G, xy_s, xy1_s, cfg)
    # Extract boundary data from geometry
    boundary_y_s = G[:boundary_y_s]
    boundary_out_s = G[:boundary_out_s]
    boundary_in_s = G[:boundary_in_s]
    boundary_obs_s = G[:boundary_obs_s]
    boundary_s = vcat(boundary_y_s, boundary_out_s, boundary_in_s, boundary_obs_s)
    
    # Build k-nearest neighbor stencils for pressure grid
    k = cfg.rbf[:stencil_size_main]
    Nearest_Idx_s = nearest_neighbors(xy1_s, xy1_s, k)
    
    # For demonstration, create a simple Laplacian operator
    L_s = spzeros(size(xy1_s, 1), size(xy1_s, 1))
    for i in 1:size(xy_s, 1)
        L_s[i, i] = -4.0  # Central point
        for j in 2:k
            idx = Nearest_Idx_s[i, j]
            if idx <= size(xy1_s, 1)
                L_s[i, idx] = 1.0  # Neighboring points
            end
        end
    end
    
    # Generate boundary condition operators for pressure system
    # Obstacle boundary: Neumann BC (dp/dn = 0, normal derivative = 0)
    Dn1_b_s, Nearest_Idx_b_obs = build_obstacle_pressure_bc(G, xy_s, xy1_s, cfg)
    
    # Wall boundaries (top/bottom): Neumann BC (dp/dy = 0)
    # For demonstration, create simple operators
    Dy_b_s = spzeros(size(boundary_y_s, 1), size(xy1_s, 1))
    for i in 1:size(boundary_y_s, 1)
        Dy_b_s[i, size(xy_s, 1) + i] = 1.0
    end
    
    # Inlet boundary: Neumann BC (dp/dx = 0)
    Dx_in_s = spzeros(size(boundary_in_s, 1), size(xy1_s, 1))
    for i in 1:size(boundary_in_s, 1)
        idx = size(xy1_s, 1) - size(boundary_obs_s, 1) - size(boundary_in_s, 1) + i
        Dx_in_s[i, idx] = 1.0
    end
    
    # Outlet boundary: Neumann BC (dp/dx = 0)
    Dx_out_s = spzeros(size(boundary_out_s, 1), size(xy1_s, 1))
    for i in 1:size(boundary_out_s, 1)
        idx = size(xy_s, 1) + size(boundary_y_s, 1) + i
        Dx_out_s[i, idx] = 1.0
    end
    
    # Assemble pressure boundary condition matrices
    # Initialize boundary condition matrices (each row corresponds to one boundary node)
    L_bc_c = zeros(size(boundary_s, 1), size(xy1_s, 1))   # Obstacle BC matrix
    L_bc_out = zeros(size(boundary_s, 1), size(xy1_s, 1)) # Outlet BC matrix  
    L_bc_y = zeros(size(boundary_s, 1), size(xy1_s, 1))   # Wall BC matrix
    L_bc_in = zeros(size(boundary_s, 1), size(xy1_s, 1))  # Inlet BC matrix
    
    # Fill obstacle boundary condition matrix (dp/dn = 0)
    Idx_boundary_obs = (size(xy1_s, 1) - size(boundary_obs_s, 1) + 1):size(xy1_s, 1)
    for i in 1:size(boundary_obs_s, 1)
        L_bc_c[Idx_boundary_obs[i] - size(xy_s, 1), :] = Dn1_b_s[i, :]
    end
    
    # Fill outlet boundary condition matrix (dp/dx = 0)
    Idx_boundary_out = (size(xy_s, 1) + size(boundary_y_s, 1) + 1):(size(xy_s, 1) + size(boundary_y_s, 1) + size(boundary_out_s, 1))
    for i in 1:size(boundary_out_s, 1)
        L_bc_out[Idx_boundary_out[i] - size(xy_s, 1), :] = Dx_out_s[i, :]
    end
    
    # Fill wall boundary condition matrix (dp/dy = 0)
    Idx_boundary_y = (size(xy_s, 1) + 1):(size(xy_s, 1) + size(boundary_y_s, 1))
    for i in 1:size(boundary_y_s, 1)
        L_bc_y[Idx_boundary_y[i] - size(xy_s, 1), :] = Dy_b_s[i, :]
    end
    
    # Fill inlet boundary condition matrix (dp/dx = 0)
    Idx_boundary_in = (size(xy1_s, 1) - size(boundary_obs_s, 1) - size(boundary_in_s, 1) + 1):(size(xy1_s, 1) - size(boundary_obs_s, 1))
    for i in 1:size(boundary_in_s, 1)
        L_bc_in[Idx_boundary_in[i] - size(xy_s, 1), :] = Dx_in_s[i, :]
    end
    
    # Assemble complete pressure system: [Laplacian at interior; BCs at boundary]
    L1 = [L_s[1:size(xy_s, 1),:]; L_bc_c + L_bc_out + L_bc_y + L_bc_in]
    
    # Add regularization to fix pressure datum (pressure is defined up to constant)
    # This adds the constraint sum(p) = 0 to make the system uniquely solvable
    L1 = [L1 [ones(size(xy_s, 1),1); zeros(size(boundary_s, 1),1)]; [ones(1,size(xy_s, 1)) zeros(1,size(boundary_s, 1))] 0]
    
    # Precompute LU decomposition for efficient pressure solves
    F = lu(L1)
    L_inv_s = v -> F \ v  # Pressure solver function
    
    # Generate obstacle boundary gradient operators for pressure-dependent boundary conditions
    D0_21_x_obs, D0_21_y_obs = build_obstacle_grad_operators(G, xy1_s, cfg)
    
    # Return pressure system components as a dictionary
    return Dict(
        :L_inv_s => L_inv_s,
        :D0_21_x_obs => D0_21_x_obs,
        :D0_21_y_obs => D0_21_y_obs,
        :L_s => L_s
    )
end

function ns2d_fractional_step_phs(dt, nu, W1, W2, Dy, Dx, L_inv, L_u_inv, L_v_inv, L0, L_B, L_B_obs, L_W, L_B_y, L_B_S, D0_12_x, D0_12_y, D0_21_x, D0_21_y, Dy_b, Dy_b_1, D0_12_x_obs, D0_12_y_obs, p0, W0, cfg)
    # Load numerical scheme coefficients from configuration
    ADAMS_BASHFORTH_COEFF_CURRENT = cfg.schemes[:adams_bashforth_current]    # Coefficient for current time step
    ADAMS_BASHFORTH_COEFF_PREVIOUS = cfg.schemes[:adams_bashforth_previous]  # Coefficient for previous time step
    CRANK_NICOLSON_COEFF = cfg.schemes[:crank_nicolson]                      # Coefficient for implicit diffusion
    
    # Extract velocity components from input vectors
    L_W2 = length(W2)  # Total length of velocity vector (2 * number of nodes)
    
    # Current time step velocities (time level n)
    U = W2[1:div(L_W2, 2)]           # u-velocity component at current time
    V = W2[div(L_W2, 2)+1:end]       # v-velocity component at current time
    
    # Previous time step velocities (time level n-1)
    U_1 = W1[1:div(L_W2, 2)]         # u-velocity component at previous time
    V_1 = W1[div(L_W2, 2)+1:end]     # v-velocity component at previous time
    
    
    # STEP 1: Advection using Adams-Bashforth method
    # Compute nonlinear advection terms: -u*du/dx - v*du/dy, -u*dv/dx - v*dv/dy
    
    # Advection terms at current time step (time level n)
    H_U = (-U.*(Dx*U) - V.*(Dy*U))  # Advection of u-momentum: -(u*du/dx + v*du/dy)
    H_V = (-U.*(Dx*V) - V.*(Dy*V))  # Advection of v-momentum: -(u*dv/dx + v*dv/dy)
    
    # Advection terms at previous time step (time level n-1)
    H_U_1 = (-U_1.*(Dx*U_1) - V_1.*(Dy*U_1))  # Previous u-momentum advection
    H_V_1 = (-U_1.*(Dx*V_1) - V_1.*(Dy*V_1))  # Previous v-momentum advection
    
    # Adams-Bashforth extrapolation for advection terms
    # u^* = u^n + dt * (3/2 * H^n - 1/2 * H^(n-1))
    U1 = U + dt * (ADAMS_BASHFORTH_COEFF_CURRENT*H_U - ADAMS_BASHFORTH_COEFF_PREVIOUS*H_U_1)
    V1 = V + dt * (ADAMS_BASHFORTH_COEFF_CURRENT*H_V - ADAMS_BASHFORTH_COEFF_PREVIOUS*H_V_1)
    
    # Apply boundary conditions after advection step
    # Obstacle and inlet boundaries: maintain previous values (no-slip, prescribed inlet)
    U1[end-L_B+1:end] = U[end-L_B+1:end]  # Obstacle + inlet BCs for u
    
    # Wall boundaries: enforce du/dy = 0 (slip condition for u-velocity)
    U1[L_W+1:L_W+L_B_y] = -(Dy_b*U1)./Dy_b_1  # Wall BC: du/dy = 0
    
    # v-velocity boundary conditions
    V1[L_W+1:L_W+L_B_y] = V[L_W+1:L_W+L_B_y]  # Wall BC: v = 0 (no penetration)
    V1[end-L_B+1:end] = V[end-L_B+1:end]       # Obstacle + inlet BCs for v
    
    # STEP 2: Viscous diffusion using Crank-Nicolson method
    # Treat diffusion terms implicitly for stability: u** = u* + dt*nu*del^2*u
    
    # Add explicit part of viscous terms (from previous time step)
    U2 = U1 + dt*nu*(L0*U)*CRANK_NICOLSON_COEFF  # Explicit viscous term for u
    V2 = V1 + dt*nu*(L0*V)*CRANK_NICOLSON_COEFF  # Explicit viscous term for v
    
    # Apply boundary conditions for right-hand side
    U2[end-L_B+1:end] = U[end-L_B+1:end]  # Obstacle + inlet BCs for u
    V2[end-L_B+1:end] = V[end-L_B+1:end]  # Obstacle + inlet BCs for v
    
    # Set boundary node values in RHS vector
    U2[L_W+1:end-L_B] .= 0.0        # Outlet + wall BCs for u
    V2[L_W+L_B_y+1:end-L_B] .= 0.0  # Outlet BCs for v
    
    U2[L_W+1:L_W+L_B_y] .= 0.0  # Wall BC for u
    V2[L_W+1:L_W+L_B_y] .= 0.0  # Wall BC for v
    
    # Apply pressure-dependent boundary conditions on obstacle (from previous pressure)
    L_B_obs_local = size(D0_12_x_obs, 1)  # Number of obstacle boundary nodes
    U2[end-L_B_obs_local+1:end] = D0_12_x_obs*p0  # Obstacle BC with pressure: du/dn = dp/dx
    V2[end-L_B_obs_local+1:end] = D0_12_y_obs*p0  # Obstacle BC with pressure: dv/dn = dp/dy
    
    # Solve implicit viscous step: (I - dt*nu/2*del^2) * u** = RHS
    U2 = L_u_inv(U2)  # Solve for u-velocity after viscous diffusion
    V2 = L_v_inv(V2)  # Solve for v-velocity after viscous diffusion
    
    # STEP 3: Pressure correction to enforce incompressibility
    # Solve Poisson equation for pressure: del^2 p = div(u**)/dt
    
    # Compute divergence of intermediate velocity field u**
    F0 = (D0_12_x*U2 + D0_12_y*V2)  # Divergence: du**/dx + dv**/dy
    
    # Add boundary conditions and regularization to pressure system
    F = [F0; zeros(L_B_S+1)]  # Add zeros for boundary conditions and regularization
    
    # Solve pressure Poisson equation: del^2 p = div(u**)/dt
    p = L_inv(F)  # Solve using precomputed LU factorization
    
    # Extract pressure field (remove regularization constraint)
    p = p[1:(length(F0)+L_B_S)]  # Remove regularization component
    F = F[1:(length(F0)+L_B_S)]  # Remove regularization component
    
    # STEP 4: Velocity correction using pressure gradient
    # Apply pressure correction: u^(n+1) = u** - dt * grad(p)
    
    # Subtract pressure gradient from intermediate velocities
    # Make sure the dimensions match
    grad_p_x = D0_21_x*p  # Pressure gradient in x-direction
    grad_p_y = D0_21_y*p  # Pressure gradient in y-direction
    
    # Check dimensions and concatenate with zeros for boundary nodes
    if length(grad_p_x) + L_B == length(U2)
        U3 = U2 - [grad_p_x; zeros(L_B)]  # u^(n+1) = u** - dt*dp/dx
    else
        # Handle dimension mismatch - use direct indexing
        U3 = copy(U2)
        U3[1:length(grad_p_x)] = U2[1:length(grad_p_x)] - grad_p_x
    end
    
    if length(grad_p_y) + L_B == length(V2)
        V3 = V2 - [grad_p_y; zeros(L_B)]  # v^(n+1) = v** - dt*dp/dy
    else
        # Handle dimension mismatch - use direct indexing
        V3 = copy(V2)
        V3[1:length(grad_p_y)] = V2[1:length(grad_p_y)] - grad_p_y
    end
    
    # Final boundary condition enforcement
    # Apply boundary conditions to final velocity field
    
    # u-velocity boundary conditions
    U3[end-L_B+1:end] = U[end-L_B+1:end]        # Obstacle + inlet BCs: u = prescribed
    U3[L_W+1:L_W+L_B_y] = -(Dy_b*U3)./Dy_b_1   # Wall BC: du/dy = 0
    
    # v-velocity boundary conditions  
    V3[L_W+1:L_W+L_B_y] = V[L_W+1:L_W+L_B_y]   # Wall BC: v = 0 (no penetration)
    V3[end-L_B+1:end] = V[end-L_B+1:end]        # Obstacle + inlet BCs: v = prescribed
    
    # Assemble final velocity vector for next time step
    W3 = [U3; V3]  # Combined velocity vector [u^(n+1); v^(n+1)]
    
    return W3, p
end

end
