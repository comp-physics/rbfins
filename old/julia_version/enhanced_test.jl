include("SimpleRBFNS.jl")
using .SimpleRBFNS
using LinearAlgebra
using SparseArrays

println("Julia RBF-FD Navier-Stokes Solver - Testing")
println("=" ^ 60)

# Test configuration
cfg = simple_config()
println("Reynolds number: $(cfg.reynolds_number)")
println("Time step: $(cfg.time_step)")
println("Crank-Nicolson parameter: $(cfg.crank_nicolson)")
println("Viscosity: $(cfg.viscosity)")

# Test nearest neighbors
println("\nTesting nearest neighbors search...")
target = [0.0 0.0; 1.0 1.0; 2.0 2.0]
source = [0.1 0.1; 1.1 1.1; 2.1 2.1; 3.0 3.0; -1.0 -1.0]
k = 3
neighbors = nearest_neighbors(target, source, k)
println("Neighbor matrix size: $(size(neighbors))")
println("First target point neighbors: $(neighbors[1, :])")

# Test RBF-FD
println("\nTesting RBF-FD operators...")
xy1 = rand(5, 2)
xy_s = rand(10, 2)
Nearest_Idx = nearest_neighbors(xy1, xy_s, 4)
Dx, Dy = simple_rbf_fd(xy1, xy_s, Nearest_Idx, 4, 3)
println("Dx size: $(size(Dx))")
println("Dy size: $(size(Dy))")

# Test simple solver
println("\nTesting simple solver...")
n = 10
L0 = spdiagm(0 => 2.0*ones(n), 1 => -ones(n-1), -1 => -ones(n-1))
solver = simple_solver(L0, cfg)
b = rand(n)
x = solver(b)
println("Solution norm: $(norm(x))")
println("Residual: $(norm(L0*x - b))")

# Test velocity solvers
println("\nTesting velocity solvers...")
n_interior = 8
n_wall = 2
n_outlet = 2
n_total = n_interior + n_wall + n_outlet

# Create mock grid
xy = rand(n_interior, 2)
boundary_y = rand(n_wall, 2)
boundary_out = rand(n_outlet, 2)

# Create mock operators - ensure non-singularity
L0 = spdiagm(0 => 2.0*ones(n_total), 1 => -ones(n_total-1), -1 => -ones(n_total-1))

# Create boundary operators with proper dimensions
Dy_b_0 = spzeros(n_wall, n_total)
Dx_b_0 = spzeros(n_outlet, n_total)

# Fill boundary operators with values that will ensure non-singularity
for i in 1:n_wall
    Dy_b_0[i, n_interior+i] = 1.0  # Diagonal term
    if i < n_wall
        Dy_b_0[i, n_interior+i+1] = -0.5  # Off-diagonal term
    end
end

for i in 1:n_outlet
    Dx_b_0[i, n_interior+n_wall+i] = 1.0  # Diagonal term
    if i < n_outlet
        Dx_b_0[i, n_interior+n_wall+i+1] = -0.5  # Off-diagonal term
    end
end

# Precompute velocity solvers
L_u_inv, L_v_inv = precompute_velocity_solvers(L0, Dy_b_0, Dx_b_0, xy, boundary_y, boundary_out, cfg)

# Test solvers
b_u = rand(n_total)
b_v = rand(n_total)
u = L_u_inv(b_u)
v = L_v_inv(b_v)

println("u-velocity solution norm: $(norm(u))")
println("v-velocity solution norm: $(norm(v))")

# Test pressure system
println("\nTesting pressure system...")

# Create mock geometry
n_interior_p = 6
n_wall_p = 2
n_outlet_p = 2
n_inlet_p = 2
n_obs_p = 3
n_total_p = n_interior_p + n_wall_p + n_outlet_p + n_inlet_p + n_obs_p

# Create mock grid for pressure
xy_s = rand(n_interior_p, 2)
xy1_s = rand(n_total_p, 2)

# Create mock geometry structure with dictionary
G = Dict{Symbol, Any}(
    :boundary_y_s => rand(n_wall_p, 2),
    :boundary_out_s => rand(n_outlet_p, 2),
    :boundary_in_s => rand(n_inlet_p, 2),
    :boundary_obs_s => rand(n_obs_p, 2)
)

# Build pressure system
P = build_pressure_system(G, xy_s, xy1_s, cfg)

# Test pressure solver
b_p = rand(n_total_p + 1)  # +1 for regularization constraint
p = P[:L_inv_s](b_p)

println("Pressure system size: $(n_total_p + 1)")
println("Pressure solution norm: $(norm(p))")
println("D0_21_x_obs size: $(size(P[:D0_21_x_obs]))")
println("D0_21_y_obs size: $(size(P[:D0_21_y_obs]))")

# Test fractional step method
println("\nTesting fractional step method...")

# Create mock velocity field
n_v = 20  # Number of velocity nodes
n_p = 15  # Number of pressure nodes
n_b = 5   # Number of boundary nodes

# Create mock velocity fields with consistent dimensions
W1 = rand(2*n_v)  # Previous velocity
W2 = rand(2*n_v)  # Current velocity
W0 = rand(2*n_v)  # Initial velocity

# Create mock pressure field
p0 = rand(n_p)

# Mock operators for the fractional step method
Dx_v = sprand(n_v, n_v, 0.2)  # Velocity x-derivative
Dy_v = sprand(n_v, n_v, 0.2)  # Velocity y-derivative
L0_v = spdiagm(0 => -4.0*ones(n_v), 1 => ones(n_v-1), -1 => ones(n_v-1))  # Laplacian

# Boundary parameters
L_B = n_b  # Number of boundary nodes
L_B_obs = 3  # Number of obstacle boundary nodes
L_W = n_v - n_b  # Number of interior nodes
L_B_y = 2  # Number of wall boundary nodes
L_B_S = n_p - (n_v - n_b)  # Boundary size for pressure

# Mock divergence and gradient operators with consistent dimensions
D0_12_x = sprand(n_p, n_v, 0.2)  # Divergence x-component
D0_12_y = sprand(n_p, n_v, 0.2)  # Divergence y-component
D0_21_x = sprand(n_v - n_b, n_p, 0.2)  # Gradient x-component (interior nodes only)
D0_21_y = sprand(n_v - n_b, n_p, 0.2)  # Gradient y-component (interior nodes only)

# Mock boundary operators
Dy_b = sprand(L_B_y, n_v, 0.2)  # Wall boundary derivative
Dy_b_1 = -ones(L_B_y)  # Wall boundary normalization
D0_12_x_obs = sprand(L_B_obs, n_p, 0.2)  # Obstacle x-derivative
D0_12_y_obs = sprand(L_B_obs, n_p, 0.2)  # Obstacle y-derivative

# Mock solvers
L_inv = x -> x ./ 2.0  # Mock pressure solver
L_u_inv_mock = x -> x ./ 1.5  # Mock u-velocity solver
L_v_inv_mock = x -> x ./ 1.5  # Mock v-velocity solver

# Run fractional step method
W3, p = ns2d_fractional_step_phs(
    cfg.time_step, cfg.viscosity, W1, W2, 
    Dy_v, Dx_v, L_inv, L_u_inv_mock, L_v_inv_mock, 
    L0_v, L_B, L_B_obs, L_W, L_B_y, L_B_S, 
    D0_12_x, D0_12_y, D0_21_x, D0_21_y, 
    Dy_b, Dy_b_1, D0_12_x_obs, D0_12_y_obs, 
    p0, W0, cfg
)

println("Velocity field size: $(size(W3))")
println("Pressure field size: $(size(p))")
println("Velocity solution norm: $(norm(W3))")
println("Pressure solution norm: $(norm(p))")
