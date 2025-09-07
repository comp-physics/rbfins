include("SimpleRBFNS.jl")
using .SimpleRBFNS
using LinearAlgebra
using SparseArrays

# Test configuration
cfg = simple_config()
println("Reynolds number: $(cfg.reynolds_number)")
println("Time step: $(cfg.time_step)")

# Test nearest neighbors
target = [0.0 0.0; 1.0 1.0; 2.0 2.0]
source = [0.1 0.1; 1.1 1.1; 2.1 2.1; 3.0 3.0; -1.0 -1.0]
k = 3
neighbors = nearest_neighbors(target, source, k)
println("Neighbor matrix size: $(size(neighbors))")
println("First target point neighbors: $(neighbors[1, :])")

# Test RBF-FD
xy1 = rand(5, 2)
xy_s = rand(10, 2)
Nearest_Idx = nearest_neighbors(xy1, xy_s, 4)
Dx, Dy = simple_rbf_fd(xy1, xy_s, Nearest_Idx, 4, 3)
println("Dx size: $(size(Dx))")
println("Dy size: $(size(Dy))")

# Test solver
n = 10
L0 = spdiagm(0 => 2.0*ones(n), 1 => -ones(n-1), -1 => -ones(n-1))
solver = simple_solver(L0, cfg)
b = rand(n)
x = solver(b)
println("Solution norm: $(norm(x))")
println("Residual: $(norm(L0*x - b))")
