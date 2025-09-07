open("example.jl", "w") do f
    write(f, """
#!/usr/bin/env julia

"""
Example usage of the Julia RBF-FD Navier-Stokes solver components.
"""

using Pkg
Pkg.activate(@__DIR__)

try
    include("RBFNS.jl")
    using .RBFNS
    using LinearAlgebra, SparseArrays

    println("Julia RBF-FD Navier-Stokes Solver - Example Usage")
    println("=" ^ 60)

    # Create configuration
    println("\\n1. Creating configuration...")
    cfg = default_config()
    println("   Reynolds number: \$(cfg.reynolds_number)")
    println("   Time step: \$(cfg.time_step)")
    println("   Geometry: \$(cfg.geometry_type)")

    # Example 1: Nearest neighbor search
    println("\\n2. Testing nearest neighbor search...")
    n_target = 10
    n_source = 20
    target_points = rand(n_target, 2) .* 10 .- 5  # Random points in [-5, 5]²
    source_points = rand(n_source, 2) .* 10 .- 5
    k = 5

    println("   Finding \$k nearest neighbors for \$n_target target points in \$n_source source points...")
    neighbors = nearest_interp(target_points, source_points, k)
    println("   Neighbor matrix size: \$(size(neighbors))")
    println("   First target point neighbors: \$(neighbors[1, :])")
catch e
    println("Error: \$e")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end
""")
end

println("Created example.jl")
