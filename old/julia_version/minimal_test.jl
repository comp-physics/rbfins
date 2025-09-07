#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

include("RBFNS.jl")
using .RBFNS

# Minimal test to isolate the bounds error
function minimal_test()
    println("🔬 Minimal Test - Isolating RBF-FD bounds error")
    
    # Simple test points
    xy1 = [0.0 0.0; 1.0 0.0; 0.0 1.0; 1.0 1.0; 0.5 0.5]
    xy_s = copy(xy1)
    
    k = 4  # Use all points as neighbors
    
    println("Test points: $(size(xy1))")
    println("Source points: $(size(xy_s))")
    println("Stencil size: $k")
    
    # Build nearest neighbors manually
    Nearest_Idx = nearest_interp(xy1, xy_s, k)
    println("Nearest_Idx size: $(size(Nearest_Idx))")
    println("Nearest_Idx:")
    display(Nearest_Idx)
    
    # Test RBF-FD call
    try
        D_all = rbf_phs_fd_all(xy1, xy_s, Nearest_Idx, 3, 3, 2)
        println("✅ RBF-FD successful")
        println("   D_x size: $(size(D_all[1]))")
        println("   D_y size: $(size(D_all[2]))")
    catch e
        println("❌ RBF-FD failed: $e")
        
        # Check bounds
        println("\nDebugging bounds:")
        println("xy1 size: $(size(xy1))")
        println("xy_s size: $(size(xy_s))")
        println("Nearest_Idx bounds: min=$(minimum(Nearest_Idx)), max=$(maximum(Nearest_Idx))")
        println("Should be in range [1, $(size(xy_s, 1))]")
    end
    
    # Test with full geometry
    println("\n--- Testing with full geometry ---")
    cfg = default_config()
    G = build_geometry(cfg)
    
    xy = G.xy[1:100, :]  # Use only first 100 points
    xy_s = G.xy_s[1:100, :]
    
    k = 10
    
    println("Full geometry test:")
    println("xy size: $(size(xy))")
    println("xy_s size: $(size(xy_s))")
    println("Stencil size: $k")
    
    try
        Nearest_Idx = nearest_interp(xy, xy_s, k)
        println("Nearest_Idx size: $(size(Nearest_Idx))")
        println("Nearest_Idx bounds: min=$(minimum(Nearest_Idx)), max=$(maximum(Nearest_Idx))")
        
        if maximum(Nearest_Idx) > size(xy_s, 1)
            println("❌ Index out of bounds! Max index $(maximum(Nearest_Idx)) > $(size(xy_s, 1))")
            return
        end
        
        D_all = rbf_phs_fd_all(xy, xy_s, Nearest_Idx, 8, 3, 1)  # Smaller order
        println("✅ Full geometry RBF-FD successful")
        
    catch e
        println("❌ Full geometry RBF-FD failed: $e")
    end
end

minimal_test()
