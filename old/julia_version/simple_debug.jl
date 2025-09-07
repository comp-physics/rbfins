#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

include("RBFNS.jl")
using .RBFNS

# Simplified debug version to identify exact location of NaN
function simple_debug()
    println("🔬 Simple Debug - Finding NaN source")
    
    # Create geometry
    cfg = default_config()
    G = build_geometry(cfg)
    
    # Get node arrays
    xy = G.xy
    xy_s = G.xy_s
    
    # Create combined grids (like simulate.jl)
    xy1 = vcat(xy, G.boundary_y, G.boundary_out, G.boundary_in, G.boundary_obs)
    xy1_s = vcat(xy_s, G.boundary_y_s, G.boundary_out_s, G.boundary_in_s, G.boundary_obs_s)
    
    println("V-grid: $(size(xy1, 1)) total, P-grid: $(size(xy1_s, 1)) total")
    
    # Simple operators 
    k = 15
    println("Building operators with k=$k...")
    
    # Try to identify if RBF-FD is the source of NaN
    try
        Nearest_Idx = nearest_interp(xy1, xy1, k)
        println("✅ Nearest neighbors built successfully")
        
        D_all = rbf_phs_fd_all(xy1[1:10, :], xy1, Nearest_Idx[1:10, :], 28, 3, 3)
        println("✅ Small RBF-FD test successful")
        
        D_all = rbf_phs_fd_all(xy1, xy1, Nearest_Idx, 28, 3, 3)
        Dx, Dy, L = D_all[1], D_all[2], D_all[3]
        
        # Check for NaN in operators
        if any(!isfinite.(Dx)) || any(!isfinite.(Dy)) || any(!isfinite.(L))
            println("❌ NaN detected in main operators")
            return
        end
        println("✅ Main operators built successfully")
        
        # Test pressure operators
        Nearest_Idx_s = nearest_interp(xy1_s, xy1_s, k)
        D_s_all = rbf_phs_fd_all(xy_s, xy1_s, Nearest_Idx_s[1:size(xy_s, 1), :], 28, 3, 3)
        L_s = D_s_all[3]
        
        if any(!isfinite.(L_s))
            println("❌ NaN detected in pressure operators")
            return
        end
        println("✅ Pressure operators built successfully")
        
        # Test a simple time step
        println("Testing simple time step...")
        
        n_total = size(xy1, 1)
        U = zeros(n_total) 
        V = zeros(n_total)
        
        # Set simple inlet BC
        inlet_start = size(xy, 1) + size(G.boundary_y, 1) + size(G.boundary_out, 1) + 1
        inlet_end = inlet_start + size(G.boundary_in, 1) - 1
        U[inlet_start:inlet_end] .= 1.0
        
        println("Initial state: max|U|=$(maximum(abs.(U))), max|V|=$(maximum(abs.(V)))")
        
        # Test advection 
        adv_U = Dx * (U .* U) + Dy * (U .* V)
        adv_V = Dx * (U .* V) + Dy * (V .* V)
        
        if any(!isfinite.(adv_U)) || any(!isfinite.(adv_V))
            println("❌ NaN in advection terms")
            return
        end
        println("✅ Advection terms OK")
        
        # Test diffusion
        diff_U = L * U
        diff_V = L * V
        
        if any(!isfinite.(diff_U)) || any(!isfinite.(diff_V))
            println("❌ NaN in diffusion terms")
            return  
        end
        println("✅ Diffusion terms OK")
        
        # Test RHS construction
        dt = 1e-3
        nu = 0.01
        U_rhs = U - dt * adv_U + dt * nu * 0.5 * diff_U
        V_rhs = V - dt * adv_V + dt * nu * 0.5 * diff_V
        
        if any(!isfinite.(U_rhs)) || any(!isfinite.(V_rhs))
            println("❌ NaN in RHS construction") 
            return
        end
        println("✅ RHS construction OK")
        
        # Test velocity solve (just copy for now)
        U2 = copy(U_rhs)
        V2 = copy(V_rhs)
        
        # Apply BCs
        U2[inlet_start:inlet_end] .= 1.0
        V2[inlet_start:inlet_end] .= 0.0
        
        if any(!isfinite.(U2)) || any(!isfinite.(V2))
            println("❌ NaN after BC application")
            return
        end
        println("✅ BC application OK")
        
        # Test divergence
        k_interp = min(10, size(xy1_s, 1))
        Nearest_Idx_12 = nearest_interp(xy_s, xy1, k_interp)
        D_12 = rbf_phs_fd_all(xy_s, xy1, Nearest_Idx_12, 8, 3, 1)
        D0_12_x, D0_12_y = D_12[1], D_12[2]
        
        if any(!isfinite.(D0_12_x)) || any(!isfinite.(D0_12_y))
            println("❌ NaN in inter-grid operators")
            return
        end
        
        div_test = D0_12_x * U2 + D0_12_y * V2
        
        if any(!isfinite.(div_test))
            println("❌ NaN in divergence calculation")
            println("   D0_12_x range: [$(minimum(D0_12_x)), $(maximum(D0_12_x))]")
            println("   D0_12_y range: [$(minimum(D0_12_y)), $(maximum(D0_12_y))]")
            println("   U2 range: [$(minimum(U2)), $(maximum(U2))]")
            println("   V2 range: [$(minimum(V2)), $(maximum(V2))]")
            return
        end
        println("✅ Divergence calculation OK")
        
        println("\n🎉 All basic operations successful - NaN likely occurs in more complex interactions")
        
    catch e
        println("❌ Error in debug: $e")
    end
end

simple_debug()
