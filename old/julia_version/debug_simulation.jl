#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)

include("RBFNS.jl")
using .RBFNS

# Run debug version with detailed diagnostics
function debug_simulation()
    println("🔬 Debug Simulation - Tracking all matrix dimensions and values")
    
    # Create config
    cfg = default_config()
    
    # Use simplified config with conservative parameters
    println("   Using conservative parameters: dt=1e-3, stencils=[15,10,8,12]")
    
    # Build geometry and operators exactly like simulate.jl
    println("1. Building geometry...")
    G = build_geometry(cfg)
    
    # Get arrays for convenience
    xy = G.xy
    xy1 = G.xy1
    xy_s = G.xy_s  
    xy1_s = G.xy1_s
    
    println("   V-grid: $(size(xy, 1)) interior + $(size(xy1, 1) - size(xy, 1)) boundary = $(size(xy1, 1)) total")
    println("   P-grid: $(size(xy_s, 1)) interior + $(size(xy1_s, 1) - size(xy_s, 1)) boundary = $(size(xy1_s, 1)) total")
    
    # Build nearest neighbors
    println("2. Building nearest neighbors...")
    k = cfg.stencil_size_main
    Nearest_Idx = nearest_interp(xy1, xy1, k)
    Nearest_Idx_s = nearest_interp(xy1_s, xy1_s, k)
    
    # Build operators with minimal RBF-FD
    println("3. Building operators...")
    D_all = rbf_phs_fd_all(xy1, xy1, Nearest_Idx, cfg.order_main, cfg.poly_degree_main, cfg.laplacian_order)
    Dx, Dy, L = D_all[1], D_all[2], D_all[3]
    
    D_s_all = rbf_phs_fd_all(xy_s, xy1_s, Nearest_Idx_s[1:size(xy_s, 1), :], cfg.order_main, cfg.poly_degree_main, cfg.laplacian_order)
    L_s = D_s_all[3]
    
    # Build simplified inter-grid operators (P-grid to V-grid gradients)
    k_interp = min(15, size(xy1_s, 1))
    Nearest_Idx_interp = nearest_interp(xy1, xy1_s, k_interp)
    D_interp = rbf_phs_fd_all(xy1, xy1_s, Nearest_Idx_interp, cfg.order_boundary, cfg.poly_degree_boundary, cfg.derivative_order)
    D0_21_x, D0_21_y = D_interp[1], D_interp[2]
    
    # Build simplified V-grid to P-grid operators  
    Nearest_Idx_21 = nearest_interp(xy_s, xy1, min(15, size(xy1, 1)))
    D_21 = rbf_phs_fd_all(xy_s, xy1, Nearest_Idx_21, cfg.order_boundary, cfg.poly_degree_boundary, cfg.derivative_order)
    D0_12_x, D0_12_y = D_21[1], D_21[2]
    
    println("   Operator sizes:")
    println("   L: $(size(L))")
    println("   L_s: $(size(L_s))")
    println("   D0_21_x: $(size(D0_21_x))")
    println("   D0_12_x: $(size(D0_12_x))")
    
    # Build SIMPLIFIED pressure system (avoid complex boundary operators)
    println("4. Building SIMPLIFIED pressure system...")
    n_interior_s = size(G.xy_s, 1)
    boundary_s_nodes = vcat(G.boundary_y_s, G.boundary_out_s, G.boundary_in_s, G.boundary_obs_s)
    n_boundary_s = size(boundary_s_nodes, 1)
    
    # Use simplified Dirichlet boundary conditions (p = 0 at boundaries)
    L_pressure = copy(L_s[1:n_interior_s, :])
    L_bc_simple = spzeros(n_boundary_s, size(xy1_s, 1))
    
    # Simple identity matrix for boundary conditions
    for i in 1:n_boundary_s
        L_bc_simple[i, n_interior_s + i] = 1.0
    end
    
    L_pressure = vcat(L_pressure, L_bc_simple)
    
    # Add regularization EXACTLY like MATLAB
    reg_col = vcat(ones(n_interior_s), zeros(n_boundary_s))
    reg_row = hcat(ones(1, n_interior_s), zeros(1, n_boundary_s))
    
    L_pressure = hcat(L_pressure, reg_col)
    L_pressure = vcat(L_pressure, hcat(reg_row, [0.0]))
    
    println("   Pressure system size: $(size(L_pressure))")
    
    # Build solver
    L_inv_s = build_pressure_system(sparse(L_pressure))
    
    # Initialize state
    println("5. Initializing simulation state...")
    n_total = size(xy1, 1)
    U = zeros(n_total)
    V = zeros(n_total)
    
    # Set inlet boundary condition: U = 1 at inlet
    inlet_start = size(xy, 1) + size(G.boundary_y, 1) + size(G.boundary_out, 1) + 1
    inlet_end = inlet_start + size(G.boundary_in, 1) - 1
    U[inlet_start:inlet_end] .= 1.0
    
    W = vcat(U, V)
    W_prev = copy(W)
    
    # Define boundary parameters
    L_B = size(G.boundary_in, 1) + size(G.boundary_obs, 1)
    L_W = size(xy, 1)
    L_B_y = size(G.boundary_y, 1)
    L_B_S = n_boundary_s
    
    println("   Boundary parameters: L_B=$L_B, L_W=$L_W, L_B_y=$L_B_y, L_B_S=$L_B_S")
    
    # Simplified boundary operators
    Dy_b = spzeros(L_B_y, n_total)
    Dy_b_1 = ones(L_B_y)
    
    # Time stepping parameters (conservative)
    dt = 1e-3  # Reduced time step
    nu = cfg.viscosity
    ab_curr = cfg.adams_bashforth_current
    ab_prev = cfg.adams_bashforth_previous
    cn = cfg.crank_nicolson
    
    println("   Time parameters: dt=$dt, nu=$nu")
    
    # Run a few time steps with detailed diagnostics
    println("\n6. Running time steps with diagnostics...")
    for step in 1:10
        println("\n--- Step $step ---")
        
        # Extract current velocity components
        U_curr = W[1:n_total]
        V_curr = W[n_total+1:end]
        U_prev_step = W_prev[1:n_total]
        V_prev_step = W_prev[n_total+1:end]
        
        println("   Initial max|U|: $(maximum(abs.(U_curr))), max|V|: $(maximum(abs.(V_curr)))")
        
        # Check for NaN/Inf in initial state
        if any(!isfinite.(U_curr)) || any(!isfinite.(V_curr))
            println("   ❌ NaN/Inf detected in initial velocity at step $step")
            break
        end
        
        # STEP 1: Adams-Bashforth advection step
        adv_U = ab_curr * (Dx * (U_curr .* U_curr) + Dy * (U_curr .* V_curr)) - 
                ab_prev * (Dx * (U_prev_step .* U_prev_step) + Dy * (U_prev_step .* V_prev_step))
        adv_V = ab_curr * (Dx * (U_curr .* V_curr) + Dy * (V_curr .* V_curr)) - 
                ab_prev * (Dx * (U_prev_step .* V_prev_step) + Dy * (V_prev_step .* V_prev_step))
        
        println("   Advection max|adv_U|: $(maximum(abs.(adv_U))), max|adv_V|: $(maximum(abs.(adv_V)))")
        
        if any(!isfinite.(adv_U)) || any(!isfinite.(adv_V))
            println("   ❌ NaN/Inf detected in advection at step $step")
            break
        end
        
        # STEP 2: Crank-Nicolson diffusion step  
        U_rhs = U_curr - dt * adv_U + dt * nu * cn * L * U_curr
        V_rhs = V_curr - dt * adv_V + dt * nu * cn * L * V_curr
        
        println("   Diffusion max|U_rhs|: $(maximum(abs.(U_rhs))), max|V_rhs|: $(maximum(abs.(V_rhs)))")
        
        if any(!isfinite.(U_rhs)) || any(!isfinite.(V_rhs))
            println("   ❌ NaN/Inf detected in diffusion at step $step")
            break
        end
        
        # Apply simplified boundary conditions
        U_rhs[inlet_start:inlet_end] .= 1.0  # Inlet: U = 1
        V_rhs[inlet_start:inlet_end] .= 0.0  # Inlet: V = 0
        
        # Simplified wall and obstacle BCs
        U_rhs[end-L_B+1:end] .= 0.0  # Obstacle: U = 0
        V_rhs[end-L_B+1:end] .= 0.0  # Obstacle: V = 0
        
        U2 = U_rhs
        V2 = V_rhs
        
        println("   After BCs max|U2|: $(maximum(abs.(U2))), max|V2|: $(maximum(abs.(V2)))")
        
        # STEP 3: Pressure correction
        div_full = D0_12_x * U2 + D0_12_y * V2
        F0 = div_full[1:end-L_B_S]  # Interior divergence only
        
        println("   Divergence: size=$(length(F0)), max|F0|: $(maximum(abs.(F0)))")
        
        if any(!isfinite.(F0))
            println("   ❌ NaN/Inf detected in divergence at step $step")
            break
        end
        
        # Solve pressure
        F = vcat(F0, zeros(L_B_S + 1))
        p_full = L_inv_s(F)
        p = p_full[1:end-1]  # Remove regularization
        
        println("   Pressure: size=$(length(p)), max|p|: $(maximum(abs.(p)))")
        
        if any(!isfinite.(p))
            println("   ❌ NaN/Inf detected in pressure at step $step")
            break
        end
        
        # STEP 4: Velocity correction
        grad_p_x = D0_21_x * p
        grad_p_y = D0_21_y * p
        
        println("   Pressure gradients: max|grad_p_x|: $(maximum(abs.(grad_p_x))), max|grad_p_y|: $(maximum(abs.(grad_p_y)))")
        
        if any(!isfinite.(grad_p_x)) || any(!isfinite.(grad_p_y))
            println("   ❌ NaN/Inf detected in pressure gradients at step $step")
            break
        end
        
        # Apply velocity correction (simplified - no padding for now)
        U3 = U2 - grad_p_x[1:length(U2)]
        V3 = V2 - grad_p_y[1:length(V2)]
        
        println("   Final velocities: max|U3|: $(maximum(abs.(U3))), max|V3|: $(maximum(abs.(V3)))")
        
        if any(!isfinite.(U3)) || any(!isfinite.(V3))
            println("   ❌ NaN/Inf detected in final velocities at step $step")
            break
        end
        
        # Update for next time step
        W_prev = copy(W)
        W = vcat(U3, V3)
        
        println("   ✅ Step $step completed successfully")
    end
    
    println("\n🏁 Debug simulation completed")
end

debug_simulation()
