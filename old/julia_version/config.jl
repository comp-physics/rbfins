"""
Configuration structure for RBF-FD Navier-Stokes simulation
"""
struct Config
    # Domain Configuration
    x_min::Float64
    x_max::Float64
    y_min::Float64
    y_max::Float64
    
    # Mesh Generation Parameters
    dist::Float64
    boundary_eps::Float64
    refine_a1::Float64
    refine_b1::Float64
    refine_a2::Float64
    refine_b2::Float64
    edge_multiplier::Int
    
    # Geometry Parameters
    geometry_type::String  # "cylinder" or "ellipse"
    obstacle_radius::Float64
    ellipse_a::Float64
    ellipse_b::Float64
    
    # RBF-FD Algorithm Parameters
    stencil_size_main::Int
    stencil_size_boundary_obstacle::Int
    stencil_size_boundary_wall::Int
    stencil_size_boundary_outlet::Int
    
    order_main::Int
    poly_degree_main::Int
    laplacian_order::Int
    
    order_boundary::Int
    poly_degree_boundary::Int
    derivative_order::Int
    
    order_interpolation_low::Int
    poly_degree_interpolation_low::Int
    
    order_interpolation_high::Int
    order_interpolation_high_poly::Int
    poly_degree_interpolation_high::Int
    
    order_near_obstacle::Int
    order_near_obstacle_poly::Int
    poly_degree_near_obstacle::Int
    
    order_near_boundary::Int
    poly_degree_near_boundary::Int
    
    # Simulation Parameters
    reynolds_number::Float64
    viscosity::Float64
    time_step::Float64
    num_time_steps::Int
    num_time_steps_ci::Int
    random_seed::Int
    show_progress::Bool
    
    # Distance Thresholds
    x_min_dist::Float64
    x_max_dist::Float64
    y_min_dist::Float64
    y_max_dist::Float64
    
    # Numerical Scheme Coefficients
    adams_bashforth_current::Float64
    adams_bashforth_previous::Float64
    crank_nicolson::Float64
    
    # Visualization Parameters
    scatter_size::Int
    plot_tick_y::Vector{Float64}
    plot_tick_x::Vector{Float64}
    color_axis_range::Float64
end

"""
    default_config()

Create default configuration matching the MATLAB config.m
"""
function default_config()
    return Config(
        # Domain Configuration
        -8.0, 24.0, -8.0, 8.0,
        
        # Mesh Generation Parameters
        0.1, 0.002, 0.05, 0.08, 0.05, 0.08, 3,
        
        # Geometry Parameters
        "cylinder", 0.5, 0.6, 0.4,
        
        # RBF-FD Algorithm Parameters (conservative for stability)
        15, 10, 8, 12,
        28, 3, 3,
        8, 3, 1,
        12, 2,
        25, 5, 3,
        30, 7, 4,
        17, 2,
        
        # Simulation Parameters  
        100.0, 0.01, 2e-3, 5000, 20, 42, true,
        
        # Distance Thresholds
        1.0, 1.0, 0.5, 0.5,
        
        # Numerical Scheme Coefficients
        3/2, 1/2, 1/2,
        
        # Visualization Parameters
        15, [-5.0, 0.0, 5.0], [-5.0, 0.0, 5.0, 10.0, 15.0], 1e-0
    )
end
