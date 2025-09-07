using Plots
using Printf

"""
Visualization utilities for RBF-FD Navier-Stokes results.

This module provides plotting functions for mesh visualization,
velocity fields, pressure fields, and streamlines.
"""

"""
    plot_mesh(G::GeometryData; title="Computational Mesh")

Plot the computational mesh with boundaries highlighted.
"""
function plot_mesh(G::GeometryData; title="Computational Mesh")
    p = scatter(G.xy[:, 1], G.xy[:, 2], 
               marker=:circle, markersize=2, markercolor=:blue, 
               label="Interior (V-grid)", alpha=0.6)
    
    scatter!(p, G.xy_s[:, 1], G.xy_s[:, 2],
            marker=:square, markersize=2, markercolor=:red,
            label="Interior (P-grid)", alpha=0.6)
    
    # Boundary nodes
    if size(G.boundary_in, 1) > 0
        scatter!(p, G.boundary_in[:, 1], G.boundary_in[:, 2],
                marker=:circle, markersize=3, markercolor=:green,
                label="Inlet")
    end
    
    if size(G.boundary_out, 1) > 0
        scatter!(p, G.boundary_out[:, 1], G.boundary_out[:, 2],
                marker=:circle, markersize=3, markercolor=:orange,
                label="Outlet")
    end
    
    if size(G.boundary_y, 1) > 0
        scatter!(p, G.boundary_y[:, 1], G.boundary_y[:, 2],
                marker=:circle, markersize=3, markercolor=:purple,
                label="Walls")
    end
    
    if size(G.boundary_obs, 1) > 0
        scatter!(p, G.boundary_obs[:, 1], G.boundary_obs[:, 2],
                marker=:circle, markersize=4, markercolor=:black,
                label="Obstacle")
    end
    
    plot!(p, title=title, xlabel="x", ylabel="y", aspect_ratio=:equal,
          legend=:outertopright)
    
    return p
end

"""
    plot_velocity_field(xy, u, v; title="Velocity Field", subsample=10)

Plot velocity field as arrows on the computational domain.
"""
function plot_velocity_field(xy::AbstractMatrix, u::AbstractVector, v::AbstractVector; 
                            title="Velocity Field", subsample=10)
    
    # Subsample for cleaner visualization
    n = size(xy, 1)
    step = max(1, n ÷ subsample)
    indices = 1:step:n
    
    x_sub = xy[indices, 1]
    y_sub = xy[indices, 2] 
    u_sub = u[indices]
    v_sub = v[indices]
    
    # Create quiver plot
    p = quiver(x_sub, y_sub, quiver=(u_sub, v_sub),
              arrow=arrow(:closed, :head, 0.3),
              color=:blue, alpha=0.7)
    
    plot!(p, title=title, xlabel="x", ylabel="y", aspect_ratio=:equal)
    
    return p
end

"""
    plot_velocity_magnitude(xy, u, v; title="Velocity Magnitude")

Plot velocity magnitude as a color contour.
"""
function plot_velocity_magnitude(xy::AbstractMatrix, u::AbstractVector, v::AbstractVector;
                                title="Velocity Magnitude")
    
    vel_mag = sqrt.(u.^2 .+ v.^2)
    
    p = scatter(xy[:, 1], xy[:, 2], 
               marker_z=vel_mag, markersize=4,
               color=:viridis, alpha=0.8)
    
    plot!(p, title=title, xlabel="x", ylabel="y", 
          aspect_ratio=:equal, colorbar_title="‖v‖")
    
    return p
end

"""
    plot_pressure_field(xy_s, p; title="Pressure Field")

Plot pressure field as a color contour.
"""
function plot_pressure_field(xy_s::AbstractMatrix, p::AbstractVector;
                            title="Pressure Field")
    
    p_plot = scatter(xy_s[:, 1], xy_s[:, 2],
                    marker_z=p, markersize=4,
                    color=:RdBu, alpha=0.8)
    
    plot!(p_plot, title=title, xlabel="x", ylabel="y",
          aspect_ratio=:equal, colorbar_title="p")
    
    return p_plot
end

"""
    plot_vorticity(xy, Dxy, u, v; title="Vorticity Field")

Plot vorticity field (ω = ∂v/∂x - ∂u/∂y).
"""
function plot_vorticity(xy::AbstractMatrix, Dx::AbstractMatrix, Dy::AbstractMatrix, 
                       u::AbstractVector, v::AbstractVector; title="Vorticity Field")
    
    # Compute vorticity: ω = ∂v/∂x - ∂u/∂y
    dvdx = Dx * v
    dudy = Dy * u
    omega = dvdx - dudy
    
    p = scatter(xy[:, 1], xy[:, 2],
               marker_z=omega, markersize=4,
               color=:RdBu, alpha=0.8)
    
    plot!(p, title=title, xlabel="x", ylabel="y",
          aspect_ratio=:equal, colorbar_title="ω")
    
    return p
end

"""
    plot_simulation_summary(xy, xy_s, u, v, p, G; time_step=0)

Create a comprehensive summary plot with multiple subplots.
"""
function plot_simulation_summary(xy::AbstractMatrix, xy_s::AbstractMatrix,
                                u::AbstractVector, v::AbstractVector, p::AbstractVector,
                                G::GeometryData; time_step=0)
    
    # Create subplot layout
    p1 = plot_mesh(G, title="Computational Mesh")
    p2 = plot_velocity_magnitude(xy, u, v, title="Velocity Magnitude")
    p3 = plot_velocity_field(xy, u, v, title="Velocity Vectors", subsample=20)
    p4 = plot_pressure_field(xy_s, p, title="Pressure Field")
    
    main_title = "RBF-FD Navier-Stokes Simulation"
    if time_step > 0
        main_title *= " (Step $time_step)"
    end
    
    summary_plot = plot(p1, p2, p3, p4, 
                        layout=(2, 2), size=(800, 600),
                        plot_title=main_title)
    
    return summary_plot
end

"""
    animate_simulation(W, xy, xy_s, G, cfg; output_file="simulation.gif")

Create an animated GIF of the simulation evolution.
"""
function animate_simulation(W::AbstractMatrix, xy::AbstractMatrix, xy_s::AbstractMatrix,
                          G::GeometryData, cfg; output_file="simulation.gif", 
                          frame_skip=10)
    
    n_velocity_nodes = size(xy, 1)
    n_time_steps = size(W, 2)
    
    # Create animation
    anim = @animate for i in 1:frame_skip:n_time_steps
        u = W[1:n_velocity_nodes, i]
        v = W[n_velocity_nodes+1:end, i]
        
        # Simple pressure approximation (zero for this demo)
        p = zeros(size(xy_s, 1))
        
        t = (i-1) * cfg.time_step
        title = @sprintf("Velocity Field (t = %.3f)", t)
        
        plot_velocity_magnitude(xy, u, v, title=title)
    end
    
    gif(anim, output_file, fps=10)
    println("Animation saved to $output_file")
    
    return anim
end
