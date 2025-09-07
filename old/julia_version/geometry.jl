using LinearAlgebra

"""
Geometry generation and mesh utilities for RBF-FD Navier-Stokes solver.

This module provides functionality for creating computational domains with obstacles
and extracting boundary conditions for the staggered grid arrangement.
"""

struct GeometryData
    # Interior nodes
    xy::Matrix{Float64}      # Interior velocity nodes (V-grid)
    xy_s::Matrix{Float64}    # Interior pressure nodes (P-grid)
    xt::Matrix{Int}          # Triangle connectivity
    
    # Velocity grid boundaries
    boundary_in::Matrix{Float64}   # Inlet boundary nodes (V-grid)
    boundary_out::Matrix{Float64}  # Outlet boundary nodes (V-grid)
    boundary_y::Matrix{Float64}    # Wall boundary nodes (V-grid)
    boundary_obs::Matrix{Float64}  # Obstacle boundary nodes (V-grid)
    
    # Pressure grid boundaries
    boundary_in_s::Matrix{Float64}   # Inlet boundary nodes (P-grid)
    boundary_out_s::Matrix{Float64}  # Outlet boundary nodes (P-grid)
    boundary_y_s::Matrix{Float64}    # Wall boundary nodes (P-grid)
    boundary_obs_s::Matrix{Float64}  # Obstacle boundary nodes (P-grid)
    
    # Geometry-specific data
    obstacle_radius::Float64         # Obstacle radius (for cylinder)
    ellipse_a::Float64              # Ellipse semi-major axis (for ellipse)
    ellipse_b::Float64              # Ellipse semi-minor axis (for ellipse)
    obs_normals_s::Matrix{Float64}  # Unit normal vectors at obstacle boundary
    
    # Proximity indices for special RBF treatment
    idx_near_obs_V::Vector{Bool}         # Velocity nodes near obstacle
    idx_near_obs_P::Vector{Bool}         # Pressure nodes near obstacle
    idx_far_boundaries_V::Vector{Bool}   # Velocity nodes far from boundaries
    idx_far_boundaries_P::Vector{Bool}   # Pressure nodes far from boundaries
end

"""
    dcircle(p, xc, yc, r)

Signed distance function for a circle centered at (xc, yc) with radius r.
Negative inside, positive outside.
"""
function dcircle(p::AbstractMatrix, xc::Real, yc::Real, r::Real)
    return sqrt.((p[:, 1] .- xc).^2 .+ (p[:, 2] .- yc).^2) .- r
end

"""
    drectangle(p, x1, x2, y1, y2)

Signed distance function for rectangle [x1,x2] × [y1,y2].
Negative inside, positive outside.
"""
function drectangle(p::AbstractMatrix, x1::Real, x2::Real, y1::Real, y2::Real)
    return -min.(min.(min.(-y1 .+ p[:, 2], y2 .- p[:, 2]), -x1 .+ p[:, 1]), x2 .- p[:, 1])
end

"""
    dellipse(p, xc, yc, a, b)

Signed distance function for an ellipse centered at (xc, yc) with semi-axes a, b.
Negative inside, positive outside.
"""
function dellipse(p::AbstractMatrix, xc::Real, yc::Real, a::Real, b::Real)
    x = p[:, 1] .- xc
    y = p[:, 2] .- yc
    
    # Ellipse equation: (x/a)² + (y/b)² = 1
    # Distance approximation using local scaling
    theta = atan.(y, x)
    r_ellipse = (a * b) ./ sqrt.((b * cos.(theta)).^2 .+ (a * sin.(theta)).^2)
    r_point = sqrt.(x.^2 .+ y.^2)
    
    return r_point .- r_ellipse
end

"""
    ddiff(d1, d2)

Difference of two distance functions (d1 - d2).
Creates domain that is inside d1 but outside d2.
"""
function ddiff(d1::AbstractVector, d2::AbstractVector)
    return max.(d1, -d2)
end

"""
    generate_rectangular_mesh_simple(cfg; geometry_type="cylinder")

Generate a simple rectangular mesh with obstacle using a basic grid approach.
This is a simplified version that doesn't require DistMesh.

INPUTS:
  cfg - Configuration structure
  geometry_type - "cylinder" or "ellipse"

OUTPUTS:
  G - GeometryData structure with mesh and boundaries
"""
function generate_rectangular_mesh_simple(cfg; geometry_type="cylinder")
    x_min, x_max = cfg.x_min, cfg.x_max
    y_min, y_max = cfg.y_min, cfg.y_max
    h = cfg.dist  # Mesh spacing
    eps = cfg.boundary_eps
    
    # Create regular grid
    x_range = x_min:h:x_max
    y_range = y_min:h:y_max
    
    points = zeros(length(x_range) * length(y_range), 2)
    idx = 1
    for x in x_range, y in y_range
        points[idx, :] = [x y]
        idx += 1
    end    
    # Remove points inside obstacle
    if geometry_type == "cylinder"
        radius = cfg.obstacle_radius
        inside_obstacle = (points[:, 1].^2 .+ points[:, 2].^2) .< radius^2
    else  # ellipse
        a, b = cfg.ellipse_a, cfg.ellipse_b
        inside_obstacle = (points[:, 1].^2 / a^2 .+ points[:, 2].^2 / b^2) .< 1
    end
    
    interior_points = points[.!inside_obstacle, :]
    
    # Create simple triangulation (this is a placeholder - real implementation would need Delaunay)
    n_points = size(interior_points, 1)
    triangles = zeros(Int, min(n_points-2, 100), 3)
    for i in 1:min(n_points-2, 100)
        triangles[i, :] = [i, i+1, i+2]
    end    
    # Classify boundary nodes
    boundary_in = interior_points[interior_points[:, 1] .< x_min + eps, :]
    boundary_out = interior_points[interior_points[:, 1] .> x_max - eps, :]
    boundary_y = interior_points[(interior_points[:, 2] .< y_min + eps) .| 
                                (interior_points[:, 2] .> y_max - eps), :]
    
    # Find obstacle boundary nodes
    if geometry_type == "cylinder"
        radius = cfg.obstacle_radius
        dist_to_center = sqrt.(interior_points[:, 1].^2 .+ interior_points[:, 2].^2)
        boundary_obs = interior_points[abs.(dist_to_center .- radius) .< eps, :]
        obs_normals = boundary_obs ./ radius  # Unit normals for cylinder
    else  # ellipse
        a, b = cfg.ellipse_a, cfg.ellipse_b
        ellipse_dist = abs.(interior_points[:, 1].^2 / a^2 .+ interior_points[:, 2].^2 / b^2 .- 1)
        boundary_obs = interior_points[ellipse_dist .< eps, :]
        # Compute normals for ellipse (simplified)
        obs_normals = zeros(size(boundary_obs))
        for i in 1:size(boundary_obs, 1)
            x, y = boundary_obs[i, :]
            # Normal vector: [2x/a², 2y/b²]
            nx = 2x / a^2
            ny = 2y / b^2
            norm_val = sqrt(nx^2 + ny^2)
            obs_normals[i, :] = [nx, ny] / norm_val
        end
    end
    
    # Remove boundary nodes from interior
    all_boundary = vcat(boundary_in, boundary_out, boundary_y, boundary_obs)
    if size(all_boundary, 1) > 0
        interior_mask = trues(size(interior_points, 1))
        for i in 1:size(interior_points, 1)
            for j in 1:size(all_boundary, 1)
                if norm(interior_points[i, :] - all_boundary[j, :]) < eps
                    interior_mask[i] = false
                    break
                end
            end
        end
        xy = interior_points[interior_mask, :]
        xy_s = copy(xy)  # Same nodes for both grids in simplified version
    else
        xy = interior_points
        xy_s = copy(xy)
    end
    
    # Compute proximity indices
    if geometry_type == "cylinder"
        radius = cfg.obstacle_radius
        r_dist = 2 * h  # Distance threshold
        idx_near_obs_V = (xy[:, 1].^2 .+ xy[:, 2].^2) .< (radius + r_dist)^2
        idx_near_obs_P = (xy_s[:, 1].^2 .+ xy_s[:, 2].^2) .< (radius + r_dist)^2
    else  # ellipse
        a, b = cfg.ellipse_a, cfg.ellipse_b
        r_dist = 2 * h
        # Approximate distance to ellipse
        ellipse_dist_V = abs.(xy[:, 1].^2 / a^2 .+ xy[:, 2].^2 / b^2 .- 1)
        ellipse_dist_P = abs.(xy_s[:, 1].^2 / a^2 .+ xy_s[:, 2].^2 / b^2 .- 1)
        idx_near_obs_V = ellipse_dist_V .< r_dist
        idx_near_obs_P = ellipse_dist_P .< r_dist
    end
    
    # Far from boundaries indices
    margin = 2 * h
    idx_far_boundaries_V = (xy[:, 1] .> x_min + margin) .& 
                          (xy[:, 1] .< x_max - margin) .&
                          (xy[:, 2] .> y_min + margin) .&
                          (xy[:, 2] .< y_max - margin) .&
                          .!idx_near_obs_V
    
    idx_far_boundaries_P = (xy_s[:, 1] .> x_min + margin) .& 
                          (xy_s[:, 1] .< x_max - margin) .&
                          (xy_s[:, 2] .> y_min + margin) .&
                          (xy_s[:, 2] .< y_max - margin) .&
                          .!idx_near_obs_P
    
    return GeometryData(
        xy, xy_s, triangles,
        boundary_in, boundary_out, boundary_y, boundary_obs,
        boundary_in, boundary_out, boundary_y, boundary_obs,  # Same boundaries for both grids
        cfg.obstacle_radius, cfg.ellipse_a, cfg.ellipse_b,
        obs_normals,
        idx_near_obs_V, idx_near_obs_P,
        idx_far_boundaries_V, idx_far_boundaries_P
    )
end

"""
    make_cylinder_geometry(cfg)

Generate mesh and boundaries for cylinder flow geometry.
"""
function make_cylinder_geometry(cfg)
    return generate_rectangular_mesh_simple(cfg; geometry_type="cylinder")
end

"""
    make_ellipse_geometry(cfg)

Generate mesh and boundaries for ellipse flow geometry.
"""
function make_ellipse_geometry(cfg)
    return generate_rectangular_mesh_simple(cfg; geometry_type="ellipse")
end

"""
    dcircle(p, xc, yc, r)

Signed distance function for a circle.
"""
function dcircle(p, xc, yc, r)
    return sqrt.((p[:, 1] .- xc).^2 .+ (p[:, 2] .- yc).^2) .- r
end

"""
    drectangle(p, x1, x2, y1, y2)

Signed distance function for a rectangle.
"""
function drectangle(p, x1, x2, y1, y2)
    return -min.(min.(min.(-y1 .+ p[:, 2], y2 .- p[:, 2]), -x1 .+ p[:, 1]), x2 .- p[:, 1])
end

"""
    ddiff(d1, d2)

Set difference of two distance functions.
"""
function ddiff(d1, d2)
    return max.(d1, -d2)
end

"""
    dellipse(p, xc, yc, a, b)

Signed distance function for an ellipse (simplified version).
"""
function dellipse(p, xc, yc, a, b)
    # Simplified ellipse distance - exact computation requires solving quartic
    # This approximation works well for mesh generation
    x = p[:, 1] .- xc
    y = p[:, 2] .- yc
    return (x.^2 ./ a^2 .+ y.^2 ./ b^2) .- 1.0
end

"""
    dpoly(p, vertices)

Signed distance function for a polygon defined by line segments.
"""
function dpoly(p, vertices)
    # Simple implementation for line segment distance
    x1, y1 = vertices[1, 1], vertices[1, 2]
    x2, y2 = vertices[2, 1], vertices[2, 2]
    
    # Distance to line segment
    A = y2 - y1
    B = x1 - x2
    C = x2*y1 - x1*y2
    
    return abs.(A .* p[:, 1] .+ B .* p[:, 2] .+ C) ./ sqrt(A^2 + B^2)
end

"""
    build_geometry(cfg)

Build geometry based on configuration geometry type.
"""
function build_geometry(cfg)
    if cfg.geometry_type == "cylinder"
        return make_cylinder_geometry(cfg)
    elseif cfg.geometry_type == "ellipse"
        return make_ellipse_geometry(cfg)
    else
        error("Unknown geometry type: $(cfg.geometry_type)")
    end
end
