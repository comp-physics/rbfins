module RBFNS

using LinearAlgebra
using SparseArrays

include("config.jl")
include("nearest_interp.jl")
include("RBF_PHS_FD_all.jl")
include("precompute_velocity_solvers.jl")
include("build_pressure_system.jl")
include("NS_2d_fractional_step_PHS.jl")
include("geometry.jl")
include("plotting.jl")

# Re-export main functions
export nearest_interp, rbf_phs_fd_all, precompute_velocity_solvers, build_pressure_system
export ns2d_fractional_step_phs, default_config, Config
export GeometryData, build_geometry, make_cylinder_geometry, make_ellipse_geometry
export plot_mesh, plot_velocity_field, plot_velocity_magnitude, plot_pressure_field
export plot_simulation_summary, animate_simulation

end
