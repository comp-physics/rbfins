function G = build_geometry(cfg)
%BUILD_GEOMETRY Generate geometry using appropriate geometry helper
%
% This function dispatches to the appropriate geometry generation function
% based on the configuration and prints a summary of the selected geometry.
%
% INPUTS:
%   cfg - Configuration structure from config()
%
% OUTPUTS:
%   G - Geometry structure containing mesh, boundaries, and geometry-specific data

% Generate geometry using appropriate geometry helper
% This supports different obstacle geometries via the config.geometry.type parameter
switch lower(cfg.geometry.type)
    case 'cylinder'
        fprintf('Using cylinder geometry (radius=%.2f)...\n', cfg.geometry.obstacle_radius);
        G = make_cylinder_geometry(cfg);
    case 'ellipse'
        fprintf('Using ellipse geometry (a=%.2f, b=%.2f)...\n', cfg.geometry.ellipse_a, cfg.geometry.ellipse_b);
        G = make_ellipse_geometry(cfg);
    otherwise
        error('Unknown geometry type: %s. Supported types: cylinder, ellipse', cfg.geometry.type);
end

end
