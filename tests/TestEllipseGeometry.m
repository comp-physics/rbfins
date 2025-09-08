classdef TestEllipseGeometry < BaseGeometryTest
    % TESTELLIPSEGEOMETRY Test ellipse geometry-specific functionality

    properties (Constant)
        GEOMETRY_TYPE = 'ellipse'
        CONFIG_FUNCTION = 'config_ellipse'
        EXPECTED_FIELDS = {'ellipse_a', 'ellipse_b'}
    end
end
