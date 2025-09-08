classdef TestCylinderGeometry < BaseGeometryTest
    % TESTCYLINDERGEOMETRY Test cylinder geometry-specific functionality

    properties (Constant)
        GEOMETRY_TYPE = 'cylinder'
        CONFIG_FUNCTION = 'config_cylinder'
        EXPECTED_FIELDS = {'obstacle_radius'}
    end
end
