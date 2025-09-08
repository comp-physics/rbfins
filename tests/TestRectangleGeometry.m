classdef TestRectangleGeometry < BaseGeometryTest
    % TESTRECTANGLEGEOMETRY Test rectangle geometry-specific functionality

    properties (Constant)
        GEOMETRY_TYPE = 'rectangle'
        CONFIG_FUNCTION = 'config_rectangle'
        EXPECTED_FIELDS = {'rect_width', 'rect_height'}
    end
end
