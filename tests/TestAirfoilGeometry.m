classdef TestAirfoilGeometry < BaseGeometryTest
    % TESTAIRFOILGEOMETRY Test NACA airfoil geometry-specific functionality

    properties (Constant)
        GEOMETRY_TYPE = 'airfoil'
        EXPECTED_FIELDS = {'naca_digits', 'chord_length', 'angle_of_attack', 'airfoil_x_center', 'airfoil_y_center'}
    end
end
