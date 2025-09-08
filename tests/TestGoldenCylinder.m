classdef TestGoldenCylinder < BaseGeometryTest
    % TESTGOLDENCYLINDER Test cylinder simulation against golden reference

    properties (Constant)
        GEOMETRY_TYPE = 'cylinder'
        CONFIG_FUNCTION = 'config_cylinder'
        EXPECTED_FIELDS = {'obstacle_radius'}
    end

    methods (Test)
        function testMatchesGolden(testCase)
            % Run the inherited golden test for cylinder
            testMatchesGolden@BaseGeometryTest(testCase);
        end
    end
end
