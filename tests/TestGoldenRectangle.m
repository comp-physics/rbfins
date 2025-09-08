classdef TestGoldenRectangle < BaseGeometryTest
    % TESTGOLDENRECTANGLE Test rectangle simulation against golden reference

    properties (Constant)
        GEOMETRY_TYPE = 'rectangle'
        CONFIG_FUNCTION = 'config_rectangle'
        EXPECTED_FIELDS = {'rect_width', 'rect_height'}
    end

    methods (Test)
        function testMatchesGolden(testCase)
            % Run the inherited golden test for rectangle
            testMatchesGolden@BaseGeometryTest(testCase);
        end
    end
end
