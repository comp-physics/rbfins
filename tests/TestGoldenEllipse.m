classdef TestGoldenEllipse < BaseGeometryTest
    % TESTGOLDENELLIPSE Test ellipse simulation against golden reference

    properties (Constant)
        GEOMETRY_TYPE = 'ellipse'
        EXPECTED_FIELDS = {'ellipse_a', 'ellipse_b'}
    end

    methods (Test)

        function testMatchesGolden(testCase)
            % Run the inherited golden test for ellipse
            testMatchesGolden@BaseGeometryTest(testCase);
        end

    end
end
