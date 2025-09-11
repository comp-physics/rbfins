classdef TestGoldenAirfoil < BaseGeometryTest
  % TESTGOLDENAIRFOIL Test NACA airfoil simulation against golden reference

  properties (Constant)
    GEOMETRY_TYPE = 'airfoil'
    EXPECTED_FIELDS = {'naca_digits', 'chord_length', 'angle_of_attack', 'airfoil_x_center', 'airfoil_y_center'}
  end

  methods (Test)

    function testMatchesGolden(testCase)
      % Run the inherited golden test for airfoil
      testMatchesGolden@BaseGeometryTest(testCase);
    end

  end
end
