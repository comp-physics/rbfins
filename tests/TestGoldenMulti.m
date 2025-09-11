classdef TestGoldenMulti < BaseGeometryTest
  %TESTGOLDENMULTI Test class for multi-obstacle geometry golden file validation
  %
  % This class tests the multi-obstacle geometry implementation against
  % a golden reference file to ensure numerical consistency and prevent
  % regressions in the simulation results.

  properties (Constant)
    GEOMETRY_TYPE = 'multi'
    EXPECTED_FIELDS = {'obstacles'}
  end

  methods

    function goldenFile = getGoldenFilePath(~)
      % Override to use multi-specific golden file path with correct time step
      goldenFile = fullfile(fileparts(mfilename('fullpath')), 'golden', ...
                            'multi_Re100_Nt20_dt0.005_seed42.mat');
    end

  end

end
