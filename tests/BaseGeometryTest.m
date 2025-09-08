classdef BaseGeometryTest < matlab.unittest.TestCase
    % BASEGEOMETRYTEST Base class for geometry-specific tests
    %
    % This class provides common functionality for testing different geometries
    % (cylinder, ellipse, rectangle) to avoid code duplication.

    properties (Abstract, Constant)
        GEOMETRY_TYPE     % e.g., 'cylinder', 'ellipse', 'rectangle'
        EXPECTED_FIELDS   % Cell array of required geometry-specific config fields
    end

    properties (Constant)
        % Golden test tolerances
        XY_TOL = 5e-3       % Tolerance for node coordinate matching
        REL_TOL = 5e-3      % Relative tolerance for velocity field comparison
        ABS_TOL = 5e-3      % Absolute tolerance for velocity field comparison
    end

    methods

        function goldenFile = getGoldenFilePath(testCase)
            % Get the path to the golden reference file for this geometry
            goldenFile = fullfile(fileparts(mfilename('fullpath')), 'golden', ...
                                  sprintf('%s_Re100_Nt20_dt0.01_seed42.mat', testCase.GEOMETRY_TYPE));
        end

    end

    methods (Test)

        function smokeRun(testCase)
            % Test that the geometry simulation runs without errors in CI environment

            % Set CI environment variable to ensure short run
            setenv('CI', 'true');
            setenv('MATLAB_TEST', 'true');

            % Set up paths
            setup_paths();

            % Test that the script runs successfully
            try
                % Run the script and capture any errors
                evalin('base', 'clear; success = false; errorMsg = "";');
                evalin('base', sprintf(['try; setenv(''CONFIG_FUNC'', ''config''); ' ...
                                        'setenv(''GEOMETRY_TYPE'', ''%s''); simulate; success = true; ' ...
                                        'catch ME; errorMsg = ME.message; end;'], testCase.GEOMETRY_TYPE));

                % Check if the script ran successfully
                success = evalin('base', 'success');

                if success
                    fprintf('✅ %s simulation executed successfully without errors!\n', ...
                            testCase.GEOMETRY_TYPE);
                    testCase.verifyTrue(true, sprintf('%s simulation completed successfully', ...
                                                      testCase.GEOMETRY_TYPE));
                else
                    errorMsg = evalin('base', 'errorMsg');
                    fprintf('❌ %s simulation failed with error: %s\n', ...
                            testCase.GEOMETRY_TYPE, errorMsg);
                    testCase.verifyTrue(false, sprintf('%s simulation failed: %s', ...
                                                       testCase.GEOMETRY_TYPE, errorMsg));
                end

            catch ME
                fprintf('❌ Test execution failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Test execution failed: %s', ME.message));
            end
        end

        function testConfig(testCase)
            % Test that the geometry configuration loads correctly
            try
                setup_paths();

                % Test that config function can be called
                cfg = config(testCase.GEOMETRY_TYPE);
                testCase.verifyTrue(isstruct(cfg), sprintf('config(''%s'') should return a structure', ...
                                                           testCase.GEOMETRY_TYPE));

                % Check geometry type
                testCase.verifyEqual(cfg.geometry.type, testCase.GEOMETRY_TYPE, ...
                                     sprintf('Geometry type should be %s', testCase.GEOMETRY_TYPE));

                % Check for geometry-specific fields
                for i = 1:length(testCase.EXPECTED_FIELDS)
                    field = testCase.EXPECTED_FIELDS{i};
                    testCase.verifyTrue(isfield(cfg.geometry, field), ...
                                        sprintf('Config should have %s field', field));

                    % Verify field is positive if it's a numeric dimension
                    if contains(field, {'radius', 'width', 'height', '_a', '_b'})
                        testCase.verifyTrue(cfg.geometry.(field) > 0, ...
                                            sprintf('%s should be positive', field));
                    end
                end

                fprintf('✅ %s configuration test passed!\n', testCase.GEOMETRY_TYPE);

            catch ME
                fprintf('❌ %s configuration test failed: %s\n', testCase.GEOMETRY_TYPE, ME.message);
                testCase.verifyTrue(false, sprintf('%s configuration test failed: %s', ...
                                                   testCase.GEOMETRY_TYPE, ME.message));
            end
        end

        function testMatchesGolden(testCase)
            % Test that current simulation matches the golden reference

            fprintf('=== Testing %s Against Golden Reference ===\n', ...
                    upper(testCase.GEOMETRY_TYPE));

            % Set up environment
            setenv('CI', 'true');
            setenv('MATLAB_TEST', 'true');
            setup_paths();

            % Check if golden file exists
            goldenFile = testCase.getGoldenFilePath();
            if ~exist(goldenFile, 'file')
                testCase.assumeFail(sprintf( ...
                                            'Golden file not found: %s\nGenerate it with: make_golden_%s()', ...
                                            goldenFile, testCase.GEOMETRY_TYPE));
            end

            % Load golden reference data
            fprintf('Loading golden reference from: %s\n', goldenFile);
            S = load(goldenFile, 'gold');
            gold = S.gold;

            fprintf('Golden reference info:\n');
            fprintf('  Algorithm: %s\n', gold.meta.algorithm);
            fprintf('  Geometry: %s\n', gold.meta.geometry_type);
            fprintf('  Reynolds: %d\n', gold.meta.reynolds_number);
            fprintf('  Time steps: %d\n', gold.meta.num_time_steps);
            fprintf('  Nodes: %d\n', length(gold.U));
            fprintf('  Created: %s\n', gold.meta.timestamp);

            % Run current simulation
            fprintf('\nRunning current %s simulation...\n', testCase.GEOMETRY_TYPE);

            try
                evalin('base', 'clear; success = false; errorMsg = "";');
                evalin('base', sprintf(['try; setenv(''CONFIG_FUNC'', ''config''); ' ...
                                        'setenv(''GEOMETRY_TYPE'', ''%s''); simulate; success = true; ' ...
                                        'catch ME; errorMsg = ME.message; end;'], testCase.GEOMETRY_TYPE));

                success = evalin('base', 'success');
                testCase.verifyTrue(logical(success), ...
                                    sprintf('Current %s simulation should complete without errors', ...
                                            testCase.GEOMETRY_TYPE));

                % Check if variables exist before fetching them
                xy1_exists = evalin('base', 'exist(''xy1'', ''var'') > 0');
                W_exists = evalin('base', 'exist(''W'', ''var'') > 0');

                testCase.verifyTrue(xy1_exists, 'xy1 variable should exist after simulation');
                testCase.verifyTrue(W_exists, 'W variable should exist after simulation');

                if xy1_exists && W_exists
                    % Fetch current results
                    xy1 = evalin('base', 'xy1');
                    W = evalin('base', 'W');
                else
                    testCase.verifyTrue(false, 'Required simulation variables not found');
                    return
                end
            catch ME
                testCase.verifyTrue(false, sprintf('Failed to run simulation: %s', ME.message));
                return
            end

            n = size(xy1, 1);
            U = W(1:n, end);
            V = W(n + 1:end, end);

            fprintf('Current simulation completed:\n');
            fprintf('  Nodes: %d\n', n);
            fprintf('  Max |U|: %.6f\n', max(abs(U)));
            fprintf('  Max |V|: %.6f\n', max(abs(V)));

            % Verify we have the same number of nodes
            testCase.verifyEqual(length(U), length(gold.U), ...
                                 'Current and golden simulations must have same number of nodes');

            % Map current results to golden node ordering
            fprintf('\nMapping current nodes to golden reference...\n');
            [U_mapped, V_mapped] = BaseGeometryTest.mapToGolden(xy1, U, V, gold.xy1, testCase.XY_TOL);

            % Compute differences
            [relU, absU] = BaseGeometryTest.relAbsDiff(U_mapped, gold.U);
            [relV, absV] = BaseGeometryTest.relAbsDiff(V_mapped, gold.V);

            fprintf('\nComparison results:\n');
            fprintf('  U field - Rel diff: %.2e, Abs diff: %.2e\n', relU, absU);
            fprintf('  V field - Rel diff: %.2e, Abs diff: %.2e\n', relV, absV);
            fprintf('  Tolerances - Rel: %.2e, Abs: %.2e\n', testCase.REL_TOL, testCase.ABS_TOL);

            % Verify tolerances
            testCase.verifyLessThanOrEqual(relU, testCase.REL_TOL, ...
                                           sprintf('U field relative difference %.2e exceeds tolerance %.2e', ...
                                                   relU, testCase.REL_TOL));
            testCase.verifyLessThanOrEqual(relV, testCase.REL_TOL, ...
                                           sprintf('V field relative difference %.2e exceeds tolerance %.2e', ...
                                                   relV, testCase.REL_TOL));
            testCase.verifyLessThanOrEqual(absU, testCase.ABS_TOL, ...
                                           sprintf('U field absolute difference %.2e exceeds tolerance %.2e', ...
                                                   absU, testCase.ABS_TOL));
            testCase.verifyLessThanOrEqual(absV, testCase.ABS_TOL, ...
                                           sprintf('V field absolute difference %.2e exceeds tolerance %.2e', ...
                                                   absV, testCase.ABS_TOL));

            fprintf('\n[PASS] %s golden file test passed!\n', testCase.GEOMETRY_TYPE);
        end

    end

    methods (Static)

        function [U_mapped, V_mapped] = mapToGolden(xy_current, U_current, V_current, xy_golden, xy_tol)
            % Map current simulation results to golden reference node ordering

            % Sort current data lexicographically to match golden ordering
            [~, idx_sorted] = sortrows(xy_current, [1, 2]);
            xy_sorted = xy_current(idx_sorted, :);
            U_sorted = U_current(idx_sorted);
            V_sorted = V_current(idx_sorted);

            % Verify that sorted coordinates match golden coordinates
            coord_diff = xy_sorted - xy_golden;
            max_coord_diff = max(sqrt(sum(coord_diff.^2, 2)));

            if max_coord_diff > xy_tol
                error('Node mapping failed: maximum coordinate difference %.2e exceeds tolerance %.2e', ...
                      max_coord_diff, xy_tol);
            end

            % If coordinates match within tolerance, use sorted ordering
            U_mapped = U_sorted;
            V_mapped = V_sorted;
        end

        function [relErr, absErr] = relAbsDiff(a, b)
            % Compute relative and absolute differences between two arrays

            diff = a - b;
            absErr = max(abs(diff));

            % Use norm of reference for relative error, with small epsilon to avoid division by zero
            denom = max(1e-12, norm(b));
            relErr = norm(diff) / denom;
        end

    end
end
