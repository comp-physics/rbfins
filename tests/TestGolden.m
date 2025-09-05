classdef TestGolden < matlab.unittest.TestCase
    % TESTGOLDEN Test current simulation results against golden reference files
    %
    % This test class compares the current simulation output with a saved
    % "golden" reference file to detect regressions in the numerical algorithm.
    %
    % The test performs the following:
    % 1. Loads a golden reference file
    % 2. Runs the current simulation
    % 3. Maps current nodes to golden nodes using nearest neighbor matching
    % 4. Compares velocity fields within specified tolerances
    %
    % Tolerances can be adjusted based on expected numerical precision
    % and platform-specific differences.
    
    properties (Constant)
        % Tolerance for node coordinate matching (geometric space)
        XY_TOL = 5e-3;
        
        % Relative tolerance for velocity field comparison
        REL_TOL = 5e-3;
        
        % Absolute tolerance for velocity field comparison (fallback)
        ABS_TOL = 5e-3;
        
        % Golden reference file path (relative to script location)
        GOLDEN_FILE = fullfile(fileparts(mfilename('fullpath')), 'golden', 'cylinder_Re100_Nt20_dt0.01_seed42.mat');
    end
    
    methods (Test)
        function testMatchesGolden(testCase)
            % Test that current simulation matches the golden reference
            
            fprintf('=== Testing Against Golden Reference ===\n');
            
            % Set up environment
            setenv('CI', 'true');
            setenv('MATLAB_TEST', 'true');
            setup_paths();
            
            % Check if golden file exists
            if ~exist(testCase.GOLDEN_FILE, 'file')
                testCase.assumeFail(sprintf(...
                    'Golden file not found: %s\nGenerate it with: make_golden_simple()', ...
                    testCase.GOLDEN_FILE));
            end
            
            % Load golden reference data
            fprintf('Loading golden reference from: %s\n', testCase.GOLDEN_FILE);
            S = load(testCase.GOLDEN_FILE, 'gold');
            gold = S.gold;
            
            fprintf('Golden reference info:\n');
            fprintf('  Algorithm: %s\n', gold.meta.algorithm);
            fprintf('  Reynolds: %d\n', gold.meta.reynolds_number);
            fprintf('  Time steps: %d\n', gold.meta.num_time_steps);
            fprintf('  Nodes: %d\n', length(gold.U));
            fprintf('  Created: %s\n', gold.meta.timestamp);
            
            % Run current simulation
            fprintf('\nRunning current simulation...\n');
            
            % Run simulation using same approach as other tests
            try
                evalin('base', 'clear; success = false; errorMsg = "";');
                evalin('base', 'try; simulate; success = true; catch ME; errorMsg = ME.message; end;');
                
                success = evalin('base', 'success');
                testCase.verifyTrue(logical(success), 'Current simulation should complete without errors');
                
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
                    return;
                end
            catch ME
                testCase.verifyTrue(false, sprintf('Failed to run simulation: %s', ME.message));
                return;
            end
            
            n = size(xy1, 1);
            U = W(1:n, end);
            V = W(n+1:end, end);
            
            fprintf('Current simulation completed:\n');
            fprintf('  Nodes: %d\n', n);
            fprintf('  Max |U|: %.6f\n', max(abs(U)));
            fprintf('  Max |V|: %.6f\n', max(abs(V)));
            
            % Verify we have the same number of nodes
            testCase.verifyEqual(length(U), length(gold.U), ...
                'Current and golden simulations must have same number of nodes');
            
            % Map current results to golden node ordering
            fprintf('\nMapping current nodes to golden reference...\n');
            [U_mapped, V_mapped] = TestGolden.mapToGolden(xy1, U, V, gold.xy1, testCase.XY_TOL);
            
            % Compute differences
            [relU, absU] = TestGolden.relAbsDiff(U_mapped, gold.U);
            [relV, absV] = TestGolden.relAbsDiff(V_mapped, gold.V);
            
            fprintf('\nComparison results:\n');
            fprintf('  U field - Rel diff: %.2e, Abs diff: %.2e\n', relU, absU);
            fprintf('  V field - Rel diff: %.2e, Abs diff: %.2e\n', relV, absV);
            fprintf('  Tolerances - Rel: %.2e, Abs: %.2e\n', testCase.REL_TOL, testCase.ABS_TOL);
            
            % Verify tolerances
            testCase.verifyLessThanOrEqual(relU, testCase.REL_TOL, ...
                sprintf('U field relative difference %.2e exceeds tolerance %.2e', relU, testCase.REL_TOL));
            testCase.verifyLessThanOrEqual(relV, testCase.REL_TOL, ...
                sprintf('V field relative difference %.2e exceeds tolerance %.2e', relV, testCase.REL_TOL));
            testCase.verifyLessThanOrEqual(absU, testCase.ABS_TOL, ...
                sprintf('U field absolute difference %.2e exceeds tolerance %.2e', absU, testCase.ABS_TOL));
            testCase.verifyLessThanOrEqual(absV, testCase.ABS_TOL, ...
                sprintf('V field absolute difference %.2e exceeds tolerance %.2e', absV, testCase.ABS_TOL));
            
            fprintf('\n[PASS] All golden file tests passed!\n');
        end
    end
    
    methods (Static)
        function [U_mapped, V_mapped] = mapToGolden(xy_current, U_current, V_current, xy_golden, xy_tol)
            % Map current simulation results to golden reference node ordering
            %
            % This function handles potential differences in node ordering between
            % the current simulation and the golden reference by using nearest
            % neighbor matching.
            
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
            %
            % relErr = ||a - b|| / max(||b||, eps)
            % absErr = ||a - b||_inf (infinity norm)
            
            diff = a - b;
            absErr = max(abs(diff));
            
            % Use norm of reference for relative error, with small epsilon to avoid division by zero
            denom = max(1e-12, norm(b));
            relErr = norm(diff) / denom;
        end
    end
end
