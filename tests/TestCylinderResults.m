classdef TestCylinderResults < matlab.unittest.TestCase
    methods (Test)
        function testBasicPhysics(testCase)
            % A simplified test that just verifies the simulation runs and produces reasonable results
            try
                % Set up environment
                setenv('CI', 'true');
                setenv('MATLAB_TEST', 'true');
                setup_paths();

                % Run the simulation in a controlled way
                evalin('base', 'clear; success = false; errorMsg = "";');
                evalin('base', 'try; simulate; success = true; catch ME; errorMsg = ME.message; end;');

                % Check if the simulation ran successfully
                success = evalin('base', 'success');
                testCase.verifyTrue(logical(success), 'Simulation should complete without errors');

                % Basic verification of simulation results
                % Check that W exists and has reasonable size
                W_exists = evalin('base', 'exist(''W'', ''var'') > 0');
                testCase.verifyTrue(W_exists, 'W variable should exist after simulation');

                if W_exists
                    % Get the size of W
                    W_size = evalin('base', 'size(W)');
                    testCase.verifyTrue(W_size(1) > 0, 'W should have positive number of rows');
                    testCase.verifyTrue(W_size(2) > 0, 'W should have positive number of columns');

                    % Check for NaN or Inf values in the final state
                    has_nan = evalin('base', 'any(isnan(W(:,end)))');
                    has_inf = evalin('base', 'any(isinf(W(:,end)))');
                    testCase.verifyFalse(has_nan, 'Final state should not contain NaN values');
                    testCase.verifyFalse(has_inf, 'Final state should not contain Inf values');

                    % Check that velocity components have reasonable values

                    % Extract U and V from the final state
                    U = evalin('base', 'W(1:floor(length(W)/2),end)');
                    V = evalin('base', 'W(floor(length(W)/2)+1:end,end)');

                    % Check velocity magnitudes
                    max_U = max(abs(U));
                    max_V = max(abs(V));
                    testCase.verifyLessThan(max_U, 10.0, 'Maximum U velocity should be reasonable');
                    testCase.verifyLessThan(max_V, 10.0, 'Maximum V velocity should be reasonable');

                    fprintf('[PASS] Basic physics test passed!\n');
                    fprintf('   Max |U|: %.3f, Max |V|: %.3f\n', max_U, max_V);
                end

            catch ME
                fprintf('[FAIL] Basic physics test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Basic physics test failed: %s', ME.message));
            end
        end

        function testVelocityDistribution(testCase)
            % A simplified test that checks basic velocity distribution properties
            try
                % Set up environment
                setenv('CI', 'true');
                setenv('MATLAB_TEST', 'true');
                setup_paths();

                % Run the simulation if it hasn't been run already
                W_exists = evalin('base', 'exist(''W'', ''var'') > 0');
                if ~W_exists
                    evalin('base', 'try; simulate; catch; end;');
                end

                % Check if W exists now
                W_exists = evalin('base', 'exist(''W'', ''var'') > 0');
                testCase.verifyTrue(W_exists, 'W variable should exist after simulation');

                if W_exists
                    % Extract U and V from the final state
                    U = evalin('base', 'W(1:floor(length(W)/2),end)');
                    V = evalin('base', 'W(floor(length(W)/2)+1:end,end)');

                    % Test 1: Most of the domain should have positive U velocity (flow to the right)
                    positive_U_fraction = sum(U > 0) / length(U);
                    testCase.verifyGreaterThan(positive_U_fraction, 0.6, 'Most of domain should have positive U velocity');

                    % Test 2: V velocities should be roughly symmetric (mean close to 0)
                    mean_V = mean(V);
                    testCase.verifyLessThan(abs(mean_V), 0.5, 'V velocities should be roughly symmetric');

                    fprintf('[PASS] Velocity distribution test passed!\n');
                    fprintf('   Positive U fraction: %.3f\n', positive_U_fraction);
                    fprintf('   Mean V: %.3f\n', mean_V);
                end

            catch ME
                fprintf('❌ Velocity distribution test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Velocity distribution test failed: %s', ME.message));
            end
        end

        function testEnergyBalance(testCase)
            % A simplified test that checks total kinetic energy is reasonable
            try
                % Set up environment
                setenv('CI', 'true');
                setenv('MATLAB_TEST', 'true');
                setup_paths();

                % Run the simulation if it hasn't been run already
                W_exists = evalin('base', 'exist(''W'', ''var'') > 0');
                if ~W_exists
                    evalin('base', 'try; simulate; catch; end;');
                end

                % Check if W exists now
                W_exists = evalin('base', 'exist(''W'', ''var'') > 0');
                testCase.verifyTrue(W_exists, 'W variable should exist after simulation');

                if W_exists
                    % Extract U and V from the final state
                    U = evalin('base', 'W(1:floor(length(W)/2),end)');
                    V = evalin('base', 'W(floor(length(W)/2)+1:end,end)');

                    % Calculate total kinetic energy
                    kinetic_energy = 0.5 * sum(U.^2+V.^2);
                    testCase.verifyGreaterThan(kinetic_energy, 0, 'Total kinetic energy should be positive');
                    testCase.verifyLessThan(kinetic_energy, 10000, 'Total kinetic energy should be reasonable');

                    % Calculate energy density
                    energy_density = kinetic_energy / length(U);
                    testCase.verifyLessThan(energy_density, 5.0, 'Energy density should be reasonable');

                    fprintf('✅ Energy balance test passed!\n');
                    fprintf('   Total kinetic energy: %.3f\n', kinetic_energy);
                    fprintf('   Energy density: %.4f\n', energy_density);
                end

            catch ME
                fprintf('❌ Energy balance test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Energy balance test failed: %s', ME.message));
            end
        end
    end
end

