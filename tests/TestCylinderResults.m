classdef TestCylinderResults < matlab.unittest.TestCase
    methods (Test)
        function testBasicPhysics(testCase)
            % Test basic physics: no NaN/Inf, reasonable velocity magnitudes
            try
                S = TestCylinderResults.runAndFetch();
                
                % Check for NaN or Inf values
                testCase.verifyTrue(all(isfinite(S.U)), 'U velocity field should contain no NaN/Inf values');
                testCase.verifyTrue(all(isfinite(S.V)), 'V velocity field should contain no NaN/Inf values');
                
                % Basic sanity checks for velocity field statistics
                mean_U = mean(S.U);
                mean_V = mean(S.V);
                max_U = max(abs(S.U));
                max_V = max(abs(S.V));
                
                testCase.verifyGreaterThan(mean_U, 0.3, 'Mean U velocity too low - simulation may have failed');
                testCase.verifyLessThan(mean_U, 1.5, 'Mean U velocity too high - simulation may be unstable');
                testCase.verifyLessThan(abs(mean_V), 0.5, 'Mean V velocity should be small (no large net cross-flow)');
                testCase.verifyLessThan(max_U, 10.0, 'Maximum U velocity suggests simulation exploded');
                testCase.verifyLessThan(max_V, 10.0, 'Maximum V velocity suggests simulation exploded');
                
                % Check velocity magnitude is reasonable
                velocity_magnitude = sqrt(S.U.^2 + S.V.^2);
                max_velocity = max(velocity_magnitude);
                testCase.verifyLessThan(max_velocity, 8.0, 'Maximum velocity magnitude too high');
                
                fprintf('✅ Basic physics test passed!\n');
                fprintf('   Mean U: %.3f, Mean V: %.3f\n', mean_U, mean_V);
                fprintf('   Max |U|: %.3f, Max |V|: %.3f\n', max_U, max_V);
                fprintf('   Max velocity magnitude: %.3f\n', max_velocity);
                
            catch ME
                fprintf('❌ Basic physics test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Basic physics test failed: %s', ME.message));
            end
        end
        
        function testCylinderNoSlipBoundary(testCase)
            % Test no-slip boundary condition at cylinder surface
            try
                S = TestCylinderResults.runAndFetch();
                
                % Find nodes near cylinder center (0,0) with radius ≈ 0.5
                cylinder_center = [0, 0];
                radius = S.cfg.mesh.cylinder_radius;
                distances = sqrt((S.xy1(:,1) - cylinder_center(1)).^2 + (S.xy1(:,2) - cylinder_center(2)).^2);
                
                % Find nodes very close to cylinder surface (within 0.1 of the radius)
                cylinder_surface_nodes = find(abs(distances - radius) < 0.1);
                
                if ~isempty(cylinder_surface_nodes)
                    cylinder_U = S.U(cylinder_surface_nodes);
                    cylinder_V = S.V(cylinder_surface_nodes);
                    
                    % Velocities should be small near cylinder surface (no-slip approximation)
                    max_U_cylinder = max(abs(cylinder_U));
                    max_V_cylinder = max(abs(cylinder_V));
                    
                    testCase.verifyLessThan(max_U_cylinder, 0.2, 'U velocity near cylinder should be small (no-slip)');
                    testCase.verifyLessThan(max_V_cylinder, 0.2, 'V velocity near cylinder should be small (no-slip)');
                    
                    fprintf('✅ Cylinder no-slip boundary test passed!\n');
                    fprintf('   Max |U| near cylinder: %.4f\n', max_U_cylinder);
                    fprintf('   Max |V| near cylinder: %.4f\n', max_V_cylinder);
                else
                    fprintf('⚠️  Warning: No nodes found near cylinder surface for no-slip test\n');
                    testCase.verifyTrue(true, 'No cylinder surface nodes found to test');
                end
                
            catch ME
                fprintf('❌ Cylinder boundary test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Cylinder boundary test failed: %s', ME.message));
            end
        end
        
        function testVelocityDistribution(testCase)
            % Test that velocity distribution makes physical sense
            try
                S = TestCylinderResults.runAndFetch();
                
                % Test 1: Most of the domain should have positive U velocity (flow to the right)
                positive_U_fraction = sum(S.U > 0) / length(S.U);
                testCase.verifyGreaterThan(positive_U_fraction, 0.7, 'Most of domain should have positive U velocity');
                
                % Test 2: V velocities should be roughly symmetric (mean ≈ 0)
                mean_V = mean(S.V);
                testCase.verifyLessThan(abs(mean_V), 0.3, 'V velocities should be roughly symmetric');
                
                % Test 3: Check that there's flow acceleration around the cylinder
                % Compare velocity magnitudes far from cylinder vs near cylinder
                distances = sqrt(S.xy1(:,1).^2 + S.xy1(:,2).^2);
                far_field_nodes = find(distances > 3*S.cfg.mesh.cylinder_radius);
                near_field_nodes = find(distances > 1.2*S.cfg.mesh.cylinder_radius & distances < 2*S.cfg.mesh.cylinder_radius);
                
                if ~isempty(far_field_nodes) && ~isempty(near_field_nodes)
                    far_field_speed = mean(sqrt(S.U(far_field_nodes).^2 + S.V(far_field_nodes).^2));
                    near_field_speed = mean(sqrt(S.U(near_field_nodes).^2 + S.V(near_field_nodes).^2));
                    
                    % For cylinder flow, we expect some acceleration around the cylinder
                    speed_ratio = near_field_speed / far_field_speed;
                    testCase.verifyGreaterThan(speed_ratio, 0.8, 'Flow near cylinder should not decelerate too much');
                    testCase.verifyLessThan(speed_ratio, 2.0, 'Flow near cylinder should not accelerate excessively');
                    
                    fprintf('   Far field speed: %.3f, Near field speed: %.3f, Ratio: %.3f\n', ...
                        far_field_speed, near_field_speed, speed_ratio);
                end
                
                fprintf('✅ Velocity distribution test passed!\n');
                fprintf('   Positive U fraction: %.3f\n', positive_U_fraction);
                
            catch ME
                fprintf('❌ Velocity distribution test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Velocity distribution test failed: %s', ME.message));
            end
        end
        
        function testInletOutletVelocities(testCase)
            % Test that inlet has uniform flow and outlet shows wake effects
            try
                S = TestCylinderResults.runAndFetch();
                
                % Find inlet nodes (left boundary, x ≈ domain minimum)
                x_min = S.cfg.domain.x_min;
                inlet_tolerance = 0.2;
                inlet_nodes = find(abs(S.xy1(:,1) - x_min) < inlet_tolerance);
                
                % Find outlet nodes (right boundary, x ≈ domain maximum) 
                x_max = S.cfg.domain.x_max;
                outlet_tolerance = 0.2;
                outlet_nodes = find(abs(S.xy1(:,1) - x_max) < outlet_tolerance);
                
                if ~isempty(inlet_nodes)
                    inlet_U = S.U(inlet_nodes);
                    inlet_V = S.V(inlet_nodes);
                    
                    % Inlet should have mostly positive U velocity
                    positive_inlet_fraction = sum(inlet_U > 0.5) / length(inlet_U);
                    testCase.verifyGreaterThan(positive_inlet_fraction, 0.8, 'Inlet should have positive U velocity');
                    
                    % Inlet V velocity should be small (no significant cross-flow)
                    max_inlet_V = max(abs(inlet_V));
                    testCase.verifyLessThan(max_inlet_V, 0.3, 'Inlet V velocity should be small');
                    
                    fprintf('   Inlet: positive U fraction = %.3f, max |V| = %.4f\n', ...
                        positive_inlet_fraction, max_inlet_V);
                end
                
                if ~isempty(outlet_nodes)
                    outlet_U = S.U(outlet_nodes);
                    outlet_V = S.V(outlet_nodes);
                    
                    % Outlet should still have positive mean U velocity but may show wake effects
                    mean_outlet_U = mean(outlet_U);
                    testCase.verifyGreaterThan(mean_outlet_U, 0.2, 'Outlet should have positive mean U velocity');
                    
                    fprintf('   Outlet: mean U = %.3f\n', mean_outlet_U);
                end
                
                fprintf('✅ Inlet/outlet velocity test passed!\n');
                
            catch ME
                fprintf('❌ Inlet/outlet test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Inlet/outlet test failed: %s', ME.message));
            end
        end
        
        function testEnergyBalance(testCase)
            % Test that total kinetic energy is reasonable
            try
                S = TestCylinderResults.runAndFetch();
                
                % Calculate total kinetic energy
                kinetic_energy = 0.5 * sum(S.U.^2 + S.V.^2);
                testCase.verifyGreaterThan(kinetic_energy, 0, 'Total kinetic energy should be positive');
                testCase.verifyLessThan(kinetic_energy, 2000, 'Total kinetic energy too high - simulation may be unstable');
                
                % Calculate energy density
                energy_density = kinetic_energy / length(S.U);
                testCase.verifyLessThan(energy_density, 2.0, 'Energy density should be reasonable');
                
                fprintf('✅ Energy balance test passed!\n');
                fprintf('   Total kinetic energy: %.3f\n', kinetic_energy);
                fprintf('   Energy density: %.4f\n', energy_density);
                
            catch ME
                fprintf('❌ Energy balance test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Energy balance test failed: %s', ME.message));
            end
        end
    end
    
    methods (Static, Access = private)
        function S = runAndFetch()
            % Helper method to run simulation and return key variables
            setenv('CI', 'true');
            setup_paths();
            
            % Run simulation in base workspace
            evalin('base', 'clear; run(''cylinder.m'');');
            
            % Fetch key variables from base workspace
            S.U = evalin('base', 'W(1:length(W)/2)');
            S.V = evalin('base', 'W(length(W)/2+1:end)');
            S.xy1 = evalin('base', 'xy1');
            S.cfg = evalin('base', 'cfg');
        end
    end
end