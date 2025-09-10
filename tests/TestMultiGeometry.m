classdef TestMultiGeometry < BaseGeometryTest
    %TESTMULTIGEOMETRY Test class for multi-obstacle geometry
    %
    % This class tests the multi-obstacle geometry implementation including:
    % - Configuration validation
    % - Geometry structure validation
    % - Boundary node classification
    % - Normal vector computation
    % - Smoke test execution

    properties (Constant)
        GEOMETRY_TYPE = 'multi'
        EXPECTED_FIELDS = {'obstacles'}  % Multi geometry has obstacles array instead of individual params
    end

    methods (Test)

        function testMultiConfig(testCase)
            % Test multi-obstacle configuration validation
            try
                setup_paths();

                % Test that config function can be called
                cfg = config(testCase.GEOMETRY_TYPE);
                testCase.verifyTrue(isstruct(cfg), 'config(''multi'') should return a structure');

                % Check geometry type
                testCase.verifyEqual(cfg.geometry.type, testCase.GEOMETRY_TYPE, ...
                                     'Geometry type should be multi');

                % Check obstacles array exists and is non-empty
                testCase.verifyTrue(isfield(cfg.geometry, 'obstacles'), ...
                                    'Config should have obstacles field');
                testCase.verifyTrue(~isempty(cfg.geometry.obstacles), ...
                                    'Obstacles array should not be empty');

                % Validate each obstacle structure
                obstacles = cfg.geometry.obstacles;
                for i = 1:length(obstacles)
                    obs = obstacles(i);
                    testCase.verifyTrue(isfield(obs, 'type'), ...
                                        sprintf('Obstacle %d should have type field', i));
                    testCase.verifyTrue(isfield(obs, 'center'), ...
                                        sprintf('Obstacle %d should have center field', i));
                    testCase.verifyTrue(isfield(obs, 'params'), ...
                                        sprintf('Obstacle %d should have params field', i));

                    % Validate center is 2D coordinate
                    testCase.verifyEqual(length(obs.center), 2, ...
                                         sprintf('Obstacle %d center should be [x, y]', i));

                    % Validate type-specific parameters
                    switch lower(obs.type)
                        case 'cylinder'
                            testCase.verifyTrue(isfield(obs.params, 'radius'), ...
                                                sprintf('Cylinder obstacle %d should have radius', i));
                            testCase.verifyTrue(obs.params.radius > 0, ...
                                                sprintf('Cylinder obstacle %d radius should be positive', i));
                        case 'ellipse'
                            testCase.verifyTrue(isfield(obs.params, 'a'), ...
                                                sprintf('Ellipse obstacle %d should have semi-major axis a', i));
                            testCase.verifyTrue(isfield(obs.params, 'b'), ...
                                                sprintf('Ellipse obstacle %d should have semi-minor axis b', i));
                            testCase.verifyTrue(obs.params.a > 0 && obs.params.b > 0, ...
                                                sprintf('Ellipse obstacle %d axes should be positive', i));
                        case 'rectangle'
                            testCase.verifyTrue(isfield(obs.params, 'width'), ...
                                                sprintf('Rectangle obstacle %d should have width', i));
                            testCase.verifyTrue(isfield(obs.params, 'height'), ...
                                                sprintf('Rectangle obstacle %d should have height', i));
                            testCase.verifyTrue(obs.params.width > 0 && obs.params.height > 0, ...
                                                sprintf('Rectangle obstacle %d dimensions should be positive', i));
                        case 'airfoil'
                            testCase.verifyTrue(isfield(obs.params, 'naca_digits'), ...
                                                sprintf('Airfoil obstacle %d should have naca_digits', i));
                            testCase.verifyTrue(isfield(obs.params, 'chord_length'), ...
                                                sprintf('Airfoil obstacle %d should have chord_length', i));
                            testCase.verifyTrue(obs.params.chord_length > 0, ...
                                                sprintf('Airfoil obstacle %d chord should be positive', i));
                        otherwise
                            testCase.verifyTrue(false, sprintf('Unknown obstacle type: %s', obs.type));
                    end
                end

                fprintf('[PASS] Multi-obstacle configuration test passed!\n');

            catch ME
                fprintf('[FAIL] Multi-obstacle configuration test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Multi-obstacle configuration test failed: %s', ME.message));
            end
        end

        function testGeometryStructure(testCase)
            % Test that multi-obstacle geometry structure is valid
            try
                setup_paths();

                % Load configuration and build geometry
                cfg = config(testCase.GEOMETRY_TYPE);
                G = build_geometry(cfg);

                % Verify basic geometry structure
                testCase.verifyTrue(isstruct(G), 'Geometry should be a structure');

                % Check required fields exist
                required_fields = {'xy', 'xy_s', 'xt', ...
                                   'boundary_in', 'boundary_out', 'boundary_y', 'boundary_obs', ...
                                   'boundary_in_s', 'boundary_out_s', 'boundary_y_s', 'boundary_obs_s', ...
                                   'obs_normals_s', 'idx_near_obs_V', 'idx_near_obs_P', ...
                                   'idx_far_boundaries_V', 'idx_far_boundaries_P'};

                for i = 1:length(required_fields)
                    field = required_fields{i};
                    testCase.verifyTrue(isfield(G, field), ...
                                        sprintf('Geometry should have %s field', field));
                end

                % Check multi-specific fields
                testCase.verifyTrue(isfield(G, 'obstacles'), 'Geometry should have obstacles field');
                testCase.verifyTrue(isfield(G, 'num_obstacles'), 'Geometry should have num_obstacles field');
                testCase.verifyTrue(isfield(G, 'boundary_obs_ids'), 'Geometry should have boundary_obs_ids field');

                % Verify node arrays are non-empty and have correct dimensions
                testCase.verifyTrue(~isempty(G.xy), 'Interior velocity nodes should not be empty');
                testCase.verifyTrue(~isempty(G.xy_s), 'Interior pressure nodes should not be empty');
                testCase.verifyEqual(size(G.xy, 2), 2, 'Velocity nodes should be 2D');
                testCase.verifyEqual(size(G.xy_s, 2), 2, 'Pressure nodes should be 2D');

                % Verify boundary nodes
                testCase.verifyTrue(~isempty(G.boundary_obs), 'Obstacle boundary nodes (V-grid) should not be empty');
                testCase.verifyTrue(~isempty(G.boundary_obs_s), 'Obstacle boundary nodes (P-grid) should not be empty');

                % Verify normal vectors
                testCase.verifyEqual(size(G.obs_normals_s, 1), size(G.boundary_obs_s, 1), ...
                                     'Normal vectors should match obstacle boundary nodes');
                testCase.verifyEqual(size(G.obs_normals_s, 2), 2, 'Normal vectors should be 2D');

                % Check that normals are unit vectors (within tolerance)
                normal_magnitudes = sqrt(sum(G.obs_normals_s.^2, 2));
                testCase.verifyTrue(all(abs(normal_magnitudes - 1) < 1e-10), ...
                                    'Normal vectors should be unit vectors');

                % Verify normal vectors are finite
                testCase.verifyTrue(all(isfinite(G.obs_normals_s(:))), ...
                                    'Normal vectors should be finite');

                % Verify obstacle IDs
                testCase.verifyEqual(length(G.boundary_obs_ids), size(G.boundary_obs_s, 1), ...
                                     'Obstacle IDs should match boundary nodes');
                testCase.verifyTrue(all(G.boundary_obs_ids >= 1 & G.boundary_obs_ids <= G.num_obstacles), ...
                                    'Obstacle IDs should be valid indices');

                % Verify proximity indices are logical
                testCase.verifyTrue(islogical(G.idx_near_obs_V), 'Near-obstacle V indices should be logical');
                testCase.verifyTrue(islogical(G.idx_near_obs_P), 'Near-obstacle P indices should be logical');
                testCase.verifyTrue(islogical(G.idx_far_boundaries_V), 'Far-boundary V indices should be logical');
                testCase.verifyTrue(islogical(G.idx_far_boundaries_P), 'Far-boundary P indices should be logical');

                % Check that some nodes are near obstacles (reasonable coverage)
                near_fraction_V = sum(G.idx_near_obs_V) / length(G.idx_near_obs_V);
                near_fraction_P = sum(G.idx_near_obs_P) / length(G.idx_near_obs_P);
                testCase.verifyTrue(near_fraction_V > 0.01, 'Some V-grid nodes should be near obstacles');
                testCase.verifyTrue(near_fraction_P > 0.01, 'Some P-grid nodes should be near obstacles');

                fprintf('[PASS] Multi-obstacle geometry structure test passed!\n');
                fprintf('  Obstacles: %d\n', G.num_obstacles);
                fprintf('  V-grid nodes: %d interior, %d boundary\n', size(G.xy, 1), ...
                        size(G.boundary_in, 1) + size(G.boundary_out, 1) + size(G.boundary_y, 1) + size(G.boundary_obs, 1));
                fprintf('  P-grid nodes: %d interior, %d boundary\n', size(G.xy_s, 1), ...
                        size(G.boundary_in_s, 1) + size(G.boundary_out_s, 1) + size(G.boundary_y_s, 1) + size(G.boundary_obs_s, 1));
                fprintf('  Obstacle boundary nodes: %d (V-grid), %d (P-grid)\n', ...
                        size(G.boundary_obs, 1), size(G.boundary_obs_s, 1));

            catch ME
                fprintf('[FAIL] Multi-obstacle geometry structure test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Multi-obstacle geometry structure test failed: %s', ME.message));
            end
        end

        function testBoundaryClassification(testCase)
            % Test that boundary nodes are correctly classified
            try
                setup_paths();

                % Load configuration and build geometry
                cfg = config(testCase.GEOMETRY_TYPE);
                G = build_geometry(cfg);

                % Check that boundary arrays are disjoint (no overlapping nodes)
                all_boundary_s = [G.boundary_in_s; G.boundary_out_s; G.boundary_y_s; G.boundary_obs_s];
                testCase.verifyEqual(size(unique(all_boundary_s, 'rows'), 1), size(all_boundary_s, 1), ...
                                     'Pressure boundary nodes should be disjoint');

                all_boundary = [G.boundary_in; G.boundary_out; G.boundary_y; G.boundary_obs];
                testCase.verifyEqual(size(unique(all_boundary, 'rows'), 1), size(all_boundary, 1), ...
                                     'Velocity boundary nodes should be disjoint');

                % Check that interior and boundary nodes are disjoint
                all_nodes_s = [G.xy_s; all_boundary_s];
                testCase.verifyEqual(size(unique(all_nodes_s, 'rows'), 1), size(all_nodes_s, 1), ...
                                     'Interior and boundary pressure nodes should be disjoint');

                all_nodes = [G.xy; all_boundary];
                testCase.verifyEqual(size(unique(all_nodes, 'rows'), 1), size(all_nodes, 1), ...
                                     'Interior and boundary velocity nodes should be disjoint');

                % Verify domain boundary nodes are at expected locations
                eps = cfg.mesh.boundary_eps;
                x_min = cfg.domain.x_min;
                x_max = cfg.domain.x_max;
                y_min = cfg.domain.y_min;
                y_max = cfg.domain.y_max;

                % Check inlet nodes (left boundary)
                if ~isempty(G.boundary_in_s)
                    testCase.verifyTrue(all(G.boundary_in_s(:, 1) < x_min + eps), ...
                                        'Inlet pressure nodes should be at left boundary');
                end
                if ~isempty(G.boundary_in)
                    testCase.verifyTrue(all(G.boundary_in(:, 1) < x_min + eps), ...
                                        'Inlet velocity nodes should be at left boundary');
                end

                % Check outlet nodes (right boundary)
                if ~isempty(G.boundary_out_s)
                    testCase.verifyTrue(all(G.boundary_out_s(:, 1) > x_max - eps), ...
                                        'Outlet pressure nodes should be at right boundary');
                end
                if ~isempty(G.boundary_out)
                    testCase.verifyTrue(all(G.boundary_out(:, 1) > x_max - eps), ...
                                        'Outlet velocity nodes should be at right boundary');
                end

                % Check wall nodes (top/bottom boundaries)
                if ~isempty(G.boundary_y_s)
                    testCase.verifyTrue(all(G.boundary_y_s(:, 2) > y_max - eps | G.boundary_y_s(:, 2) < y_min + eps), ...
                                        'Wall pressure nodes should be at top/bottom boundaries');
                end
                if ~isempty(G.boundary_y)
                    testCase.verifyTrue(all(G.boundary_y(:, 2) > y_max - eps | G.boundary_y(:, 2) < y_min + eps), ...
                                        'Wall velocity nodes should be at top/bottom boundaries');
                end

                fprintf('[PASS] Multi-obstacle boundary classification test passed!\n');

            catch ME
                fprintf('[FAIL] Multi-obstacle boundary classification test failed: %s\n', ME.message);
                testCase.verifyTrue(false, sprintf('Multi-obstacle boundary classification test failed: %s', ME.message));
            end
        end

    end

    methods

        function goldenFile = getGoldenFilePath(~)
            % Override to use multi-specific golden file path
            goldenFile = fullfile(fileparts(mfilename('fullpath')), 'golden', ...
                                  'multi_Re100_Nt20_dt0.01_seed42.mat');
        end

    end

end
