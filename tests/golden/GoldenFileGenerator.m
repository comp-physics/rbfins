classdef GoldenFileGenerator
    % GOLDENFILEGENERATOR Utility class for generating golden reference files
    %
    % This class provides static methods to generate golden files for any geometry
    % without code duplication.

    methods (Static)
        function generateGoldenFile(geometryType, configFunction)
            % Generate a golden file for the specified geometry
            %
            % Parameters:
            %   geometryType   - String: 'cylinder', 'ellipse', or 'rectangle'
            %   configFunction - String: 'config_cylinder', 'config_ellipse', or 'config_rectangle'

            fprintf('=== Creating %s Golden File ===\n', upper(geometryType));

            % Set up environment
            setenv('CI', 'true');
            setenv('MATLAB_TEST', 'true');
            setup_paths();

            % Run the simulation
            fprintf('Running %s simulation...\n', geometryType);
            setenv('CONFIG_FUNC', configFunction);
            simulate;

            % Load config for metadata
            cfg = feval(configFunction);

            % Extract the final state
            n = size(xy1, 1);
            U = W(1:n, end);
            V = W(n+1:end, end);

            fprintf('Extracting simulation results...\n');
            fprintf('  Number of nodes: %d\n', n);
            fprintf('  Max |U|: %.6f\n', max(abs(U)));
            fprintf('  Max |V|: %.6f\n', max(abs(V)));

            % Sort nodes for reproducible ordering
            [~, idx] = sortrows(xy1, [1, 2]);
            xy1_sorted = xy1(idx, :);
            U_sorted = U(idx);
            V_sorted = V(idx);

            % Create metadata
            meta = GoldenFileGenerator.createMetadata(geometryType, cfg);

            % Create golden data structure
            gold.xy1 = xy1_sorted;
            gold.U = U_sorted;
            gold.V = V_sorted;
            gold.meta = meta;

            % Save golden file
            outDir = fullfile('tests', 'golden');
            if ~exist(outDir, 'dir')
                mkdir(outDir);
            end

            outFile = fullfile(outDir, sprintf('%s_Re%d_Nt%d_dt%g_seed%d.mat', ...
                geometryType, meta.reynolds_number, meta.num_time_steps, meta.time_step, meta.random_seed));

            save(outFile, 'gold', '-v7');

            fprintf('=== %s Golden File Created Successfully ===\n', upper(geometryType));
            fprintf('  File: %s\n', outFile);
            fprintf('  Nodes: %d\n', length(gold.U));
            fprintf('  Max |U|: %.6f\n', max(abs(gold.U)));
            fprintf('  Max |V|: %.6f\n', max(abs(gold.V)));
            fprintf('  Mean U: %.6f\n', mean(gold.U));
            fprintf('  Mean V: %.6f\n', mean(gold.V));
        end

        function generateAllGoldenFiles()
            % Generate golden files for all supported geometries

            fprintf('=== GENERATING ALL GOLDEN FILES ===\n');

            % Define geometries and their config functions
            geometries = {
                {'cylinder', 'config_cylinder'};
                {'ellipse', 'config_ellipse'};
                {'rectangle', 'config_rectangle'}
            };

            for i = 1:length(geometries)
                geomType = geometries{i}{1};
                configFunc = geometries{i}{2};

                try
                    fprintf('\n=== GENERATING %s GOLDEN FILE ===\n', upper(geomType));
                    GoldenFileGenerator.generateGoldenFile(geomType, configFunc);

                    % Clear variables for next iteration
                    evalin('base', 'clear xy1 W U V xy1_sorted U_sorted V_sorted meta gold;');

                catch ME
                    fprintf('❌ Failed to generate %s golden file: %s\n', geomType, ME.message);
                end
            end

            fprintf('\n=== ALL GOLDEN FILE GENERATION COMPLETE ===\n');
        end

        function meta = createMetadata(geometryType, cfg)
            % Create metadata structure for golden file

            meta.algorithm = sprintf('RBF-FD Navier-Stokes %s flow', geometryType);
            meta.geometry_type = cfg.geometry.type;
            meta.reynolds_number = cfg.simulation.reynolds_number;
            meta.time_step = cfg.simulation.time_step;
            meta.num_time_steps = cfg.simulation.num_time_steps_ci;
            meta.random_seed = cfg.simulation.random_seed;
            meta.domain_size = [cfg.domain.x_min, cfg.domain.x_max, cfg.domain.y_min, cfg.domain.y_max];
            meta.timestamp = datestr(now);
            meta.matlab_version = version;

            % Add geometry-specific metadata
            switch geometryType
                case 'cylinder'
                    meta.cylinder_radius = cfg.geometry.obstacle_radius;
                case 'ellipse'
                    meta.ellipse_a = cfg.geometry.ellipse_a;
                    meta.ellipse_b = cfg.geometry.ellipse_b;
                case 'rectangle'
                    meta.rect_width = cfg.geometry.rect_width;
                    meta.rect_height = cfg.geometry.rect_height;
                otherwise
                    warning('Unknown geometry type: %s', geometryType);
            end
        end
    end
end
