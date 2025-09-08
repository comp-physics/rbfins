classdef GoldenFileGenerator
% GOLDENFILEGENERATOR Utility class for generating golden reference files
%
% This class provides static methods to generate golden files for any geometry
% without code duplication.

methods (Static)
    function generateGoldenFile(geometryType)
    % Generate a golden file for the specified geometry
    %
    % Parameters:
    %   geometryType   - String: 'cylinder', 'ellipse', or 'rectangle'

    fprintf('=== Creating %s Golden File ===\n', upper(geometryType));

    % Set up environment
    setenv('CI', 'true');
    setenv('MATLAB_TEST', 'true');
    setup_paths();

    % Store parameters before simulation (simulate may clear variables)
    geomType = geometryType;

    % Load config for metadata first
    % Store geometry type before simulate clears variables
    geometryType = geomType;

    cfg = config(geometryType);

    % Run the simulation
    fprintf('Running %s simulation...\n', geometryType);
    setenv('CONFIG_FUNC', 'config');
    setenv('GEOMETRY_TYPE', geometryType);
    simulate;

    % Restore variables that may have been cleared by simulate
    geometryType = getenv('GEOMETRY_TYPE');
    cfg = config(geometryType);

    % Extract the final state
    % Ensure variables exist from simulate script
    if ~exist('xy1', 'var') || ~exist('W', 'var')
        error('Required variables xy1 and W not found after simulation');
    end
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
    meta.algorithm = sprintf('RBF-FD Navier-Stokes %s flow', geometryType);
    meta.geometry_type = cfg.geometry.type;
    meta.reynolds_number = cfg.simulation.reynolds_number;
    meta.time_step = cfg.simulation.time_step;
    meta.num_time_steps = cfg.simulation.num_time_steps_ci;
    meta.random_seed = cfg.simulation.random_seed;

    % Add geometry-specific metadata
    switch geometryType
    case 'cylinder'
        meta.obstacle_radius = cfg.geometry.obstacle_radius;
    case 'ellipse'
        meta.ellipse_a = cfg.geometry.ellipse_a;
        meta.ellipse_b = cfg.geometry.ellipse_b;
    case 'rectangle'
        meta.rect_width = cfg.geometry.rect_width;
        meta.rect_height = cfg.geometry.rect_height;
        meta.rect_x_center = cfg.geometry.rect_x_center;
        meta.rect_y_center = cfg.geometry.rect_y_center;
    end

    meta.domain_size = [cfg.domain.x_min, cfg.domain.x_max, cfg.domain.y_min, cfg.domain.y_max];
    meta.timestamp = char(datetime('now'));
    meta.matlab_version = version;

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
% generateAllGoldenFiles Generates golden files for all supported geometries
fprintf('=== GENERATING ALL GOLDEN FILES ===\n');
setup_paths(); % Ensure paths are set for all calls

geometries = {'cylinder', 'ellipse', 'rectangle'};

for i = 1:length(geometries)
    try
        GoldenFileGenerator.generateGoldenFile(geometries{i});
    catch ME
        fprintf('❌ Failed to generate %s golden file: %s\n', geometries{i}, ME.message);
    end
    % Clear variables for next iteration
    clear xy1 W U V xy1_sorted U_sorted V_sorted meta gold;
end
fprintf('=== ALL GOLDEN FILES GENERATION COMPLETE ===\n');
end
end
end
