% Simple golden file generator that extracts data from an existing simulation run
function make_golden_simple()

fprintf('=== Creating Simple Golden File ===\n');

% Set up environment
setenv('CI', 'true');
setenv('MATLAB_TEST', 'true');
setup_paths();

% Run the simulation
fprintf('Running simulation...\n');
simulate;

% Load config for metadata
cfg = config();

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
meta.algorithm = 'RBF-FD Navier-Stokes cylinder flow';
meta.reynolds_number = cfg.simulation.reynolds_number;
meta.time_step = cfg.simulation.time_step;
meta.num_time_steps = cfg.simulation.num_time_steps_ci;
meta.random_seed = cfg.simulation.random_seed;
meta.cylinder_radius = cfg.mesh.cylinder_radius;
meta.domain_size = [cfg.domain.x_min, cfg.domain.x_max, cfg.domain.y_min, cfg.domain.y_max];
meta.timestamp = datestr(now);
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

outFile = fullfile(outDir, sprintf('cylinder_Re%d_Nt%d_dt%g_seed%d.mat', ...
    meta.reynolds_number, meta.num_time_steps, meta.time_step, meta.random_seed));

save(outFile, 'gold', '-v7');

fprintf('=== Golden File Created Successfully ===\n');
fprintf('  File: %s\n', outFile);
fprintf('  Nodes: %d\n', length(gold.U));
fprintf('  Max |U|: %.6f\n', max(abs(gold.U)));
fprintf('  Max |V|: %.6f\n', max(abs(gold.V)));
fprintf('  Mean U: %.6f\n', mean(gold.U));
fprintf('  Mean V: %.6f\n', mean(gold.V));

end
