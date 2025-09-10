function [Nearest_Idx] = nearest_interp(xy, xy_s, k, varargin)
  %NEAREST_INTERP Find k nearest neighbors for RBF-FD stencil construction
  %
  % This function finds the k nearest neighbors in the source node set {xy_s}
  % for each target node in {xy}. This is essential for building local RBF-FD
  % stencils that provide good approximations to differential operators.
  %
  % INPUTS:
  %   xy      - Target node coordinates [N1 x 2] where stencils are needed
  %   xy_s    - Source node coordinates [N2 x 2] where function values are known
  %   k       - Number of nearest neighbors to find for each target node
  %   R0      - (Optional) Maximum search radius for initial neighbor filtering
  %             Default: R0 = 1e-5
  %
  % OUTPUTS:
  %   Nearest_Idx - Matrix [N1 x k] containing indices of k nearest neighbors
  %                 for each target node. Nearest_Idx(i,:) contains the indices
  %                 in xy_s of the k closest nodes to xy(i,:)
  %
  % ALGORITHM:
  %   1. For each target node, compute squared distances to all source nodes
  %   2. Filter source nodes within search radius R0 (or expand if insufficient)
  %   3. Use knnsearch to find k nearest neighbors within filtered set
  %   4. Return global indices in xy_s
  %
  % REFERENCE:
  %   [1] T. Chu, O. T. Schmidt, "RBF-FD discretization of the Navier-Stokes
  %       equations on scattered but staggered nodes", Journal of Computational
  %       Physics 474, 111756, 2023
  %   [2] T. Chu, O. T. Schmidt, "Mesh-free hydrodynamic stability",
  %       Submitted to Journal of Computational Physics
  %
  % AUTHORS: T. Chu (tic173@ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
  %% Initialize variables and parameters
  X1 = xy(:, 1); % Extract x-coordinates (not used but kept for consistency)
  NxNy1 = length(X1); % Number of target nodes

  % Parse optional search radius parameter
  if nargin == 4
    R0 = varargin{1}; % User-specified search radius
  else
    R0 = 1e-5; % Default search radius (very small for coincident nodes)
  end
  %% Initialize output array and begin neighbor search
  Nearest_Idx = zeros(NxNy1, k); % Preallocate output matrix

  % Load configuration for progress display control
  try
    cfg = config();
    show_progress = cfg.simulation.show_progress;
  catch
    show_progress = true; % Default to showing progress if config unavailable
  end

  % Check if running in CI/test environment
  isCI = strcmpi(getenv('CI'), 'true');
  isTest = strcmpi(getenv('MATLAB_TEST'), 'true');

  if show_progress && ~isCI && ~isTest
    disp('Finding nearest neighbors for RBF-FD stencils...');
    disp(['Target: ', num2str(k), ' neighbors per node for ', num2str(NxNy1), ' nodes']);
  end
  %% Main loop: Find k nearest neighbors for each target node
  for j = 1:NxNy1
    % Progress indicator for large problems
    if mod(j, 100) == 0 && show_progress && ~isCI && ~isTest
      disp(['Processing node ', num2str(j), ' of ', num2str(NxNy1)]);
    end

    % Current target node coordinates
    Node_j = xy(j, :);

    % Compute squared distances from target node to all source nodes
    r2_j = (xy_s - Node_j).^2; % Element-wise squared differences
    r2_j = r2_j(:, 1) + r2_j(:, 2); % Sum to get squared Euclidean distances
    %% Filter source nodes within search radius
    % Find candidate neighbors within initial search radius
    Idx_n = find(r2_j < R0);

    % If insufficient neighbors found, expand search radius significantly
    if length(Idx_n) < k
      Idx_n = find(r2_j < 10^2); % Use large radius to include all nodes if needed
    end
    %% Find k nearest neighbors within candidate set
    % Use efficient knnsearch on reduced candidate set
    Idx_0 = knnsearch(Node_j, xy_s(Idx_n, :), k);

    % Convert local indices back to global indices in xy_s
    Idx = Idx_n(Idx_0(1:k));

    % Store result for current target node
    Nearest_Idx(j, :) = Idx;
  end

  if show_progress && ~isCI && ~isTest
    disp('Nearest neighbor search completed.');
  end

end
