function [D_all] = RBF_PHS_FD_all(xy1, xy_s, Nearest_Idx, k, m, d)
%RBF_PHS_FD_ALL Generate global RBF-FD differentiation matrices using PHS + polynomials
%
% This function creates global differentiation matrices using Radial Basis Function
% Finite Difference (RBF-FD) method with Polyharmonic Spline (PHS) basis functions
% augmented with polynomial terms for enhanced accuracy.
%
% INPUTS:
%   xy1         - Target node coordinates [N1 x 2] where derivatives are evaluated
%   xy_s        - Source node coordinates [N2 x 2] where function values are known
%   Nearest_Idx - Index matrix [N1 x k] of k nearest neighbors for each target node
%   k           - Stencil size (number of nearest neighbors used in local stencils)
%   m           - Order of PHS-RBF (r^m, typically m=3 for 2D problems)
%   d           - Degree of polynomial augmentation (-1=none, 0=constant, 1=linear, 2=quadratic)
%
% OUTPUTS:
%   D_all       - Cell array containing sparse differentiation matrices:
%                 {Dx, Dy, L, Dxx, Dyy, Dxy} where:
%                 Dx, Dy  = first derivatives (d/dx, d/dy)
%                 L       = Laplacian (d^2/dx^2 + d^2/dy^2)
%                 Dxx,Dyy = second derivatives (d^2/dx^2, d^2/dy^2)
%                 Dxy     = mixed derivative (d^2/dxdy)
%
% METHOD:
%   For each target node, a local stencil of k source nodes is used to construct
%   RBF-FD weights. The method combines:
%   1. PHS-RBF interpolation: phi(r) = r^m (polyharmonic splines)
%   2. Polynomial augmentation: adds polynomial terms up to degree d
%   3. Local coordinate scaling for numerical stability
%
% REFERENCE:
%   [1] T. Chu, O. T. Schmidt, "RBF-FD discretization of the Navier-Stokes
%       equations on scattered but staggered nodes", Journal of Computational
%       Physics 474, 111756, 2023
%   [2] T. Chu, O. T. Schmidt, "Mesh-free hydrodynamic stability",
%       Submitted to Journal of Computational Physics
%
% AUTHORS: T. Chu (tic173@ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)

%% Initialize variables and storage arrays
X1 = xy1(:, 1); % Extract x-coordinates of target nodes
NxNy1 = length(X1); % Number of target nodes (where derivatives are evaluated)
X2 = xy_s(:, 1); % Extract x-coordinates of source nodes
NxNy2 = length(X2); % Number of source nodes (where function values are known)

% Preallocate weight matrices for each derivative operator
weight_x = zeros(NxNy1, k); % Weights for d/dx operator
weight_y = zeros(NxNy1, k); % Weights for d/dy operator
weight_L = zeros(NxNy1, k); % Weights for Laplacian operator
weight_xx = zeros(NxNy1, k); % Weights for d^2/dx^2 operator
weight_yy = zeros(NxNy1, k); % Weights for d^2/dy^2 operator
weight_xy = zeros(NxNy1, k); % Weights for d^2/dxdy operator
%% Main loop: Compute RBF-FD weights for each target node
% Process each target node to build local RBF-FD stencils

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

for m1 = 1:NxNy1
    % Progress indicator for large problems
    if mod(m1, 100) == 0 && show_progress && ~isCI && ~isTest
        disp(['Computing RBF-FD weights: node ', num2str(m1), ' of ', num2str(NxNy1)])
    end
    %% Extract local stencil and evaluation point
    xk1 = xy_s(Nearest_Idx(m1, 1:k), :); % Local stencil: k nearest source nodes
    xe1 = xy1(m1, :); % Current target node (evaluation point)
    %% Coordinate transformation: translate and scale for numerical stability

    % Translate stencil so evaluation point is at origin
    x = xk1(:, 1) - xe1(:, 1); % Relative x-coordinates
    y = xk1(:, 2) - xe1(:, 2); % Relative y-coordinates

    % Scale coordinates using distance to farthest stencil node
    % This improves numerical conditioning of the RBF interpolation matrix
    scale = sqrt(x(end)^2+y(end)^2) / 1; % Scaling factor based on stencil size

    scale_x = scale; % Same scaling for both x and y directions
    scale_y = scale;

    x = x / scale_x; % Scaled x-coordinates
    y = y / scale_y; % Scaled y-coordinates
    %% Build RBF interpolation matrix and derivative operators

    % Construct PHS-RBF interpolation matrix: A(i,j) = ||x_i - x_j||^m
    A = hypot(bsxfun(@minus, x, x'), bsxfun(@minus, y, y')).^m;

    % Build right-hand side matrix for derivative operators
    % Each column corresponds to a different derivative: [dx, dy, Laplacian, dxx, dyy, dxy]
    % Using analytical derivatives of PHS-RBF: d/dx(r^m) = m*r^(m-2)*x, etc.
    r = hypot(x, y); % Distance from origin for each stencil point
    r_m2 = r.^(m - 2); % r^(m-2) appears in all derivative formulas

    L0 = m * bsxfun(@times, r_m2, ...
        [-x, ... % d/dx:     -m*r^(m-2)*x
        -y, ... % d/dy:     -m*r^(m-2)*y
        m * ones(k, 1), ... % Laplacian: m*r^(m-2)*m = m^2*r^(m-2)
        1 + (m - 2) * x.^2 ./ r.^2, ... % d^2/dx^2: m*r^(m-2)*(1+(m-2)*x^2/r^2)
        1 + (m - 2) * y.^2 ./ r.^2, ... % d^2/dy^2: m*r^(m-2)*(1+(m-2)*y^2/r^2)
        (m - 2) * x .* y ./ r.^2]); % d^2/dxdy: m*r^(m-2)*(m-2)*x*y/r^2

    % Handle special case at origin (avoid 0/0)
    for j = 1:k
        if x(j) == 0 && y(j) == 0
            L0(j, 4) = 0; % d^2/dx^2 = 0 at origin
            L0(j, 5) = 0; % d^2/dy^2 = 0 at origin
            L0(j, 6) = 0; % d^2/dxdy = 0 at origin
        end
    end
    %% Polynomial augmentation for enhanced accuracy

    if d == -1 % Pure RBF case (no polynomial augmentation)
        A_aug = A; % Use RBF matrix as-is
        L = L0; % Use RBF derivative operators as-is

    else % Augment RBF with polynomial terms up to degree d

        % Build polynomial basis: 1, x, y, x^2, xy, y^2, x^3, x^2*y, xy^2, y^3, ...
        X = x(:, ones(1, d+1));
        X(:, 1) = 1;
        X = cumprod(X, 2); % Powers of x: [1, x, x^2, ...]
        Y = y(:, ones(1, d+1));
        Y(:, 1) = 1;
        Y = cumprod(Y, 2); % Powers of y: [1, y, y^2, ...]

        np = (d + 1) * (d + 2) / 2; % Total number of polynomial terms for degree d
        XY = zeros(k, np); % Polynomial matrix block
        col = 1;

        % Assemble polynomial terms: x^i * y^j for i+j <= d
        for j = 0:d
            XY(:, col:col+j) = X(:, j+1:-1:1) .* Y(:, 1:j+1); % [x^j, x^(j-1)*y, ..., y^j]
            col = col + j + 1;
        end

        % Create polynomial matching conditions (derivatives of polynomial terms)
        L1 = zeros(np, 6); % RHS for polynomial constraints
        if d >= 1 % Linear terms: derivatives of x and y
            L1(2, 1) = 1; % d/dx(x) = 1
            L1(3, 2) = 1; % d/dy(y) = 1
        end
        if d >= 2 % Quadratic terms: derivatives of x^2, xy, y^2
            L1(4, 3) = 2;
            L1(4, 4) = 2; % Laplacian of x^2: d^2/dx^2(x^2) = 2, d^2/dy^2(x^2) = 0 -> sum = 2
            L1(6, 3) = 2;
            L1(6, 5) = 2; % Laplacian of y^2: d^2/dx^2(y^2) = 0, d^2/dy^2(y^2) = 2 -> sum = 2
            L1(5, 6) = 1; % d^2/dxdy(xy) = 1
        end

        % Assemble augmented system: [RBF + polynomials; polynomial constraints]
        A_aug = [A, XY; XY', zeros(col-1)]; % Augmented interpolation matrix
        L = [L0; L1]; % Combined RHS (RBF + polynomial derivatives)

    end

    %% Solve linear system for RBF-FD weights

    % Ensure matrix is symmetric (numerical symmetry)
    A_aug = (A_aug + A_aug') / 2;

    % Solve augmented system: A_aug * weights = L
    W = A_aug \ L; % Direct solve using MATLAB's backslash operator

    % Alternative SVD-based solve for ill-conditioned cases (commented out)
    % [U,S,V] = svd(A_aug);
    % S = diag(S);
    % id = (S<1e-10);  % Remove small singular values
    % U(:,id)=[]; V(:,id)=[]; S(id)=[];
    % S = diag(1./S);
    % A_inv = V*S*(U');
    % W = A_inv*L;

    % Extract RBF weights (ignore polynomial constraint multipliers)
    w = W(1:k, :); % Only first k rows correspond to RBF weights
    %% Scale weights back to original coordinate system
    % Apply inverse coordinate scaling to get weights in physical coordinates
    weight_x(m1, 1:k) = w(:, 1) / scale_x; % d/dx weights
    weight_y(m1, 1:k) = w(:, 2) / scale_y; % d/dy weights
    weight_L(m1, 1:k) = w(:, 3) / (scale^2); % Laplacian weights
    weight_xx(m1, 1:k) = w(:, 4) / scale_x^2; % d^2/dx^2 weights
    weight_yy(m1, 1:k) = w(:, 5) / scale_y^2; % d^2/dy^2 weights
    weight_xy(m1, 1:k) = w(:, 6) / scale_x / scale_y; % d^2/dxdy weights


end
%% Assemble global sparse differentiation matrices
% Convert local RBF-FD weights into global sparse matrices

% Create index arrays for sparse matrix construction
Idx_x = reshape((1:NxNy1).*ones(k, 1), 1, []); % Row indices (repeated for each stencil)
Idx_y = reshape(Nearest_Idx(:, 1:k)', 1, []); % Column indices (stencil neighbors)

% Flatten weight matrices into vectors
Wx = reshape(weight_x', 1, []); % d/dx weights
Wy = reshape(weight_y', 1, []); % d/dy weights
WL = reshape(weight_L', 1, []); % Laplacian weights
Wxx = reshape(weight_xx', 1, []); % d^2/dx^2 weights
Wyy = reshape(weight_yy', 1, []); % d^2/dy^2 weights
Wxy = reshape(weight_xy', 1, []); % d^2/dxdy weights

% Construct sparse global differentiation matrices
Dx = sparse(Idx_x, Idx_y, Wx, NxNy1, NxNy2); % First derivative in x
Dy = sparse(Idx_x, Idx_y, Wy, NxNy1, NxNy2); % First derivative in y
L = sparse(Idx_x, Idx_y, WL, NxNy1, NxNy2); % Laplacian operator
Dxx = sparse(Idx_x, Idx_y, Wxx, NxNy1, NxNy2); % Second derivative in x
Dyy = sparse(Idx_x, Idx_y, Wyy, NxNy1, NxNy2); % Second derivative in y
Dxy = sparse(Idx_x, Idx_y, Wxy, NxNy1, NxNy2); % Mixed derivative

% Package all operators into output cell array
D_all = {Dx, Dy, L, Dxx, Dyy, Dxy};

end
