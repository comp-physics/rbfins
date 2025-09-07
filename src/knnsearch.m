function [idx, D] = knnsearch(varargin)
%KNNSEARCH Linear k-nearest neighbor (KNN) search for RBF-FD applications
%
% This function performs efficient k-nearest neighbor search using linear
% distance computation. While typically considered suitable only for small
% datasets, the vectorized MATLAB implementation with JIT acceleration
% makes it competitive with tree-based methods for moderate-sized problems
% common in RBF-FD applications.
%
% SYNTAX:
%   IDX = knnsearch(Q,R,K) - Find K nearest neighbors in R for each point in Q
%   IDX = knnsearch(Q,R)   - Find 1 nearest neighbor (K=1 default)
%   IDX = knnsearch(Q)     - Self-search: find neighbors in Q for each point in Q
%   [IDX,D] = knnsearch(...) - Also return distances
%
% INPUTS:
%   Q - Query points [m x d] where m = number of queries, d = dimension
%   R - Reference points [n x d] where n = number of reference points
%   K - Number of nearest neighbors to find (default: K=1)
%
% OUTPUTS:
%   IDX - Indices of nearest neighbors [m x K]
%   D   - Distances to nearest neighbors [m x K] (if requested)
%
% ALGORITHM:
%   Uses brute-force distance computation: d(i,j) = ||Q(i,:) - R(j,:)||
%   For each query point, computes distances to all reference points,
%   sorts, and returns K smallest. Self-distance is excluded when Q=R.
%
% NOTE:
%   Optimized for RBF-FD stencil construction where high accuracy is needed
%   and dataset sizes are moderate (typically < 10^5 points).
%
% AUTHOR: Yi Cao, Cranfield University, 25 March 2008
% REFERENCE: Used in RBF-FD method for mesh-free discretization

%% Input parsing and initialization
[Q, R, K, fident] = parseinputs(varargin{:}); % Parse and validate inputs
error(nargoutchk(0, 2, nargout)); % Validate number of outputs

% Problem dimensions
[N, M] = size(Q); % N = number of query points, M = dimension
L = size(R, 1); % L = number of reference points

% Preallocate output arrays
idx = zeros(N, K); % Indices of nearest neighbors
D = idx; % Distances to nearest neighbors (same size)
%% Main search algorithm
if K == 1
    %% Single nearest neighbor search (optimized)
    for k = 1:N
        d = zeros(L, 1); % Initialize squared distance array

        % Compute squared Euclidean distances to all reference points
        for t = 1:M
            d = d + (R(:, t) - Q(k, t)).^2; % Accumulate (x_t - q_t)^2
        end

        % Exclude self-distance if query and reference sets are identical
        if fident
            d(k) = inf; % Make self-distance infinite to exclude it
        end

        % Find single nearest neighbor
        [D(k), idx(k)] = min(d); % Minimum distance and its index
    end
else
    %% Multiple nearest neighbors search (requires sorting)
    for k = 1:N
        d = zeros(L, 1); % Initialize squared distance array

        % Compute squared Euclidean distances to all reference points
        for t = 1:M
            d = d + (R(:, t) - Q(k, t)).^2; % Accumulate (x_t - q_t)^2
        end

        % Exclude self-distance if query and reference sets are identical
        if fident
            d(k) = inf; % Make self-distance infinite to exclude it
        end

        % Sort distances and extract K nearest neighbors
        [s, t] = sort(d); % Sort distances (s) and get indices (t)
        idx(k, :) = t(1:K); % Store indices of K nearest neighbors
        D(k, :) = s(1:K); % Store corresponding distances
    end
end
%% Convert squared distances to actual distances if requested
if nargout > 1
    D = sqrt(D); % Take square root to get Euclidean distances
end

function [Q, R, K, fident] = parseinputs(varargin)
%PARSEINPUTS Parse and validate input arguments for knnsearch
%
% Handles the various input syntax options and validates parameters
% Validate number of input arguments
error(nargchk(1, 3, nargin));

% Parse query points (required)
Q = varargin{1};

% Parse reference points (optional, defaults to Q for self-search)
if nargin < 2
    R = Q; % Self-search: find neighbors within Q
    fident = true; % Flag indicating identical query and reference sets
else
    fident = false; % Different query and reference sets
    R = varargin{2};
end

% Handle empty reference set (equivalent to self-search)
if isempty(R)
    fident = true;
    R = Q;
end

% Check if query and reference sets are actually identical
if ~fident
    fident = isequal(Q, R);
end

% Parse number of neighbors K (optional, defaults to 1)
if nargin < 3
    K = 1; % Single nearest neighbor
else
    K = varargin{3}; % User-specified number of neighbors
end
