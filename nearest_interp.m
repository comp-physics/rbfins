function [Nearest_Idx] = nearest_interp(xy,xy_s,k,varargin)

% NEAREST_INTERP  for each node in xy, find its k nearest nerighbors in xy_s

% [NEAREST_IDX] = NEAREST_INTERP(XY,XY_S,K,VARARGIN) returns the indices of
% the nearest k neighbors in the node set {xy_s} for each node in {xy}.
% R0 is the maximum radius for local searching.
% If not specified, R0 = 1e-5.


%   Reference:
%     [1] T. Chu, O. T. Schmidt, RBF-FD discretization of the Navier-Stokes 
%     equations on scattered but staggered nodes,
%     Journal of Computational Physics 474, 111756, 2023
%     [2] T. Chu, O. T. Schmidt, Mesh-free hydrodynamic stability,
%     Submitted to Journal of Computational Physics
%
%  T. Chu (tic173@ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)



X1            =           xy(:,1);

NxNy1         =           length(X1);


if nargin==4
    R0 = varargin{1};
else
    R0 = 1e-5;
end

%%  Find nearest nodes
 
Nearest_Idx     =       zeros(NxNy1,k);



disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Find nearest ' num2str(k) ' nodes'])




for j=1:NxNy1
    
    if mod(j,100)==0
        disp(['j= ' num2str(j)])
    end
    
    Node_j                       =       xy(j,:);
    
    
    r2_j                         =       (xy_s-Node_j).^2;
    r2_j                         =       r2_j(:,1)+r2_j(:,2);
    
    
    
    % Neighbor region
    
    Idx_n                        =       find(r2_j<R0);
    
    if length(Idx_n)<k
        
        
        Idx_n                        =       find(r2_j<10^2);
        
    end
    
    
    % search k nearest nodes in the neighbor region
    
    
    Idx_0                        =       knnsearch(Node_j,xy_s(Idx_n,:),k);    
    
    Idx                          =       Idx_n(Idx_0(1:k));  
    
    Nearest_Idx(j,:)             =       Idx;    
     
    
    
    
    
end










end
