function [D_all] = RBF_PHS_FD_all(xy1,xy_s,Nearest_Idx,k,m,d)

% RBF_PHS_FD_ALL  RBF-FD-based (PHS+poly) global differentiation matrices
%   [D_ALL] = RBF_PHS_FD_ALL(XY1,XY_S,NEAREST_IDX,K,M,D,VARARGIN) returns
%   the global differentiation matrices from the node set {xy_s} to node set
%   {xy1}, including Dx, Dy, DL, Dxx, Dyy, Dxy.
%   NEAREST_IDX is the index matrix return by NEAREST_INTERP.
%   k is local RBF stencil size, m is the order of PHS-RBF, and d is the
%   degree of the polynomial Augmentation.
%   Reference:
%     [1] T. Chu, O. T. Schmidt, RBF-FD discretization of the Navier-Stokes 
%     equations on scattered but staggered nodes,
%     Journal of Computational Physics 474, 111756, 2023
%     [2] T. Chu, O. T. Schmidt, Mesh-free hydrodynamic stability,
%     Submitted to Journal of Computational Physics
%
%   T. Chu (tic173@ucsd.edu), O. T. Schmidt (oschmidt@ucsd.edu)
%%
X1            =           xy1(:,1);
NxNy1         =           length(X1);
X2            =           xy_s(:,1);
NxNy2         =           length(X2);
weight_x      =           zeros(NxNy1,k);
weight_y      =           zeros(NxNy1,k);
weight_L      =           zeros(NxNy1,k);
weight_xx     =           zeros(NxNy1,k);
weight_yy     =           zeros(NxNy1,k);
weight_xy     =           zeros(NxNy1,k);

 %%  Local RBF-FDs for each grid   
 
 
for m1 = 1:NxNy1  
    
    if mod(m1,100)==0
        disp(['j= ' num2str(m1)])
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    xk1       =   xy_s(Nearest_Idx(m1,1:k),:);     %  local stencil
    
    xe1       =   xy1(m1,:);                       %  evaluation point
    
    
%     [w] =  RBF_FD_PHS_pol_weights (xk1(:,1),xk1(:,2),xe1(:,1),xe1(:,2),m,d);                  % RBF-PHS method

    
    % move center to the origin
    
    x         =   xk1(:,1)-xe1(:,1);
    y         =   xk1(:,2)-xe1(:,2);

    
    % scale the local stencil
    
    scale     =   sqrt(x(end)^2+y(end)^2)/1;
    
    scale_x   =   scale; 
    scale_y   =   scale;
    
    x         =   x/scale_x;
    y         =   y/scale_y;
    
    
    
    % ------ RBF part --------------------------------------------------------
    
   
    A = hypot(bsxfun(@minus,x,x'),bsxfun(@minus,y,y')).^m; % RBF matrix
    
    L0 = m*(bsxfun(@times,(hypot(x,y)).^(m-2),[-x,-y,m*ones(k,1),1+(m-2)*x.^2./(hypot(x,y)).^(2),1+(m-2)*y.^2./(hypot(x,y)).^(2),(m-2)*x.*y./(hypot(x,y)).^(2)])); % RHSs
    
    
    for j = 1:k
        
        if x(j) == 0 && y(j) == 0
            L0(j,4) =0;
            L0(j,5) =0;
            L0(j,6) =0;
        end
        
    end
    
    
    % ------ Polynomial augmentation -------------------------------------------------
    
    if d == -1 % Special case with no polynomial terms,
        
        A_aug = A; L = L0; % i.e. pure RBF
        
    else % Create matrix with polynomial terms and matching constraints
        
        X = x(:,ones(1,d+1)); X(:,1) = 1; X = cumprod( X,2);
        Y = y(:,ones(1,d+1)); Y(:,1) = 1; Y = cumprod( Y,2);
        np = (d+1)*(d+2)/2; % Number of polynomial terms
        XY = zeros(k,np); col = 1; % assemble polynomial matrix block
        
        for j = 0:d
            XY(:,col:col+j) = X(:,j+1:-1:1).*Y(:,1:j+1);
            col = col+j+1;
        end
        
        L1 = zeros(np,6); % Create matching RHSs
        if d >= 1; L1(2,1) = 1; L1(3,2) = 1; end
        if d >= 2; L1(4,3) = 2; L1(4,4)=2; L1(6,3) = 2; L1(6,5)=2; L1(5,6)=1;end
        
        A_aug = [A,XY;XY',zeros(col-1)]; % Augmented interpolation matrix
        
        L = [L0;L1]; % assemble RHSs
        
    end
    % ------ Solve for weights -----------------------------------------------
    
    A_aug = (A_aug+A_aug')/2;

      W = A_aug\L;
    
% [U,S,V]   =   svd(A_aug);
% 
% S  =  diag(S);
% 
% id =  (S<1e-10);
% 
% U(:,id)=[]; V(:,id)=[]; S(id)=[];
% S      =    diag(1./S);
% A_inv =    V*S*(U');
% 
% W = A_inv*L;
    w = W(1:k,:);      % Extragct the RBF-FD weights

    weight_x(m1,1:k)   =  w(:,1)/scale_x;

    weight_y(m1,1:k)   =  w(:,2)/scale_y;

    weight_L(m1,1:k)   =  w(:,3)/(scale^2);

    weight_xx(m1,1:k)  =  w(:,4)/scale_x^2;

    weight_yy(m1,1:k)  =  w(:,5)/scale_y^2;

    weight_xy(m1,1:k)  =  w(:,6)/scale_x/scale_y;

        
            
end
%%  Assigning RBF-FD weights to global matrices
Idx_x                                 =       reshape((1:NxNy1).*ones(k,1),1,[]);

Idx_y                                 =       reshape(Nearest_Idx(:,1:k)',1,[]);
Wx                                    =       reshape(weight_x',1,[]);
Wy                                    =       reshape(weight_y',1,[]);
WL                                    =       reshape(weight_L',1,[]);
Wxx                                   =       reshape(weight_xx',1,[]);
Wyy                                   =       reshape(weight_yy',1,[]);
Wxy                                   =       reshape(weight_xy',1,[]);
Dx                                    =       sparse(Idx_x,Idx_y,Wx,NxNy1,NxNy2);
Dy                                    =       sparse(Idx_x,Idx_y,Wy,NxNy1,NxNy2);
L                                     =       sparse(Idx_x,Idx_y,WL,NxNy1,NxNy2);
Dxx                                   =       sparse(Idx_x,Idx_y,Wxx,NxNy1,NxNy2);
Dyy                                   =       sparse(Idx_x,Idx_y,Wyy,NxNy1,NxNy2);
Dxy                                   =       sparse(Idx_x,Idx_y,Wxy,NxNy1,NxNy2);

D_all                                 =       {Dx,Dy,L,Dxx,Dyy,Dxy};
end