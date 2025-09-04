clc;
clear;

% Load Configuration
cfg = config();

% Check if running in CI environment or headless mode
isCI = strcmpi(getenv('CI'), 'true');
isTest = strcmpi(getenv('MATLAB_TEST'), 'true');

% Completely disable plotting for CI/test environments
doPlot = ~isCI && ~isTest;

if ~doPlot
    % Aggressively disable all graphics for tests
    set(0, 'DefaultFigureVisible', 'off');
    fprintf('Plotting completely disabled for testing/CI\n');
    hasDisplay = false;
else
    % Only check display if we might want to plot
    try
        % Try to create a figure to test if display is available
        f = figure('Visible', 'off');
        close(f);
        hasDisplay = true;
    catch
        hasDisplay = false;
        doPlot = false;
    end
end

% Set random seed for DistMesh reproducibility (DistMesh uses rand() for rejection method)
if isfield(cfg.simulation, 'random_seed') && ~isempty(cfg.simulation.random_seed)
    rng(cfg.simulation.random_seed);
    fprintf('Random seed set to %d for DistMesh reproducibility\n', cfg.simulation.random_seed);
end

% Ensure DistMesh is available
if ~exist('distmesh2d', 'file')
    % Try to add distmesh from the standard location
    if exist('distmesh', 'dir')
        addpath('distmesh');
    else
        error('distmesh directory not found. Please run setup_paths() or add distmesh to your path manually.');
    end
end

% Extract domain parameters
x_min = cfg.domain.x_min;
x_max = cfg.domain.x_max;
y_min = cfg.domain.y_min;
y_max = cfg.domain.y_max;

% Extract mesh parameters
dist = cfg.mesh.dist;
radius = cfg.mesh.cylinder_radius;
eps = cfg.mesh.boundary_eps;
a1 = cfg.mesh.refine_a1;
b1 = cfg.mesh.refine_b1;
a2 = cfg.mesh.refine_a2;
b2 = cfg.mesh.refine_b2;

%  local grid refinement
fd    = @(p) ddiff(drectangle(p,x_min, x_max,y_min,y_max),dcircle(p,0,0,radius));
fd1   = @(p) min(a1+b1*abs(dcircle(p,0,0,radius)),a2+b2*abs(dpoly(p,[radius 0; x_max 0])));
 
fix    = [x_min,y_min;x_min,y_max; x_max, y_max; x_max, y_min];

% Generate nodes using DistMesh
[xy_s,xt] = distmesh2d(fd,fd1,dist, [x_min,y_min; x_max, y_max],fix );

nodes ={xy_s,xt};

xy_s = nodes{1}; xt = nodes{2};

xy = zeros(length(xt)*cfg.mesh.edge_multiplier,2);
for j = 1:length(xt)
    xy((j-1)*cfg.mesh.edge_multiplier+1,:) = (xy_s(xt(j,1),:)+xy_s(xt(j,2),:))/2;
    xy((j-1)*cfg.mesh.edge_multiplier+2,:) = (xy_s(xt(j,1),:)+xy_s(xt(j,3),:))/2;
    xy((j-1)*cfg.mesh.edge_multiplier+3,:) = (xy_s(xt(j,2),:)+xy_s(xt(j,3),:))/2;
end
xy = unique(xy,'rows');

% Determine boundary nodes
idx_corners =  (xy_s(:,1)<x_min+eps & xy_s(:,2)<y_min+eps) | (xy_s(:,1)<x_min+eps & xy_s(:,2)>y_max-eps)  |  (xy_s(:,1)>x_max-eps & xy_s(:,2)>y_max-eps)  | (xy_s(:,1)>x_max-eps & xy_s(:,2)<y_min+eps) ;
xy_s(idx_corners,:)= [];
idx_b_in = xy_s(:,1)<x_min+eps;
boundary_in_s = xy_s(idx_b_in,:);
xy_s(idx_b_in,:)=[];
idx_b_y =  xy_s(:,2)>y_max-eps | xy_s(:,2)<y_min+eps ;
boundary_y_s = xy_s(idx_b_y,:);
xy_s(idx_b_y,:)=[];
idx_b_out =  xy_s(:,1)>x_max-eps ;
%
boundary_out_s = xy_s(idx_b_out,:);
xy_s(idx_b_out,:)=[];
idx_b_c =  xy_s(:,1).^2+xy_s(:,2).^2<(radius+eps)^2 ;
boundary_c_s = xy_s(idx_b_c,:);
xy_s(idx_b_c,:)=[];

boundary_s = [boundary_y_s;boundary_out_s;boundary_in_s;boundary_c_s];

% Process boundary nodes for velocity grid
idx_b_in = xy(:,1)<x_min+eps;
boundary_in = xy(idx_b_in,:);
xy(idx_b_in,:)=[];
idx_b_y =  xy(:,2)>y_max-eps | xy(:,2)<y_min+eps ;
boundary_y = xy(idx_b_y,:);
xy(idx_b_y,:)=[];
idx_b_out =  xy(:,1)>x_max-eps ;
%
boundary_out = xy(idx_b_out,:);
xy(idx_b_out,:)=[];
idx_b_c =  xy(:,1).^2+xy(:,2).^2<(radius+eps)^2 ;
boundary_c = xy(idx_b_c,:);
xy(idx_b_c,:)=[];

boundary = [boundary_y;boundary_out;boundary_in;boundary_c];

if doPlot
    figure;

    scatter(xy(:,1),xy(:,2),'k.'); hold on; axis square;
    scatter(boundary_in(:,1),boundary_in(:,2),'b+');
    scatter(boundary_y(:,1),boundary_y(:,2),'r+');
    scatter(boundary_out(:,1),boundary_out(:,2),'b+');
    scatter(boundary_c(:,1),boundary_c(:,2),'m+');

    axis equal;
end

% Combine interior and boundary nodes: xy1 (V-grid), xy1_s (P-grid)
xy1 = [xy; boundary];
xy1_s = [xy_s; boundary_s];
x1 = xy1(:,1);
y1 = xy1(:,2);
x1_s = xy1_s(:,1);
y1_s = xy1_s(:,2);
x0 = xy(:,1);
y0 = xy(:,2);
x0_s = xy_s(:,1);
y0_s = xy_s(:,2);

x_min_dist = cfg.distances.x_min;
x_max_dist = cfg.distances.x_max;
y_min_dist = cfg.distances.y_min;
y_max_dist = cfg.distances.y_max;
r_dist     = a2*cfg.mesh.edge_multiplier;

%% Determine local stencils
k = cfg.rbf.stencil_size_main;
[Nearest_Idx] = nearest_interp(xy1,xy1,k);

k = cfg.rbf.stencil_size_main;
[Nearest_Idx_s] = nearest_interp(xy1_s,xy1_s,k);

%%   differentiation: P-grid to P-grid
[D_s_all] = RBF_PHS_FD_all(xy_s,xy1_s,Nearest_Idx_s(1:length(xy_s),:),cfg.rbf.order_main,cfg.rbf.poly_degree_main,cfg.rbf.laplacian_order);
L_s = D_s_all{3};

%%  B.C.s for P-grid
[Nearest_Idx_b_c] = nearest_interp(boundary_c_s,xy_s,cfg.rbf.stencil_size_boundary_cylinder);
Nearest_Idx_b_c = [(length(xy1_s)-length(boundary_c_s)+1:length(xy1_s))' Nearest_Idx_b_c];

D = RBF_PHS_FD_all(boundary_c_s,xy1_s,Nearest_Idx_b_c, cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dn1_b_s = ((boundary_c_s(:,1)).*D{1}+(boundary_c_s(:,2)).*D{2})./radius;

[Nearest_Idx_b_y] = nearest_interp(boundary_y_s,[xy_s; xy1_s(length(xy_s)+length(boundary_y_s)+1,:)],cfg.rbf.stencil_size_boundary_wall);
Nearest_Idx_b_y = [(length(xy_s)+1:length(xy_s)+length(boundary_y_s))' Nearest_Idx_b_y];

D = RBF_PHS_FD_all(boundary_y_s,xy1_s,Nearest_Idx_b_y,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dy_b_s = D{2};

[Nearest_Idx_b_x] = nearest_interp(boundary_in_s,xy1_s(1:length(xy_s)+length(boundary_y_s),:),cfg.rbf.stencil_size_boundary_wall);
Nearest_Idx_b_x = [(length(xy1_s)+1-length(boundary_c_s)-length(boundary_in_s):length(xy1_s)-length(boundary_c_s))' Nearest_Idx_b_x];

D = RBF_PHS_FD_all(boundary_in_s,xy1_s,Nearest_Idx_b_x,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dx_in_s = D{1};

[Nearest_Idx_b_x] = nearest_interp(boundary_out_s,xy1_s(1:length(xy_s)+length(boundary_y_s),:),cfg.rbf.stencil_size_boundary_outlet);
Nearest_Idx_b_x = [(length(xy_s)+1+length(boundary_y_s):length(xy_s)+length(boundary_y_s)+length(boundary_out_s))' Nearest_Idx_b_x];

D = RBF_PHS_FD_all(boundary_out_s,xy1_s,Nearest_Idx_b_x, cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dx_out_s = D{1};

% Cylinder boundary conditions
L_bc_c = zeros(length(boundary_s), length(xy1_s));
Idx_boundary_c = length(xy1_s)-length(boundary_c_s)+[1:length(boundary_c_s)];
L_bc_c(Idx_boundary_c-length(xy_s),:) = Dn1_b_s;

% Outlet boundary conditions  
L_bc_out = zeros(length(boundary_s), length(xy1_s));
Idx_boundary_out = length(xy_s) + length(boundary_y_s) + [1:length(boundary_out_s)];
L_bc_out(Idx_boundary_out-length(xy_s),:) = Dx_out_s;

% Wall boundary conditions
L_bc_y = zeros(length(boundary_s), length(xy1_s));
Idx_boundary_y = length(xy_s) + [1:length(boundary_y_s)];
L_bc_y(Idx_boundary_y-length(xy_s),:) = Dy_b_s;

% Inlet boundary conditions
L_bc_in = zeros(length(boundary_s), length(xy1_s));
Idx_boundary_in = length(xy1_s) - length(boundary_c_s) - length(boundary_in_s) + [1:length(boundary_in_s)];
L_bc_in(Idx_boundary_in-length(xy_s),:) = Dx_in_s;

% Assemble pressure system matrix
L1 = [L_s(1:length(xy_s),:); L_bc_c+L_bc_out+L_bc_y+L_bc_in];

% Add regularization for pressure system
L1 = [L1 [ones(length(xy_s),1); zeros(length(boundary_s),1)]; [ones(1,length(xy_s)) zeros(1,length(boundary_s))] 0];

[LL,UU,pp,qq,rr]    = lu(L1);
L_inv_s = @(v) (qq*(UU\(LL\(pp*(rr\(v))))));

%%
[Nearest_Idx_interp_21_c] = nearest_interp(boundary_c,xy1_s,cfg.rbf.stencil_size_boundary_outlet);
[D0_21_all_c] = RBF_PHS_FD_all(boundary_c,xy1_s,Nearest_Idx_interp_21_c,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);

D0_21_x_c =  D0_21_all_c{1};
D0_21_y_c =  D0_21_all_c{2};

%%  differentiation : P-grid to V-grid
[Nearest_Idx_interp_21] = nearest_interp(xy1(1:length(xy)+length(boundary_y)+length(boundary_out),:),xy1_s,cfg.rbf.stencil_size_boundary_outlet);

% low-order PHS-RBFs and polynomials near boundaries
[D0_21_all] = RBF_PHS_FD_all(xy1(1:length(xy)+length(boundary_y)+length(boundary_out),:),xy1_s,Nearest_Idx_interp_21,cfg.rbf.order_interpolation_low,cfg.rbf.poly_degree_main,cfg.rbf.poly_degree_interpolation_low);

D0_21_x = D0_21_all{1};
D0_21_y = D0_21_all{2};

Idx = (1:length(xy))';

Nearest_Idx_nc = xy(:,1).^2+xy(:,2).^2>(radius+r_dist)^2  & xy(:,1)>x_min+x_min_dist & xy(:,2)<y_max-y_max_dist & xy(:,2)>y_min+y_min_dist & xy(:,1)< x_max-x_max_dist  ;
Nearest_Idx_nc = Idx(Nearest_Idx_nc);
xy_nc = xy1(Nearest_Idx_nc,:);

[D0_21_all_nc] = RBF_PHS_FD_all(xy_nc,xy1_s,Nearest_Idx_interp_21(Nearest_Idx_nc,:),cfg.rbf.order_interpolation_high,cfg.rbf.order_interpolation_high_poly,cfg.rbf.poly_degree_interpolation_high);

D0_21_x(Nearest_Idx_nc,:) =  D0_21_all_nc{1};
D0_21_y(Nearest_Idx_nc,:) =  D0_21_all_nc{2};

% Differentiation operators: V-grid to P-grid
[Nearest_Idx_interp] = nearest_interp(xy_s,xy1,cfg.rbf.stencil_size_boundary_outlet);

% low-order PHS-RBFs and polynomials near boundaries
D0_12_all = RBF_PHS_FD_all(xy_s,xy1,Nearest_Idx_interp,cfg.rbf.order_boundary,cfg.rbf.order_interpolation_high_poly,cfg.rbf.derivative_order);

D0_12_x =  D0_12_all{1};
D0_12_y =  D0_12_all{2};

Idx = (1:length(xy_s))';

Nearest_Idx_nc = xy_s(:,1).^2+xy_s(:,2).^2>(radius+r_dist)^2  & xy_s(:,1)>x_min+x_min_dist & xy_s(:,2)<y_max-y_max_dist & xy_s(:,2)>y_min+y_min_dist & xy_s(:,1)< x_max-x_max_dist;

Nearest_Idx_nc = Idx(Nearest_Idx_nc);
xy_nc = xy1_s(Nearest_Idx_nc,:);

D0_12_all_nc = RBF_PHS_FD_all(xy_nc,xy1,Nearest_Idx_interp(Nearest_Idx_nc,:),cfg.rbf.order_interpolation_high,cfg.rbf.order_interpolation_high_poly,cfg.rbf.poly_degree_interpolation_high);

D0_12_x(Nearest_Idx_nc,:)  =  D0_12_all_nc{1};
D0_12_y(Nearest_Idx_nc,:)  =  D0_12_all_nc{2};

%   B.C.s for V-grid
[Nearest_Idx_b_y] = nearest_interp(boundary_y,xy,cfg.rbf.stencil_size_boundary_cylinder);
Nearest_Idx_b_y = [(length(xy)+1:length(xy)+length(boundary_y))' Nearest_Idx_b_y];

D = RBF_PHS_FD_all(boundary_y,xy1,Nearest_Idx_b_y,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dy_b_0 = D{2};

Dy_b = Dy_b_0;
Dy_b_1 = diag(Dy_b(:,Nearest_Idx_b_y(:,1)));
Dy_b(:,Nearest_Idx_b_y(:,1)) = zeros(length(boundary_y),length(boundary_y));

[Nearest_Idx_b_out] = nearest_interp(boundary_out,xy,cfg.rbf.stencil_size_boundary_cylinder);
 Nearest_Idx_b_out = [(length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out)))' Nearest_Idx_b_out];

D =  RBF_PHS_FD_all(boundary_out,xy1,Nearest_Idx_b_out,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dx_b_0 = D{1};

Dx_b = Dx_b_0;
Dx_b_1 = diag(Dx_b(:,Nearest_Idx_b_out(:,1)));
Dx_b(:,Nearest_Idx_b_out(:,1)) = zeros(length(boundary_out),length(boundary_out));

D_all = RBF_PHS_FD_all(xy1,xy1,Nearest_Idx,cfg.rbf.order_main,cfg.rbf.poly_degree_main,cfg.rbf.laplacian_order);

Dx = D_all{1};
Dy = D_all{2};
L0 = D_all{3};

Nearest_Idx_nb =  xy(:,1).^2+xy(:,2).^2<(radius+r_dist)^2 ;
Idx = (1:length(xy))';
Nearest_Idx_nb = Idx(Nearest_Idx_nb);

xy_nb = xy1(Nearest_Idx_nb,:);

D_all_nb = RBF_PHS_FD_all(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),cfg.rbf.order_near_cylinder,cfg.rbf.order_near_cylinder_poly,cfg.rbf.poly_degree_near_cylinder);

Dx(Nearest_Idx_nb,:)  =  D_all_nb{1};
Dy(Nearest_Idx_nb,:)  =  D_all_nb{2};
L0(Nearest_Idx_nb,:)  =  D_all_nb{3};

% low-order PHS-RBFs and polynomials near boundaries
Nearest_Idx_nb =  xy(:,1)>x_max-x_max_dist |xy(:,1)<x_min+x_min_dist*2| xy(:,2)> y_max-y_max_dist  | xy(:,2)<y_min+y_min_dist ;
Idx = (1:length(xy))';
Nearest_Idx_nb = Idx(Nearest_Idx_nb);
xy_nb = xy1(Nearest_Idx_nb,:);

D_all_nb = RBF_PHS_FD_all(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),cfg.rbf.order_near_boundary,cfg.rbf.poly_degree_main,cfg.rbf.poly_degree_near_boundary);

Dx(Nearest_Idx_nb,:)  =  D_all_nb{1};
Dy(Nearest_Idx_nb,:)  =  D_all_nb{2};
L0(Nearest_Idx_nb,:)  =  D_all_nb{3};

%% Time-stepping parameters
nu = cfg.simulation.viscosity;
dt = cfg.simulation.time_step;

% Set reduced number of time steps for CI environment
Nt = cfg.simulation.num_time_steps;
if isCI
    Nt = cfg.simulation.num_time_steps_ci; % Keep tests fast in CI
end

% 
L = L0;
L(length(xy)+1:end,:) = zeros(length(boundary),length(xy1));

L_I = speye(length(xy1))-dt*nu*L*cfg.schemes.crank_nicolson;

L_u = L_I;  L_v =L_I;

L_u((length(xy)+1:length(xy)+length(boundary_y)),:) =   Dy_b_0;
L_u ((length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out))),:) = Dx_b_0;
L_v ((length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out))),:) = Dx_b_0;

clear L_I L

[LL,UU,pp,qq,rr]    = lu(L_u);
L_u_inv = @(v) (qq*(UU\(LL\(pp*(rr\(v))))));

[LL,UU,pp,qq,rr]    = lu(L_v);
L_v_inv = @(v) (qq*(UU\(LL\(pp*(rr\(v))))));

clear L_u L_v

%% Initialize simulation
V0 = zeros(length(xy1),1);  % Initial v-velocity
U0 = ones(length(xy1),1);   % Initial u-velocity
U0(end-length(boundary_c)+1:end) = zeros(length(boundary_c),1);  % Zero on cylinder
W0 = [U0; V0];  % Combined velocity vector

W = zeros(length(xy1)*2,Nt+1);
W(:,1) = W0;

p0 = zeros(length(xy1_s),1);

L_B = length(boundary_c)+length(boundary_in);
L_B_y =  length(boundary_y);
L_W = length(xy);
L_B_c= length(boundary_c);

Idx_y_s = length(xy_s)+1: length(xy_s)+length(boundary_y_s);
Idx_in_s = length(xy1_s)+1-length(boundary_c_s)-length(boundary_in_s): length(xy1_s)-length(boundary_c_s);
Idx_out_s = length(xy_s)+1+length(boundary_y_s): length(xy_s)+length(boundary_out_s)+length(boundary_y_s);

W = zeros(length(xy1)*2,Nt+1);
W(:,1) = W0;

for j = 1 :Nt
    disp(['j = ' num2str(j)])
    if j <3
        %
        W(:,j+1)=W(:,j);
    else
        [W(:,j+1),p0] = NS_2d_cylinder_PHS(dt,nu,W(:,j-1),W(:,j),Dy,Dx,L_inv_s,L_u_inv,L_v_inv,L0,L_B,L_B_c,L_W,L_B_y,length(boundary_s),D0_12_x,D0_12_y,D0_21_x,D0_21_y,Dy_b,Dy_b_1,D0_21_x_c,D0_21_y_c,p0,W(:,1));
    end
    if isnan(W(1,j+1))
        break;
    end
end

W0 = W(:,end);

if doPlot
    figure('Name','1/Re = 1e-2');
    colormap(jet)

    j = Nt;
    U = W(1:length(xy1),(j-1)*1+1);
    V = W(length(xy1)+1:end,(j-1)*1+1);
    
    % Plot u-velocity field
    subplot(2,1,1);
    scatter(x1,y1,cfg.visualization.scatter_size*ones(length(xy1),1),1*(U-1),'.');
    axis equal, axis tight, hold on;
    xlim([x_min x_max]);
    ylim([y_min y_max]);
    yticks(cfg.visualization.plot_tick_y)
    xticks(cfg.visualization.plot_tick_x)
    title(['u']);
    ylabel('y');
    xlabel('x');
    shading interp;
    caxis([-cfg.visualization.color_axis_range cfg.visualization.color_axis_range]);

    % Plot v-velocity field
    subplot(2,1,2);
    scatter(x1,y1,cfg.visualization.scatter_size*ones(length(xy1),1),1*V,'.');
    axis equal, axis tight, hold on;
    shading interp;

    xlim([x_min x_max]);
    ylim([y_min y_max]);
    yticks(cfg.visualization.plot_tick_y)
    xticks(cfg.visualization.plot_tick_x)
    title(['v' ]);

    xlabel('x');
    ylabel('y');
    set(gca,'Ytick',[]);
    caxis([-cfg.visualization.color_axis_range cfg.visualization.color_axis_range]);

    drawnow;
end