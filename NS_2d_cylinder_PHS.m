function [W3,p] = NS_2d_cylinder_PHS(dt,nu,W1,W2,Dy,Dx,L_inv,L_u_inv,L_v_inv,L0,L_B,L_B_c,L_W,L_B_y,L_B_S,D0_12_x,D0_12_y,D0_21_x,D0_21_y,Dy_b,Dy_b_1,D0_12_x_c,D0_12_y_c,p0,W0)


L_W2 =  length(W2);

U    =  W2(1:L_W2/2);
V    =  W2(L_W2/2+1:end);

U_1    =  W1(1:L_W2/2);
V_1    =  W1(L_W2/2+1:end);


% U0    =  W0(1:L_W2/2);
% V0    =  W0(L_W2/2+1:end);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% advection: adam-bashforth method




H_U     =  (-U.*(Dx*(U))-V.*(Dy*(U)));

H_V     =  (-U.*(Dx*(V))-V.*(Dy*(V)));


H_U_1   =  (-U_1.*(Dx*(U_1))-V_1.*(Dy*(U_1)));

H_V_1   =  (-U_1.*(Dx*(V_1))-V_1.*(Dy*(V_1)));




U1      =  U + (1*dt) * (3*H_U-H_U_1)/2;

V1      =  V + (1*dt) * (3*H_V-H_V_1)/2;




% enforce B.C.s for velocities



 U1(end-L_B+1:end-0) = U(end-L_B+1:end-0);


U1(L_W+1:L_W+L_B_y) = -(Dy_b*U1)./Dy_b_1;


% V1(L_W+1:end-0) = V(L_W+1:end-0);


 V1(L_W+1:L_W+L_B_y) =  V(L_W+1:L_W+L_B_y);

 V1(end-L_B+1:end-0) = V(end-L_B+1:end-0);
%  
%  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% viscousity: Crank-Nicolson

% 
U2   =  U1+(1*dt)*nu*(L0*U)/2;

V2   =  V1+(1*dt)*nu*(L0*V)/2;


U2(end-L_B+1:end) = U(end-L_B+1:end);
V2(end-L_B+1:end) = V(end-L_B+1:end);


% % 
U2(L_W+1:end-L_B) = zeros(L_W2/2-L_W-L_B,1);
V2(L_W+L_B_y+1:end-L_B) = zeros(L_W2/2-L_W-L_B-L_B_y,1);


U2(L_W+1:L_W+L_B_y) =  zeros(L_B_y,1);
V2(L_W+1:L_W+L_B_y) =  zeros(L_B_y,1);


[L_B_c,~]  =  size(D0_12_x_c);

 U2(end-L_B_c+1:end)   =  D0_12_x_c*p0;
 V2(end-L_B_c+1:end)   =  D0_12_y_c*p0;


U2   =  L_u_inv(U2);
V2   =  L_v_inv(V2);

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pressure correction 


   F0    =  (D0_12_x*(U2-0)+D0_12_y*(V2-0));

   
  F = [F0; zeros(L_B_S+1,1)];

  

  p = L_inv(F);

  
  p = p(1:(length(F0)+L_B_S));  % regularization
  F = F(1:(length(F0)+L_B_S));  % regularization
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    U3   =   (U2-0) -  [(D0_21_x*p ); zeros(L_B,1)];
    
    V3   =   (V2-0) -  [(D0_21_y*p) ; zeros(L_B,1)];
    
    
% enforce B.C.s for velocities


U3(end-L_B+1:end)   = U(end-L_B+1:end);
U3(L_W+1:L_W+L_B_y) = -(Dy_b*U3)./Dy_b_1;


V3(L_W+1:L_W+L_B_y) = V(L_W+1:L_W+L_B_y);
V3(end-L_B+1:end)   = V(end-L_B+1:end);
   
% V3(L_W+1:end) = V(L_W+1:end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% velocities at the following time-step


W3   =   [U3;V3];



end
