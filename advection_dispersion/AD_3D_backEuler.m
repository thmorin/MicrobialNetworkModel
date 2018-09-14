function [C_out]=AD_3D_backEuler(C_in,U,D,Q,dt,X,Y,Z,Q)
% C - 3D - L X M X N
% U - 3D - fluid speed in X direction
% V - 3D - fluid speed in Y direction
% W - 3D - fluid speed in Z direction
% X - 3D - center of grid cell, X direction
% Y - 3D - center of grid cell, Y direction
% Z - 3D - center of grid cell, Z direction
% dt - time step, scalar
% D  - 2X1 - rate of diffusion. 1 in vertical, 2 in lateral


%Note: q not handled correctly yet. Assume 1st order kinetics? Might need
%this per compound

%Solving:
% dC/dt = -U*dC/dX - V*dC/dY - W*dC/dZ + D_L*d2C/dX2 + D_L*d2C/dY2 + d_V*d2C/dZ2 + Q

for i=2:
    for j=2:
        for k=2:
            K=C(i,j,k)/dt %<---problem here
            
            M= -U/(X(i+1,j,k)-X(i,j,k)) + D(1)/(X(
            
        end
    end
end

% Impose boundary conditions



% Solve