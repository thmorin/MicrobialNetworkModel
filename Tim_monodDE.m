%   Written by Steve Colquhoun 12/01/2018
%   Updated by Steve Colquhoun 02/05/2019
%   Spring 2019 Research with Dr. Tim Morin
%   SUNY College of Environmental Science and Forestry

%   Script for calculating microbial populations using Monod Kinetics

%   Required parameters:
%   Ks     = half saturation constant       Units: mass/volume
%   Kd     = rate of decay                  Units: 1/time
%   mu_max = maximum growth rate            Units: 1/time
%   Yf     = Yield Factor                   Units: ---
%   kN     = Nitrogen constant              Units: 
%   kP     = Phosphorus constant            Units:

%   Function for calling into ode solver
%   Inputs: t  = time
%           IC = vector of inital conditions, order shown in function 
%           m  = microbial parameters
%           s  = substrate related parameters
%           p  = model domain parameters
function DEs = Tim_monodDE(t,IC,m,s,p)

IC=reshape(IC,5,8);
% Zero flux top boundary condition
IC=[IC(1,:);IC];
dz=[p.dz(1);p.dz];
% Zero flux top boundary condition
IC=[IC;IC(end,:)];
dz=[dz;p.dz(end)];

% Starting concentrations of substrate, microbial mass 1, nitrogen,
% phosphorus, Acetate, and microbial mass 2
DOC = IC(:,1);
N   = IC(:,2);
P   = IC(:,3);
A   = IC(:,4);
CH4 = IC(:,5);

M1  = IC(:,6);
M2  = IC(:,7);
M3  = IC(:,8);


%% Substrates
for i=2:length(DOC)-1 %Loop down water column
    % Dissolved Organic Carbon (DOC)
    S1_by_M1  = -M1(i)*m(1).mu_max*DOC(i)/(m(1).Ks+DOC(i))*N(i)/(m(1).kN+N(i))*P(i)/(m(1).kP+P(i));
    A_by_M1   = -S1_by_M1;
    
    S1_by_M2  = -M2(i)*m(2).mu_max*DOC(i)/(m(2).Ks+DOC(i))*N(i)/(m(2).kN+N(i))*P(i)/(m(2).kP+P(i));
    CO2_by_M2 = -S1_by_M2;

    DOC_diff  = Diffusion(DOC(i-1:i+1),dz(i-1:i),s(1).D);
    DEs(i-1,1) = S1_by_M1 + S1_by_M2 + DOC_diff; %DOC
    
    % Acetate
    A_by_M3 = -M3(i)*m(3).mu_max*A(i)/(m(3).Ks+A(i))*N(i)/(m(3).kN+N(i))*P(i)/(m(3).kP+P(i));
    A_diff     = Diffusion(A(i-1:i+1),dz(i-1:i),s(4).D);
    
    % Nitrogen
    N_diff     = Diffusion(N(i-1:i+1),dz(i-1:i),s(2).D);
    DEs(i-1,2) =  1/5*(S1_by_M1 + S1_by_M2 + A_by_M3) + N_diff; %Based on C:N ratio of biomass

    % Phosphorus
    P_diff     = Diffusion(P(i-1:i+1),dz(i-1:i),s(3).D);
    DEs(i-1,3) = 1/20*(S1_by_M1 + S1_by_M2 + A_by_M3) + P_diff; %Based on C:P ratio of biomass

    DEs(i-1,4) = A_by_M1 + A_by_M3 + A_diff;

    %CH4
    CH4_diff     = Diffusion(CH4(i-1:i+1),dz(i-1:i),s(5).D);
    DEs(i-1,5) = -A_by_M3 + CH4_diff;
    %% Microbial populations
    DEs(i-1,6) =  -S1_by_M1*m(1).Y-m(1).Kd*M1(i);  %Heterotrophs
    DEs(i-1,7) =  -S1_by_M2*m(2).Y-m(2).Kd*M2(i);  %Acetogens
    DEs(i-1,8) =  - A_by_M3*m(3).Y-m(3).Kd*M3(i); %Acetoclastic methanogens
end
DEs=reshape(DEs,5*8,1);



