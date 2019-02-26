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
%           p  = parameters for each microbial population stored in
%                'params'
function DEs = monodDE(t,IC,p)

% Starting concentrations of substrate, microbial mass 1, nitrogen,
% phosphorus, Acetate, and microbial mass 2
S = IC(1);
M1 = IC(2);
N = IC(3);
P = IC(4);
A = IC(5);
M2 = IC(6);
M3 = IC(7);
CH4= IC(8);

% Substrate Equation Diff EQ 
S1_by_M1=-M1*p(1).mu_max*S/(p(1).Ks+S)*N/(p(1).kN+N)*P/(p(1).kP+P);
S1_by_M2=-M2*p(2).mu_max*S/(p(2).Ks+S)*N/(p(2).kN+N)*P/(p(2).kP+P);
DEs(1,1) = S1_by_M1 + S1_by_M2; 
A_by_M1 = -S1_by_M1;
% Microbial population 1 Diff EQ
% DEs(2,1) =  mu*yield*S*M/(Ks+S)*N/(kN+N)-Kd*M;
DEs(2,1) =  -S1_by_M1*p(1).Y-p(1).Kd*M1;

% Nitrogen concentration Diff EQ
DEs(3,1) =  1/5*DEs(1,1); %Based on C:N ratio of biomass

% Phosphorus concentration Diff EQ
DEs(4,1) = 1/20*DEs(1,1);   %Assumed ratio of phosphorus biomass

% Acetate my M3
A_by_M3 = -M3*p(3).mu_max*A/(p(3).Ks+A)*N/(p(3).kN+N)*P/(p(3).kP+P);

% Acetate concentration Diff EQ
DEs(5,1) = A_by_M1 + A_by_M3;

% Microbial population 2 Diff EQ
DEs(6,1) = -S1_by_M2*p(2).Y-p(2).Kd*M2;
% DEs(6,1) = ; %New microbial population
% DE(7,1) = -DE(6,1); %CO2
DEs(7,1) = -A_by_M3*p(3).Y - p(3).Kd*M3;
DEs(8,1) = -A_by_M3;


% Microbial population 3 Diff EQ


end




