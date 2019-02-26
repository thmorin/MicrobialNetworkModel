%   Written by Steve Colquhoun 10/21/2018
%   Fall 2018 Research for Dr. Tim Morin
%   SUNY College of Environmental Science and Forestry

%   General form of speciation function

%   Inputs: Total_conc = total concentration            Units: mol/L
%                   pH = pH of system                   Units: ---
%             pka_vect = vector of pka values in ascending order           

%   Outputs: concs_out = species concentrations         Units: mol/L

function[concs_out] = poly_species(Total_conc, pH, pKa_vect)
% clear;close;Total_conc=100;pH=7;pKa_vect=[6.33;10.7;4.32];
%   Convert pKa to Ka
Ka_vect = 10.^-(pKa_vect);

%   Hydrogen ion concentration                          Units: mol/L
H = 10^-pH;

n = length(pKa_vect);
%   Preallocation of output
concs_out = nan(n+1,1);

buildingSum=1;
step1=1;
for i=1:n
    step1=step1*Ka_vect(i)/H;
    buildingSum = buildingSum+step1;
end
concs_out(1)   = Total_conc/buildingSum;
for i = 2:n+1
    concs_out(i) = Ka_vect(i-1)*concs_out(i-1)/H;
end
end






