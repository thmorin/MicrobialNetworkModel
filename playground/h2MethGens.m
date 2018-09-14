function [co2V,h2V,ch4V,mu]=h2MethGens(h2Vmax,H2,H2_Km,CO2,CO2_Km,Y,ft,fm,fph)

%% HYDROGENOTROPHIC METHANOGENS
% Catabolic stoichiometry: CO2 + 4H2 -> CH4 + 2*H2O
h2V =h2Vmax * H2/(H2 + H2_Km) * CO2/(CO2 + CO2_Km)*ft*fm*fph;
co2V=h2V/4;
ch4V=h2V/4;

mu=Y * co2V;

% mu =  mu_max * H2/(H2 + H2_Km) * CO2/(CO2 + CO2_Km)*ft*fm*fph;
% co2V =   mu/Y;
% h2V  =   co2V;%Note I took out some stoichiometry here. Maybe that's not the right way to handle it?
% ch4V =   co2V;