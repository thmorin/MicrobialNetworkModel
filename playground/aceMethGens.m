function [aceV,HV,CH4V,CO2V,acemu]=aceMethGens(aceVmax,ace,aceY,aceKm,ft,fm,fph)

%% ACETOCLASTIC METHANOGENS
%  Catabolic stoichiometry:  C2H3O2- + H+ -> 0.5*CH4 + 0.5*CO2

aceV=aceVmax * ace / (ace + aceKm)*ft*fm*fph;
  HV=aceV;               %Based on stoichiometry
CH4V=aceV*0.5;        %Based on stoichiometry
CO2V=aceV*0.5;        %Based on stoichiometry

acemu=aceY*2 * aceV;  %Adjusting for Y units of [mol C (mol-ace-C)^-1] with the /2