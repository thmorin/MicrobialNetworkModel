function [aceV,HV,CH4V,CO2V,acemu,Yapparent]=aceMethGensTEST ...
    (aceVmax,ace,aceKm,ft,fm,fph,GfCH4,GfAce,GfCO2,Gfbio,Rm,Rw,...
    CO2,CH4,NH3,PO4,H,R,T,CNrat,CPrat,Xch4)

%% ACETOCLASTIC METHANOGENS
%  Catabolic stoichiometry:  C2H3O2- + H+ -> 0.5*CH4 + 0.5*CO2

aceV=aceVmax * ace / (ace + aceKm)*ft*fm*fph;


dGf_cat = GfCO2 + GfCH4 - GfAce;
dGf_an = 2*Gfbio - GfAce;

if CO2~=0 && CH4~=0 && ace~=0 && H~=0
    dG_cat = dGf_cat + R*T*log(CO2*CH4/ace/H);
    dG_an  = dGf_an  + R*T*log(1/ace/H/NH3^CNrat/PO4^CPrat);
elseif CO2==0 || CH4==0
    dG_cat = dGf_cat;
    dG_an  = dGf_an  + R*T*log(1/ace/H);
elseif ace==0
    dG_cat = dGf_cat;
    fprintf('Ran out of acetate\n');
elseif H==0
    error('H hit 0. Something big is wrong\n');
end

% N=(Rm + Rw + aceV*dG_cat)/...
%     (aceV*dG_cat - aceV*dG_an);

N=Rm*Xch4+aceV*Xch4*(dG_cat - Rw) / ( aceV*Xch4*(dG_cat - Rw) - aceV*dG_an);
acemu=N*aceV;
Yapparent=N;

  HV=(1-N)*aceV;            %Based on stoichiometry
CH4V=(1-N)*aceV*0.5;        %Based on stoichiometry
CO2V=(1-N)*aceV*0.5;        %Based on stoichiometry

% acemu=aceY*2 * aceV;  %Adjusting for Y units of [mol C (mol-ace-C)^-1] with the /2


%Apparent Y? Should be able to comapre to measured values