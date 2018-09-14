function []=MM_CH4_A(m,ace,ch4)
%Written by Tim Morin, 6/2018
%
%Accepts ohjects from main for the acetoclastic methanogenic microbial population
%the acetate substrate pool, the methane substrate pool
%
%What is the N and P source for this microbial population?

% Acetate

% Q (I assume monomers) -> 0.67*ace + 0.33*CO2 + 0.11*H2    <--Grant old supplement, G.2



%P -> 0.67*A + 0.33*CO2 + 0.11*H2
%CH4 + 4*O2 -> CO2 + 1.5H2O + 0167*H+
%CH4 + 1.33*O2 -> CH2O + 0.167*H+ <- formation of biomass
%CH2O + 2.67*O2 -> CO2+1.5*H2O%

%CH3COO- + SO4 2- + 3*H + -> 2* Co2 + H2S + 2*H2O ----- degG=-57.5 kJ --- pp 387 Brocks

%4*H2 + So4 2- + H+ -> HS- + 4*H2O ------delG = -152 kJ --- pp387 Brock's

%4*HPO3 - + SO42- + H+ -> 4*HPO4 2- + HS- --- delG=-364 kJ ----Brock's 388

S0 + 2H -> H2S% Sulfate reducers do this? ----- pp388\\

%3*H2 + H+ + 2*HCO3 - -> Ch3CoO- + 4*H2O --- degG=-105 kJ --pp388

% v---sugars
%C6H12O6 -> 3CH3COO- + 3H+

%2*HCO3- + 4*H2 + H+ -> CH3COO- + 4*H2O

% Fermentation - anarobic catabolism where organic compound is both electron donor and electron acceptor
% respiration - compound is oxidized with o2 or an o2 subs as TEA

%2*pyruvate- -> 2*acetate- + 2*CO2 + 4H (note, not an H+, just an H), pp 389

%Glutamate, Glatamine
% Glutamine - Gf = -529.7 kJ/mol


%Sugars, amino acids, nucleotides, fatty acids


%Amino acids - 4 families - glutamat, aspartate, alanine, serine, armomatic
Glutamine
Proline



%And there are two separate cycles that handle them
%But they can get them preformed from the environment, too?
%Could roughly come up with a representative biomass

%Anabolism - pp108 anabolism description begins

2*Pyruvate + 4*H ._ 3*acetate- + H+

% Anabolism
sugarsamino acids
nucleotides
fatty acids
---storage - glycogen/starch
---polysaccharides
---activated glucose?


%Oxidation of methane to CO2
%CH4 + 4*O2 -> CO2 + 1.5 H2O + 0.167*H+                                       (G.22) - doesn't balance?
%CH4 + 2*O2 -> CO2 + 2*H2O   <---Tim's balance of it. Why is this not the correct version? Maybe anabolic/catabolic?

%Oxidation of CH4 to microbial storage C
% CH4 + 1.33*O2 -> CO2O + 0.167*H+

% Oxidatin of microbial storage C to CO2
% v--microbial storage? good biomass formula?
% CH2O-C + 2.67*O2  ->  CO2-C + 1.5*H2O                                         (G.24)



%Cellulose
%Starch
%Lactate
%Succinate
%Pectin



% Additinal oxidatin of C by denitrificers
%R = 0.429*R_NO3 + 0.429*R_NO2 + 0.214*RN2O


%Relevant free energies of formation                                            (kJ mol^-1)
G.CH4 =  -50.75;
G.ace = -369.41;
G.glu = -917.22;

G_fuckyou = Gf + R*T*log(ch4.d / ace.d)

%Energetic yield





if G_fuckyou>0
  disp('Reaction unfavorable')
end


endfunction
