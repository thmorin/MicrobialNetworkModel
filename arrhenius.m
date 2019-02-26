%   Written by Steve Colquhoun 11/01/2018
%   Fall 2018 Research for Dr. Tim Morin
%   SUNY College of Environmental Science and Forestry

%   Arrhenius function for modeling temperature dependence of microbial
%   populations

%   Inputs:  k0 = rate constant at a reference temperature  Units: Variable
%            T0 = reference temperature                     Untis: K
%            T  = temperature of system                     Units: K
%            Ea = activation energy                         Units: kJ/mol

%   Outputs: kt = rate constant at system temperautre       Units: Variable

function[kt] = arrhenius(k0,T0,T,Ea)

%   Universal gas constant                                  Units: kJ/mol*K
R  = 8.314e-3;

%   Arrhenius equation
kt = k0*exp(Ea./R.*((T0^-1)-(T.^-1)));

end