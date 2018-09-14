%% State variables (initial conditions
v.Ts             = 300.0;   %[K] - Soil temperature
v.pH             =   7.0;   %[-] %-5.69 per pH unit <- Energetics? This got unlinked from something when I was shifting things around
v.pOH            = 14-v.pH;

%% Substrates
% Initial conditions
s.H.tot =10^-(   v.pH);
s.OH.tot=10^-(14-v.pH);
s.ace.tot=   10; %[mmol m^-3] - Initial acetate concentration
s.NH4.tot=  100; %[mol m^-3] - Initial bioavilable nitrogen
s.PO3.tot=  100; %[mol m^-3] - Initial bioavailable phosphorus
s.CH4.tot=    0; %[mol m^-3] - Initial methane
s.H2 .tot=  1.8; %[mol m^-3] - Initial H2 concentration
s.CO2.tot=    0; %[mol m^-3] - Initial CO2 concentration

%% Microbial initial conditions
m.CH4_ace.tot    =  350;              %[mmol m^-3] %Pulled 350 from one of Kelly's cores for OWC
m.CH4_H2 .tot    =  350;        %[mol-CH2O-C m^-3]
m.CH4_tro.tot    =  2.0*10^-3;        %[mol-CH2O-C m^-3]