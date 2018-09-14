%% General parameters
c.t0      =     0.00;                       %[d]
c.dt      =     10/60/60/24;                %[d^-1] - Essentially 1 [s] converted to [d^-1] units. Time step - trying to make sure it is not dt=1 since that might make bugs
c.out_rate=      1      /24;                %[d] - Output interval. 1 hour by default
% c.t_final =     2.00;                     %[d]
c.t_final = 2;
c.nt      = (c.t_final-c.t0)/c.dt;          %Final size of all variables
c.time    = linspace(c.t0,c.t_final,c.nt);

c.F       =    96.48; %[kJ/V]
c.R       =     8.29; %[J/mol/Kelvin]

c.t.y  = 0.081; %Grant supplement  - said it's to give a Q10 for f_tm of 2.25. USed in temperature adjustment formula

c.CNrat=106/16;
c.CPrat=106/ 1;

%Pulled these from Grant's 2002 Methane paper
c.m.Rm=0.0115; %gC g-biomass^-1 h^-1 - Specific maintenance respiration at 25C

c.m.Rm=0; %Maintenance energy. Setting to 0 for test purposes
c.m.Rw=0; %Waste energy. Setting to 0 for test purposes.

c.m.D   = 0.1;   %Microbial generic death rate. Value not specified in Grant 2002 or supplement. Probably use the more specific ones when available

%Energetic constants - pulled form elsewhere. Liu 2007 maybe
c.e.drG0x = -3500; %kJ C-mol^-1 - Dissipitive energy - May actually be a function of the maximum growth rate of an organism. High growth rate organisms may grow less efficeintly in order to achieve a higher mu value. More applicable to high substrate environments
c.e.Gfa = - 515; %kJ C-mol^-1 - Energy released upon combusting dry biomass? Probably in reference to CO2 - Liu et al 2007, section "3. Experimental"

%I am treating the above like a free energy of formation, and you'd then
%want to look at how it converts different parts of its food sources for elements
%Citation for the above: "E.H. Battley, R.L. Putnam, J. Boerio-Goates, Thermochim. Acta 298 (1/2) (1997) 37-46."
c.e.dcGx  = - 474; %kJ c-mol^-1 - Adjusting the above for the conversion of NH4+ into biomass. This is also not clear how this adjustment was made
c.e.dcGx  = - 500; %What Liu actually used
c.e.dGw_normal     = - 500; %kJ C-mol^-1 - From Liu et al 2007, section "4. Results". May be improved with a degree of reduction calculation
c.e.dGw_autotrophs = -3500; %kJ C-mol^-1 - Same as above. Originally given by Heijnen and Dijken model


% Brock - Biology of Microorganisms, 2012 13th ed, appx 1 -----    C2H3O2- + H+ -> 0.5*CH4 + 0.5*CO2 <--Acetoclastic methanogenesis
s.CH4.Gf =  -50.75;   %kJ/mol
s.CO2.Gf = -394.4;    %kJ/mol
s.ace.Gf = -369.41;   %kJ/mol
s.H2.Gf  =    0.00;   %kJ/mol
% Below are pulled from random internet sources - probably need a stronger source
s.bio.Gf = -110.00;   %kJ/mol - Biomass - CH2O
s.NH3.Gf =  -16.50;   %kJ/mol
s.PO4.Gf = -1119.2;   %kJ/mol

%Pulled these from Grant's 2002 Methane paper
s.dGc    = -37.50; % kJ g CH2O-C^-1
s.dGf    =  -4.43; % kJ g glucose-C^-1 - Brock and Madigan [1991], Schink [1997]
s.dGh    =  -0.27; % kJ g CO2-C^-1     - Brock and Madigan [1991]
s.dGm    =  -1.03; % kJ g acetate-C^-1 - Brock and Madigan [1991], Schink [1997]
s.dGt    =  -9.45; % kJ g CH4-C^-1     - Brock and Madigan [1991]
s.CO2.EC =  75.00; % kJ g CO-C2^-1 Anthony [1982]
%% Microbes
m.bio.EG      = 23.5;    %kJ g CH4-C^-1 Anthony [1982]
m.CH4.EM      = 25;      %kJ g org.C^-1 Anthony [1982]

%% Aerobic heterotrophs

%% Acetogenesis from DOC
% Constants pulled from Xu et al 2015 JGR:B, Microbial functional group-based CH4 model, which are themselves pulled form other sources - Table 1
m.ace_DOC.ace.Vmax = 0.035; %[mmol m^-3 h^-1] - Max rate of production of acetic acid from available carbon - Smith and Mah [1966]
%^--The above may actually be for ace_CO2

%% Homoacetogenesis
% Constants pulled from Xu et al 2015 JGR:B, Microbial functional group-based CH4 model, which are themselves pulled form other sources - Table 1
m.ace_CO2.H2. Vmax = 0.01;  %[mmol-ace g^-1 h^-1] - Max rate of conversion of H2 and CO2 to acetic acid - Conrad [1989]
m.ace_CO2.CO2.Km   = 0.040; %[mmol m^-3] - Half-saturation const. of acetic acide miner. to Co2 - Stoichiometry theory (no other reference given?)
m.ace_CO2.ace.Km     = 0.0165; %[umol m^-3] - Half saturation coefficient of conversin of H2 and CO2 to acetic acid - Conrad [1989]
m.ace_CO2.CO2.Km     = 0.00825; %[umol m^-3] - Assumed to be half  of H2 basec on stoichiometry - No other reference given. Notation wrong?
%^--- This last one doesn't make sense because there are two of these. Must have misunderstood

% m.ace_CO2.H2. Vmax_PDF = 1.1:.1:12.4; %Range specified in Xu. I do not understand why it doesn't bracket the value above yet
% m.ace_CO2.CO2.Km_PDF = 0.02:0.1:4 ; %Range specified in Xu
% m.ace_CO2.ace.Km_PDF = 6:27 ; %Range specified in Xu. Do not understand why it doesn't bracket the above

%% Aerobic methanotrophy
%Pulled these from Grant's 2002 Methane paper
m.CH4_tro.CH4.Km  = 3*10^-3; %[gC m^-3] - Conrad [1984] ???
m.CH4_tro.CH4.k2  = m.CH4_tro.CH4.Km * 0.25;
m.CH4_tro.mu_max  =-9999;    %no value given in Grant

% Constants pulled from Xu et al 2015 JGR:B, Microbial functional group-based CH4 model, which are themselves pulled form other sources - Table 1
% m.CH4_tro.mu_max  = 0.15;    %[d^-1] - Growth rate of aerobic methanotroph - Servais et al [1985]
% m.CH4_tro.Kd      = 0.005;   %[d^-1] - Death rate of aerobic methanotroph - Servais et al [1985]
% m.CH4_tro.Y       = 0.4;     %[mol C (mol-ace-C)^-1] - Growth efficiency of aerobic methanotroph - Kettunen [2003]
% m.CH4_tro.CH4.Km  = 0.0025;  %[mmol L^-1] - Half-saturation coefficient of CH4 oxidation for CH4 concentration - Kettunen [2003]
% m.CH4_tro.O2.Km   = 0.5;     %[mmol L^-1] - Half saturation coefficient of CH4 oxidation for O2 concentration - Kettunen [2003]

%% Acetoclastic methanogens - NEEDED FOR FIRST PASS
%  Catabolic stoichiometry:  C2H3O2- + H+ -> 0.5*CH4 + 0.5*CO2
%Pulled these from Grant's 2002 Methane paper
% m.CH4_ace.ace.Km  = 12.00;  %[gC m^-3] - Smith and Mah [1978] and Zehnder et al [1980]
% m.CH4_ace.mu_max  =  0.20;  %[g C g microbial C h^-1] - Smith and Mah [1980]
%The below is a rquired parameter in ECA kinetics, but so far as I know has no published values - Good question to ask Kelly and Bill
% m.CH4_ace.ace.k2  = m.CH4_ace.ace.Km * 0.25;

% Constants pulled from Xu et al 2015 JGR:B, Microbial functional group-based CH4 model, which are themselves pulled form other sources - Table 1
m.CH4_ace.ace.Km  =  5.000; %[mmol m^-3] Half-saturation coefficient - [Kettunen 2003]
m.CH4_ace.mu_max  =  0.035; %[d^-1] - Growth rate of acetoclastic methanogens - Servais et al [1985]
% m.CH4_ace.Kd      =  0.035; %[d^-1] - Death rate of acetoclastic methanogens - Servais et al [1985]
m.CH4_ace.Kd      =  0.01; %[d^-1] - TIM: Changed this because I couldn't justify it looking over the equations in Xu which made it so that the cells could only die off, never grow
m.CH4_ace.Y       =  0.040; %[mol C (mol-ace-C)^-1] - Growth efficiency of acetoclastic methanogens - Kettunen [2003]
m.CH4_ace.CH4.R   =  0.50;  %[mol-CH4/(mol-ace)^-1] - Rate of CH4 production - Kettunen [2003] - NOTE: This is not Monod kinetics, so may be wiser to ignore
% m.CH4_ace.ace.Km_PDF  =  4-700; %[mmol m^-3] Half-saturation coefficient - [Kettunen 2003]

%TIM MADE UP VALUES THAT MUST BE LOOKED UP - Based these off of mu_max. Might be better to get Vmax and determine mu from that, which is what I've done in the actual code
m.CH4_ace.ace.Vmax = m.CH4_ace.mu_max/(m.CH4_ace.Y*2); %FIXME - Tim made up - How fast they will consume acetate
% m.CH4_ace.CH4.Vmax = m.CH4_ace.ace.Vmax/2;             %FIXME - Tim made up - How fast they will produce CH4 - These should be based on stoichiometry
% m.CH4_ace.CO2.Vmax = m.CH4_ace.ace.Vmax/2;
% m.CH4_ace.H  .Vmax = m.CH4_ace.ace.Vmax;

% m.CH4_ace.Vmax_PDF = 0.03:0.05:2.88; %Range specified in Xu
% m.CH4_ace.Kd_PDF = 0.01:0.05:0.79; %Range specified in Xu
% m.CH4_ace.Y_PDF = 0.037:0.05:0.3; %Range specified in Xu
%% Low acetate methanogens (Jordan's bug)

%% Hydrogenotrophic methanogens
% Pulled from Robert Grant's 2002 paper
% m.CH4_H2 .H2 .Km  = 0.01;    %[gH m^-3] - Mosey [1983] and Robinson and Tiedje [1982]
% m.CH4_H2 .CO2.Km  = 0.12;    %[gC m^-3] - no reference given in Grant, but the value is here
% m.CH4_H2 .mu_max    = 0.12;    %[g C  g microbial C h^-1] - Shea et al. [1968], Zehnder and Wuhrmann [1977]

% Constants pulled from Xu et al 2015 JGR:B, Microbial functional group-based CH4 model, which are themselves pulled form other sources - Table 1
m.CH4_H2.mu_max   = 0.25;       %[d^-1] - Growth rate of H2-CO2-dependent methanogens - Servais et al [1985]
m.CH4_H2.Kd     = 0.01;       %[d^-1] - Death rate of H2-CO2-dependent methanogens - Servais et al [1985]
m.CH4_H2.Y      = 0.20;       %[mol-C (mol-CO2-C)^-1] - Growth efficiency of H2-CO2-dependent methanogens - Grant [1998]
m.CH4_H2.CO2.Vmax  = m.CH4_H2.mu_max/m.CH4_H2.Y; %[mol-CO2 d^-1] % Shouldn't this be noramlized to microbe?
m.CH4_H2.H2 .Vmax  = m.CH4_H2.CO2.Vmax*4;        %[mol-H2  d^-1]


% m.CH4_H2.H2.Km  = 7.75*10^-6; %[mmol m^-3] - Half coefficient of H2 for CH4 production from H2 - Fennell and Gossett [1998]
% m.CH4_H2.CO2.Km = 3.10*10^-8; %[mmol m^-3] - Half coefficient of CO2 for CH4 production from H2 - Stoichiometry, no other ref given
m.CH4_H2.H2.Km  = .775; %[mmol m^-3] - Half coefficient of H2 for CH4 production from H2 - Fennell and Gossett [1998]
m.CH4_H2.CO2.Km = .31; %[mmol m^-3] - Half coefficient of CO2 for CH4 production from H2 - Stoichiometry, no other ref given


% m.CH4_H2.Vmax_PDF = 0.005:0.1:0.041; %Range specified in Xu. Again, the value given is not in this range, which is odd
% m.CH4_H2.Kd_PDF = 0.01:0.005:0.031; %Range specified in Xu

%% Orphaned constants and such - Things I don't yet know what to do with
% m.CH4_ace.KPO4  = 0.125;    %????
% m.CH4_ace.KPNO3 = 0.350;    %????
% m.CH4_ace.KPNH4 = 0.400;    %????

m.CH4_prod.ace.Km = 5; %[mmol m^-3] - Half-saturation coefficient (TIM: not clear on this one) - Kettunen [2003]
% m.CH4_prod.ace.Km = 4:10:700; %Range specified in Xu

m.CH4_prod.O2.Km = 0.01; %[mmol m^-3] - Half-saturation coefficient (TIM: I am unclear on what organism this is meant to go to. Can't figure it out) - Kettunen [2003]
% m.CH4_prod.O2.Km_PDF = 0.0002:0.0005:0.040 ; %Range specified in Xu


% pH=0:0.1:14;
% plot(pH,Gf_H_ion);
% dG = -n*F*dE0;







% % %     s.H.Gf = -5.69*pH;
% % %     m.CH4_ace.dG0_cat = 0.5*s.CO2.Gf + 0.5*s.CH4.Gf - s.H.Gf - s.ace.Gf;
% % %     m.CH4_ace.dG0_ana = 
% % % %     s.CH4.Gf =  -50.75;   %kJ/mol - Brock - Biology of Microorganisms, 2012 13th ed, appx 1
% % % %     s.CO2.Gf = -394.4;    %kJ/mol
% % % %     s.ace.Gf = -369.41;   %kJ/mol







%% The below is from Liu roughly. Good reference to start to try to model based on Gibb's. First pass will not take that approach
% log10(gamma) = -A*z^2*mu^0.5/((1+mu^0.5)-0.3*mu)
% dGa = dGc - dGm - dGw => anabolism = catabolism - maintenance - waste
% % % Y = c.dGcat/(c.dGw - c.dGa_avail - c.dGm) ;%Observed yield - based on Liu et al 2007 eqn 7 adding a term for maintenance
% % % dGa/ %Conversion from J to mol, or mol/J conversion factor. What term would fit that description?
% % % dX=dGa/dcGx; %Maybe something like this?
% Then yield would be:
% % % Y = -dX/dS; %where dS is provided by the maximum uptake rate
%So you MUST get the Gibb's terms for each reaction, the uptake rate for
%each reaction, sensitivity to abundance of elements, etc



