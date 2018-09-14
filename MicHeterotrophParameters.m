m.Het.O2 .Y  = 1;

%-------------------------------------------------------------------------------
%                     Half saturatino constants
m.Het.O2   .Km = 1; 
m.Het.NH4  .Km = 1;
m.Het.NH3  .Km = 1;
m.Het.NH4OH.Km = 1;
m.Het.glu  .Km = 1;
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
%                       Inhibition constants
m.Het.O2   .Ki = inf; %No inhibition? - no, doesn't work like this. Just don't include if no inhibition
m.Het.NO2  .Ki = inf; %FIXME
m.Het.NH4  .Ki = inf; %FIXME
m.Het.NH3  .Ki = inf; %FIXME
m.Het.NH4OH.Ki = inf; %FIXME
m.Het.NO3  .Ki = inf; %FIXME
m.Het.NO2  .Ki = inf; %FIXME
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
%                       Inhibition constants
m.Het.O2   .Y = 1;
m.Het.NH4  .Y = 1;
m.Het.NH3  .Y = 1;
m.Het.NH4OH.Y = 1;
m.Het.glu  .Y = 1; %FIXME
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
%             Prescribed microbial stoichiometry ratios
m.Het.rcnMax = 13.2; %From Bouskill et al 2012
m.Het.rcnMin =  6.6; %From Bouskill et al 2012
m.Het.rcpMax = 1;    %FIXME
m.Het.rcpMin = 1;    %FIXME
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
%                        Biomass quotas
m.Het.Qc     = 0.85; %FIXME BC/Bt - biomass c to biomass t
m.Het.Qn     = 0.10; %FIXME
m.Het.Qp     = 0.05; %FIXME



m.Het.V_max = 1; %FIXME

m.Het.