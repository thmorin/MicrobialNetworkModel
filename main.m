% % % % function []=gibbsFreeEnergyShit(fucker,fucker2)
cd C:/Users/thmorin/Desktop/myModel/
clear;close;clc
%Written by Tim Morin, 06/2018
%
%
%This code will use a simple Euler function to dynamically represent
%a micrboial population's growth over time, limiting for substrate, holding capcity,
%base on bioenergetics described in Grant's papers, Bouskill's paper, and Brock's microbiology

%The comments below lay out the general structure of the model, defining pools
%and assigning important variable names that need to be tracked explicitly


%-------------------------------------------------------------------------------
%--------------------------INITAL CONDITIONS------------------------------------
%             Initial conditions - microbial populations (mol L^-1)
m.Het.tot=100;      m.Pla.tot=100;    m.Fun.tot=100;    m.Fer.tot=100;    m.Ace.tot=100;
m.Syn.tot=100;
m.CH4_A.tot=100;    m.CH4_H.tot=100;  m.CH4_M.tot=100;  m.CH4_TA.tot=100; m.CH4_AOM.tot=100;
m.names=fieldnames(m);
for i=1:length(m.names)

%Distribute initial weight then distribute according to the general distribution
%of biomass described by: C5 H7 O2 N P(1/12)
    eval(['m.' m.names{i} '.C=m.' m.names{i} '.tot*0.85']);      eval(['m.' m.names{i} '.N=m.' m.names{i} '.tot*0.10']);      eval(['m.' m.names{i} '.P=m.' m.names{i} '.tot*0.05']);
end
%                    Initial soil properties
s.por=0.3;     s.Ksat=100;

%                Initial compound concentrations (mol L^-1)
p.NO3.d=100;  p.NO2.d=100;  p.NH4.d=100;  p.NH3.d=100;  p.NH4OH.d=100; 
p.CH4.d=100;  p.CO2.d=100;  p.but.d=100;  p.mon.d=100;  p.eth.d=100;   p.meth.d=100;
p.ace.d=100;  p.H2.d =100;  p.O2.d =100;
p.PO4.d=100;  p.PO3.d=100;  p.PO2.d=100;
p.mic.C.d=100;  p.mic.N.d=100;  p.mic.P.d=100;

p.lig.C.p=100;  p.lig.N.p=100;  p.lig.P.p=100;
p.pec.C.p=100;  p.pec.N.p=100;  p.pec.P.p=100;
p.carb.C.p=100; p.carb.N.p=100; p.carb.P.p=100;

%                    Intial water conditions
w.pH=7; w.temp=273.15+20;   w.thet=1;   w.sat=1;
%------------------------END INITIAL CONDITIONS---------------------------------
%-------------------------------------------------------------------------------


%-------------------------------------------------------------------------------
%--------------------------MODEL  MAIN ROUTINE----------------------------------
stopCond=0;
while stopCond~=1
  %-----------------------------------------------------------------------------
  %                 Refigure states of pools
  for i=1:length(p.names)
    %               Adjust for temperature and pH
    %The functions below were structured to iterate through the compounds
    %I've written what the general form should look like in non-eval format for readability
  
    %            Adjust constants for temperature, pH
    %General form for line below looks like:
    %p.CH4.kH_adj =  kHAdjust(p.CH4.kH ,w.temp,w.pH)
    eval(['p.' p.names{i} '.kH_adj =  kHAdjust(p.' p.names{i} '.kH ,w.temp,w.pH);']);
    
    %General form for line below looks like:
    %p.CH4.kH_adj = kOCAdjust(p.CH4.kOC,w.temp,w.pH)
    %FIXME: This should use Freundlich isotherm to come up with equilibrium
    eval(['p.' p.names{i} '.kH_adj = kOCAdjust(p.' p.names{i} '.kOC,w.temp,w.pH);']);
    
    %General form for line below looks like:
    %p.CH4.kH_adj =  kPAdjust(p.CH4.kP ,w.temp,w.pH)
    eval(['p.' p.names{i} '.kH_adj =  kPAdjust(p.' p.names{i} '.kP ,w.temp,w.pH);']);
    
    %General form for lien below looks like:
    %p.CH4.sat_adj= satAdjust(p.CH4.sat,w.temp,w.pH)
    eval(['p.' p.names{i} '.sat_adj= satAdjust(p.' p.names{i} '.sat,w.temp,w.pH);']);
    %---------------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------------
    %             Calculate equilibrium states  
    %General form for line below looks like:
    %[p.CH4.d, p.CH4.g]=HenrysLaw   (p.CH4.d, p.CH4.g, c.kH_adj. CH4);
    eval(['[p.' p.names{i} '.d, p.' p.names{i} '.g]=HenrysLaw   (p.' p.names{i} ...
       '.d, p.' p.names{i} '.g, c.kH_adj. ' p.names{i} ');']);   %Equilibrate gas/dissolved according to Henry's Law
    
    %General form of line below looks like:
    %[p.CH4.d, p.CH4.s]=OctanolWater(p.CH4.d, p.CH4.s, c.kOC_adj.CH4);
    eval(['[p.' p.names{i} '.d, p.' p.names{i} '.s]=OctanolWater(p.' p.names{i} ...
       '.d, p.' p.names{i} '.s, c.kOC_adj.' p.names{i} ');']);   %Equilibrate sorbed/dissolved according to Octanol Water coefficient
    
    %General form of line below looks like:
    %[p.CH4.d, p.CH4.p]=PartDiss    (p.CH4.d, p.CH4.p, c.kP_adj. CH4);
    eval(['[p.' p.names{i} '.d, p.' p.names{i} '.p]=PartDiss    (p.' p.names{i} ...
       '.d, p.' p.names{i} '.p, c.kP_adj. ' p.names{i} ');']);   %Equilibrate particulate/dissovled FIXME - not sure what this is
    %---------------------------------------------------------------------------
    
    
    %---------------------------------------------------------------------------
    %               Remove supersatured dissolved gaseous species    
    %General form of the following if statement
    % if p.CH4.g > p.CH4.sat_adj
    %   p.CH4.flux.g = p.CH4.g - p.CH4.sat_adj;
    %   p.CH4.g      = p.CH4.sat_adj;
    % end
    if eval(['p.' p.names{i} '.g > p.' p.names{i} '.sat_adj;'])
      eval( ['p.' p.names{i} '.flux.g = p.' p.names{i} '.g - p.' p.names{i} '.sat_adj;']);
      eval( ['p.' p.names{i} '.g      = p.' p.names{i} '.sat_adj;']);
    end
    %---------------------------------------------------------------------------
  end
  
  %Al(OH)3(s)     <->   (Al 3+) + 3*(OH -)
  %Fe(OH)3(s)     <->   (Fe 3+) + 3*(OH -)
  
  %CaSO4(s)       <->   (Ca 2+) +   (SO4 2-)
  %AlPO4(s)       <->   (Al 3+) +   (PO4 3-)
  %FePO4(s)       <->   (Fe 3+) +   (PO4 3-)
  
  %CaCO3(s)       <->   (Ca 2+) +   (CO3 2-)
  
  %Ca(H2PO4)2(s ) <->   (Ca 2+) + 2*(H2PO4-)
  %Ca(HPO4)(s)    <->   (Ca 2+) +   (HPO4 2-)
  %Ca5(PO4)3OH(s) <-> (5(Ca 2+) + 3*(PO4 3-)  + (OH -)
  
  %X2-Ca + 2(NH4 +) <-> 2*X-NH4 +  (Ca 2+)
  %3X-Ca + 2*(Al 3+) <-> 2*X-Al + 3(Ca 2+)
  
  %-----------------------------------------------------------------------------
  %             What is the demand for each substrate
  %Catabolic - Providing energy for reaction
  % Glucose + O2 <-> CO2 + H2O
  
  %Anabolic - Adding to biomass
  % Carb + O2 <-> Mic
  MM_CH4_A(m.CH4_A,p.ace)
  %-----------------------------------------------------------------------------
  
  
  %% Acetate
  %What is S_1_T?
  %What is E_j_T?
  
  
  % Find who uses the substrate
  % Calculate E_T for that substrate   - eqn 17 Tang and Riley
  %Then go through and do the summation of Ek_T/Ks_1k for that substrate - eqn 17 Tang and Riley
  %Unless I want to do inhibition? - eqn 18 Tang and Riley
  
  
  
  S_1_T = S_1;
  for k=1:how_many_pops_use_substrate
    S_1_T = S_1_T + C_1_k
  end
  
  E_1_T = E_1;
  for k=1:how_many_complexes_use_enzyme
    E_1_T = E_1_T + C_k_1
  end
  
  bot=1;
  for i=1:how_many_pops use_substrate
    bot = bot + E_k_T/Ks_1k;
  end
  
  C_1_j=(alpha_j*S_1_T*E_j_T) / (Ks_1j*bot + S_1_T);                            %Formulation from Tang and Riley 2013, Biogeosciences
  dP_i_j_dt = C_i_j * k_i_j_2;
  dS_i_dt   = C_i_k * k_i_k_2;
  
  
  % Autotrophic ammonia oxidizing bacteria/archaea - formulations from Bouskill 2012
  m.AOB.QBmin = min([m.AOB.QC m.AOB.QN m.AOB.QP])
  m.AOB.dBC   = max(1 - m.AOB.QBmin/m.AOB.QC);                                  %Bouskill 2012 eqn 1
  m.AOB.dBN   = max(1 - m.AOB.QBmin/m.AOB.QN);                                  %Bouskill 2012 eqn 1
  m.AOB.dBP   = max(1 - m.AOB.QBmin/m.AOB.QP);                                  %Bouskill 2012 eqn 1
  
  m.AOB.dBmin = min([m.AOB.dBC m.AOB.dBN m.AOB.dBP]);
  m.AOB.DB    = m.AOB.mu_max * m.AOB.dBmin * m.AOB.tot;                         %Bouskill 2012 eqn 2

  m.AOB.NH3.V = ...
      m.AOB.NH3.Vmax * ...
      p.NH3.d/(m.AOM.NH3.Km + p.NH3.d*(1 + p.NH3.d/m.AOB.NH3.Ki)) *...
      p.O2.d/(m.AOM.O2.Km + p.O2.d) * ...
      m.AOM.tot;                                                                %Bouskill 2012 eqn. 3
  m.AOB.CO2.Vmax = ...
      m.AOB.NH3.Vmax * m.AOB.CO2.Yn / m.AOB.QN * ...  
      max(1- (m.AOB.rCN - m.AOB.rCN_min)/(m.AOB.rCN_max - m.AOB.rCN_min) , 0 ); %Bouskill 2012 eqn. 5
  m.AOB.CO2.V = ...
      m.AOB.CO2.Vmax * p.CO2.d / ( m.AOB.CO2.Km + p.CO2.d );                    %Bouskill 2012 eqn. 4
  
  m.AOB.dBdt = m.AOB.mu_max * m.AOB.dBmin * m.AOB.tot ...
%Could this be m.AOB.dBdt = m.AOB.DB for the first three terms?  
       - m.AOB.Kd * m.AOB.tot - 0.25*(m.???.DA_NO2 + m.???.DA_NO); %<--I don't know what these terms are
  
  
  % Nitrogen oxidizing bacteria
  m.NOB.NO2.V  = m.NOB.NO2.Vmax * p.NO2.d / ( m.NOB.NO2.Km + p.NO2.d ) *...
      p.O2.d / ( m.NOB.O2.Km + p.O2.d ) * m.NOB.tot;                            %Bouskill 2012 eqn 7
  
  m.NOB.dBC    = max(1 - m.NOB.QBmin/m.NOB.QC);                                 %Bouskill 2012 eqn 1
  m.NOB.dBN    = max(1 - m.NOB.QBmin/m.NOB.QN);                                 %Bouskill 2012 eqn 1
  m.NOB.dBP    = max(1 - m.NOB.QBmin/m.NOB.QP);                                 %Bouskill 2012 eqn 1
  m.NOB.dBmin  = min([m.NOB.dBC m.NOB.dBN m.NOB.dBP]);
  
  m.NOB.dBdt   = (m.NOB.mu_max*m.NOB.dBmin*m.AOB.tot) - (m.AOB.Kd*m.AOB.tot);   %Bouskill 2012 eqn 8
  
  % dBdt = mu_max*B - k_d*B - k_r*B where k_r is respiration
  
  m.NOB.NH3.Km_adj = m.NOB.NH3.Km    *(1 + m.NOB.NO2.Kd_max * p.NO2.d/p.O2.d);  %Bouskill 2012 eqn 9
  m.NOB.NH3.Km_adj = m.NOB.NH3.Km_adj*(1 + m.NOB.NO .Kd_max * p.NO .d/p.O2.d);  %Bouskill 2012 eqn 9
 
  p.NH3.dNH3dt     = -( m.AOM.NH3.Ve + m.AOM.NH3.Vb ) ...
      - m.NOB.NH3.V ...
      + 0.25 * ( m.???.DA_NO2 + m.???.DA_NO );                                  %Bouskill 2012 eqn 10
      
  p.NO2.dNO2dt     = m.AOB.NH3.Ve - m.AOB.NO2.Ve - m.???.DA_NO2;                %Bouskill 2012 eqn 11
 

  
  
  %Nitrification
  % NH4+ + 1.863*O2 + 0.098*CO2 -> 0.0196*C5H7NO2 + 0.98*NO3- + 0.0941*H2O + 1.98*H+   - Eqn 11.19 from Eco Eng textbook
  
  %Denitrification
  % C10H19O3N + 10*NO3- -> 5*N2(g) + 10*CO2 + 3*H2O + NH3 + 10OH-
  
  
  
  % E_j_T - can I sub microbial population in here?
  % 
  
  
  
  % k_i_k_2  - Forward reaction coefficient
  % K_s_1j   - substrate affinity coefficient
  % S_1_T    - Total substrate, free and binded
  % lowercase i - Substrate
  % lowercase j - Enzyme
  % lowercase k - Also enzyme, but used for the summation to differential b/w general j
  % C_1_j - total substrate binded with enzyme 1
  
  
  % Arrhenius function
  V_mod=V0*exp(delE/R*(1/T0 - 1/T));
  
  % pH sensitivity of parameters and biomass
  f=(pH-pH_min)*(pH-pH_max)/( (pH-pH_min)*(pH-pH_max) - (pH-pH_opt)^2); %Bell curve from CLM-Microbe, MEM, TEM
  f=10^(-0.2335*pH^2 + 2.7727*pH - 8.6); %Bell curve from CLM4Me
  
  if pH>4 && pH<7;          f=1.02/(1 + 1000000*exp(-2.5*pH));       %Bell curve from DLEM
  elseif pH>=7 && pH<10;    f=1.02/(1 + 1000000*exp(-2.5*(14-pH)));  end
  
  % Moisture effects
  f=(M_V-M_min)*(M-M_max)/( (M_V-M_min)*(M_V-M_max) - (M_V-M_opt)^2 );  %Bell curve from TEM
  
  %-----------------------------------------------------------------------------
  %             Gibbs Free Energy
  dG=dGf + R*T*log(P1*P2*P3/R1/R2/R3);  %Activity coefficients?
  
  
  %-----------------------------------------------------------------------------
  %     Modify K values for spheric cell - aerobic heterotrohps and shit
  Ks_mod_sphere = Ks*(1 + Vmax/(4*pi*D*r*n*Ks) );
  
  %-----------------------------------------------------------------------------
  %     Modify K values inhibition
  bot=0;
  for i=1:totInhib;      bot=bot+I/KI;       end
  Ks_mod_inhib = Ks/(bot+1);
  %-----------------------------------------------------------------------------
  %           Compare demand to available
  % Maybe not needed with ECA?  
  
  
  
  
  %-----------------------------------------------------------------------------
  %           ECA
  bot1=0;       bot2=0;
  for i=1:subsThatMicWantsToEat;    bot1=bot1+S/Ks_mod_inhib;        end
  for j=1:micsThatWantIt;           bot2=bot2+M/Ks_mod_inhib;        end
  
  complex = S*M / (Ks * (1 + bot1 + bot2));
  dSdt = V*complex
  
  
  
  
  
  
  
  
  %Generalized working - probably trash below
  p.O2.ET=0; %Initialize total enzyme population
  p.O2.ST=p.O2.d / p.O2.MW; %Convert it to mol L^-1
  for i=1:length(p.names)
    for j=1:length(m.names)
      p.O2.ET=p.O2.ET + m.O2.Het.tot;  %mol L^-1  --  only add microbial group that use these substrates
    end
    p.O2.
    for j=1:length(m.names)
      if isfield(p.names{i},eval(['m.' m.names{j}]))
        if isfield('O2',m.Het)
        p.O2.ST / m.Het.O2.Km
        
      end
      p.O2.Het.=   % mol m^-3
      p.O2 =    %mol m^-3
    end
  end
  m.Het.O2.dem =p.O2. ST*p.O2. ET/(m.Het.O2. Km*(1+p.O2. ST/m.Het.O2. Km+p.O2. ST/m.Ace.O2 .Km)+p.O2. ST); %From eqn. 13 of Tang and Riley
  m.Het.glu.dem=p.glu.ST*p.glu.ET/(m.Het.glu.Km*(1+p.glu.ST/m.Het.glu.Km+p.glu.ST/m.Ace.glu.Km)+p.glu.ST);
  m.Het.NH4.dem=p.NH4.ST*p.NH4.ET/(m.Het.NH4.Km*(1+p.NH4.ST/m.Het.NH4.Km+p.NH4.ST/m.Ace.NH4.Km)+p.NH4.ST);
  
  
  p.O2.ST/m.Het.O2.Km
  
  %-----------------------------------------------------------------------------
  
  
  %-----------------------------------------------------------------------------
  %       Are microbes getting what they need for maintenance
  
  %-----------------------------------------------------------------------------
  
  
  %-----------------------------------------------------------------------------
  %       Calculate fluxes associated with microbeial biomass
  
  %-----------------------------------------------------------------------------
  
  
  %-----------------------------------------------------------------------------
  %          Apply fluxes to all pools and populations
  
  %-----------------------------------------------------------------------------
  
end
%-----------------------END MODEL  MAIN ROUTINE---------------------------------
%-------------------------------------------------------------------------------