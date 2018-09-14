%-------------------------------------------------------------------------------
%---------------------DESCRIPTION OF VARIABLES----------------------------------
%c=constants/parameters
%  c.dt = Model time step                                                       (s)
%  c.nCol = number of columns                                                   (--)
%  c.nRow = number of rows                                                      (--)
%  c.nLyr = number of soil layers                                               (--)
%  c.R    = Universal gas constatn                                              FIXME

%p=pool of a compound
%  p.NO3 = Nitrate
%  p.NO2 = Nitrite
%  p.NH4 = Ammonia
%  p.CH4 = Methane
%  p.CO2 = Carbon dioxide
%  p.but = Butyrate
%  p.mon = Generic, undefined monomers
%  p.ace = Acetate
%  p.H2  = Elemental hydrogen
%  p.O2  = Elemental oxygen
%  p.PO3 = Phosphate
%  p.PO2 = Phosphite
%  p.PO4 = ???
%  p.eth = Ethanol
%  p.meth= Methanol
%  p.mic.C = Dead microbial Carbon
%  p.mic.N = Dead microbial Nitrogen
%  p.mic.P = Dead microbial Phosphorus
%    p.XXX.d      - dissolved                                                   (mol L^-1)
%    p.XXX.g      - gaseous                                                     (mol L^-1)
%    p.XXX.s      - sorbed                                                      (mol L^-1)
%    p.XXX.p      - particulate/solid/suspended                                 (mol L^-1)
%    p.XXX.kH     - Henry's constant for compound XXX                           (  gas mol L^-1 / water mol L^-1)
%    p.XXX.kOC    - Sorption constant to organic carbon                         ( soil mol L^-1 / water mol L^-1)
%    p.XXX.kPa    - Particulate/dissolved constant                              (part. mol L^-1 / water mol L^-1)
%    p.XXX.sat    - Saturated concentrtion in water at T=20C, pH=7              (mol L^-1)
%    p.XXX.MW     - Molecular weight                                            (g mol^-1)
%    p.XXX.flux.d - dissolved flux
%    p.XXX.flux.g - gaseous flux
%    p.XXX.flux.s - sorbed flux
%    p.XXX.flux.p - particulate flux

%  p.lig = Lignin
%  p.pec = Pectin
%  p.carb= Generic, undefined carbohydrates
%    p.XXX.C - Carbon for the complex molecules above                           (gC L^-1)
%    p.XXX.N - Nitrogen for the complex molecules above                         (gN L^-1)
%    p.XXX.P - Phosphorus for the complex molecules above                       (gP L^-1)

%m=microbial population - Vectorize or matricize these to make the fully distributed trait model
%  m.XXX                                                                        (g L^-1)
%  m.Het = Aerobic heterotophs
%  m.Pla = Plants
%  m.Fun = Fungi
%  m.Fer = Fermenters
%  m.Ace = Acetogens
%  m.Syn = Syntrophs
%  m.CH4_A   = Acetoclastic methanogens
%  m.CH4_H   = Hydrogenotrophic methanogens
%  m.CH4_M   = Methyolotrophic methanogens
%  m.CH4_TA  = Aerobic methanotrophs
%  m.CH4_AOM = Anaerobic methanotrophy - may need to be expanded
%  m.Fe      = Iron reducers
%  m.S       = Sulfate reducers
%  m.Mn      = Manganese reduers
%  m.NitF    = Nitrogen fixers
%  m.DeNO3   = Denitrifiers
%  m.DeNO2   = Nitrosifiers
%  m.AOM     = Ammonia oxidizing bacteria
%  m.NOB     = Nitrogen oxidizing bacteria
%  m.Anmox   = Annamox
%    m.XXX.tot        = Carbon in biomass                                       (gC L^-1)
%    m.XXX.C          = Carbon in biomass                                       (gC L^-1)
%    m.XXX.N          = Nitrogen in biomass                                     (gN L^-1)
%    m.XXX.P          = Phosphorus in biomass                                   (gP L^-1)
%    m.XXX.QC         = BC/BT - biomass C to total biomass ratio                (--)
%    m.XXX.QN         = BN/BT - biomass N to total biomass ratio                (--)
%    m.XXX.QP         = BP/BT - biomass P to total biomass ratio                (--)
%    m.XXX.YYY.Km     = Half saturation constant for substrate YYY              (umol L^-1)
%    m.XXX.YYY.Ki     = Inhibition constant for substrate YYY                   (umol L^-1)
%    m.XXX.YYY.Vmax   = Maximum substrate uptake rate for substrate YY          (M s^-1)
%    m.XXX.YYY.Yn     = Substrate use efficiency constant for YYY               (--)
%    m.XXX.YYY.dem    = Microbial demand for substrate YYY
%    m.XXX.mu_max     = Maximum specific growth rate                            (day^-1)

%    m.XXX.flux.g     = Growth of microbial population                          (g L^-1)
%    m.XXX.flux.d     = Death of microbial population                           (g L^-1)
%    m.XXX.flux.m     = Maint. of microbial population                          (g L^-1)

%s=soil properties
%  s.por = porosity                                                             (--)
%  s.Ksat= saturated hydraulic conductivity                                     (m s^-1)

%w=water properties
%  w.pH   = pH of the pore water                                                (--)
%  w.thet = water potential                                                     (Pa)
%  w.sat  = saturated water ___                                                 (FIXME)
%  w.temp = Temperature                                                         (K)

%--------------------END DESCRIPTION OF VARIABLES-------------------------------
%-------------------------------------------------------------------------------