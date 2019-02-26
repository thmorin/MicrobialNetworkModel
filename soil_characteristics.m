%   Written by Steve Colquhoun 10/06/2018
%   Fall 2018 Research for Dr. Tim Morin
%   SUNY College of Environmental Science and Forestry

%   Function for calculating soil characteristics

%   Inputs: rho_bulk  = bulk density, Ms/V_tot         Units: Mass/Length^3
%           n         = porosity, Vv/V_tot             Untis: %
%           V_tot     = total volume                   Units: Length^3
%           M_water   = mass of water in sample        Units: Mass
%           rho_water = denisty of water in units bulk density

%   Outputs: rho_samp = wet sample desnity             Units: Mass/Length^3
%            S        = degree of saturation           Units: %
%            e        = void ratio                     Units: --
%            Va       = volume of air                  Units: Length^3
%            Vw       = volume of water                Units: Length^3
%            Vs       = volume of soil                 Units: Length^3
%            w        = water content                  Units: %
%            Ms       = mass of dry soil               Units: Mass


function[rho_samp,S,e,V_air,V_soil,V_water,w,M_soil] = soil_characteristics(rho_bulk,n,V_tot,M_water,rho_water)
                                      
V_water = M_water/rho_water;

%   Volume of voids, V_voids = V_air+V_water           Units: Length^3
V_voids = n*V_tot;
V_soil = V_tot-V_voids;

V_air = V_voids-V_water;
M_soil = V_tot*rho_bulk;
w = (M_water/M_soil)*100;
rho_samp = (M_soil+M_water)/V_tot;
e = V_voids/V_soil;
S = (V_water/V_voids)*100;

end
