%Written by Steve Colquhoun 09/14/2018
%Fall 2018 Research with Dr. Tim Morin
%SUNY College of Environmental Science and Forestry

%Function to partition chemical concentrations between liquid, solid, and
%gaseous phases

%Inputs: conc_aq = aquesous concentration     Units: Mass/Length^3
%        Kh      = Henry's Law Constant       Units: ---
%        Kf      = Freundlich Constant        Units: ---
%        n       = Freundlich Exponent        Units: N/A
%        phase   = bianry to indicate desired phase concentration
            %% phase = 0 --> solid
            %% phase = 1 --> gas
%Outputs: conc_D   = desired concentration      Units: g/m^3 

function [conc_D] = partitioning(conc_aq, Kh, n, Kf, phase)

if phase == 0
    %Solid phase desired, apply Freundlich Isotherm
    conc_D = Kf*conc_aq^n;
else
    %Gas phase desired, apply Henry's Law
    conc_D = Kh*conc_aq;
end

end
