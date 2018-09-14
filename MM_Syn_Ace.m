function []=MM_Syn_Ace()
% Taken from Brock's Microbiology of Microorganisms, version 13, pp.381
% Holding capacity
% Max consumption rate
% Km
% Yield

%Gibb's Free Energy of Formation                                                kJ
Gf = +48.2;                          %Butyrate- + 2*H2O -> 2*acetate- + H+ + 2H2

K = (H*H2^2*acetate^2)/(Butyrate)

G=Gf + R*T*ln(K)

% Can the energetics sustain the maintenance?

demand=population*
endfunction
