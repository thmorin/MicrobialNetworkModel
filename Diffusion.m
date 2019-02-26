function [dCdt]=Diffusion(C,dz,D)
% Written by Tim Morin
% 2/15/2019
% Evaluate diffusion according to Fick's 2nd Law of Diffusion
% C should be 3 elements, dz sould be 2
% D is diffusion constant of the solute in the solvent

dCdt=D*(C(1)-C(2))/dz(1) - (C(2)-C(3))/dz(2)/mean(dz);