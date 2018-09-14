function [fm]=f_h2o(K_thet,thet_op,thet,a,n_s,phi,b)
% TIM NOTE: I read this paper and coded this function but I don't think
% it's ovely useful for the type of model I'm developing. I don't represent
% the separation of carbon stocks from microbes, but I will represent the
% sorption/desorption of it to the soil and I'll definitely represent the
% O2 availability.
%
% Used to describe the sensitivity to soil moisture of microbial processes
% Based on Yan et al 2018 - Nature Communications - A moisture function of soil heterotrophic respiration that incorporates microscale processes
% Concept of model is that at low moisture contents microbes can have
% trouble accessing organic carbon while high moisture contents will limit
% the oxygen available to them.
%
% Variables:
% fm      - the adjustment formula for heterotrophic respiration (HR)
% K_thet  - moisture constant relfecting imapct of water content on soil-adsorbed organic carbon desporption
% thet    - water content
% thet_op - optimal water content at which HR rate peaks
% a       - integrating parameter, 0 to 1 - describes the separation of microbes from their substrate as a function of soil moisture
% b       - describes oxygen sensitivity. if used at all, 0.17 seems like a good default value.

if thet<thet_op
    fm = (K_thet + thet_op)/(K_thet + thet)*(thet/thet_op)^(a*n_s);
elseif thet>thet_op
    fm = (phi - thet)/(phi-thet_op)^b; %Describes oxygen availability. Can probably be turned off if you are going to explicitly track oxygen pool
end