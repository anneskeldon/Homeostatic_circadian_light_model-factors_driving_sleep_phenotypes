function [Cparam,param] = default_HCL_parameters(age)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default sleep and circadian parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%
% Circadian model parameters (Forger)
%
  Cparam.alpha0  = 0.16; 
  Cparam.I0      = 9500;
  Cparam.p       = 0.6; 
  Cparam.G       = 19.9; 
  Cparam.b       = 0.4; 
  Cparam.gam     = 0.23;
  Cparam.kappa   = (12/pi)*60*60; 
  Cparam.f       = 0.99669; 
  Cparam.k       = 0.55;
  Cparam.lambda  = 60/(60*60); 
  Cparam.betaval = 0.013;
  Cparam.tau_c   = 24.2;  
%
% Sleep parameters
%
  param.Xs   = 45*60*60; 
  param.Xw   = 45*60*60;
  param.H0p  = 13.5;
  param.H0m  = 12.5;
%
% Age 30
  if age == 30
    param.a    = 2.72; 
    param.mu   = 18.67; 
  end
  if age == 65
    param.a    = 1.72; 
    param.mu   = 17.87; 
  end
  param.onoff = 0;  %If onoff=0 then in sleep state.  IF onoff=1 then in wake state.
%
% Light parameters for synthetic light profile
%
  Cparam.daylight     = 700; 
  Cparam.eveninglight = 40;
  Cparam.dawn         = 7.5*60*60; 
  Cparam.dusk         = 16.5*60*60;
%
  Cparam.gated = 1; % 1 means that light is self-selected. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
