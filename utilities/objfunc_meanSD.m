function [residual] = objfunc_meanSD(param,Cparam,setup,x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the difference between predicted and observed mean sleep duration.
% Designed to work with fzero - so the sign is important.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if setup.fit_param == 1 | setup.fit_param == 2
    param.mu     = x0(1)
  end
  xc           = setup.Y0(2);
  yc           = setup.Y0(3);
  H0           = param.H0p+circadian(xc,yc,param);
  setup.Y0     = [H0 setup.Y0(2:4)];
%
  [Y0,convergence_measure] = remove_transient(param,Cparam,setup);
  setup.Y0    = Y0(1:4);
  param.onoff = Y0(5);
%
  [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup);
%
  sleep_on_days  = sleep_on/(24*60*60);
  sleep_off_days = sleep_off/(24*60*60);
  sleepparam = find_sleep_param_from_on_off_times(sleep_on_days,sleep_off_days,setup);
%
  residual =  sleepparam.mean_sleep_duration*24-setup.diary_SD;
%
end
