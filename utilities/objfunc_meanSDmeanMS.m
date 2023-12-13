function [residual] = objfunc_meanSDmeanMS(x0,param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates the difference between predicted and observed mean sleep onset and offset times.
% The residual is given as 
%    residual = sqrt(residual_SD.^2 + residual_MS.^2) 
% where
% residual_SD is the residual from the sleep duration and 
% residual_MS is the residual from the mid-sleep.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  if setup.fit_param == 1
    param.mu     = x0(1)
    Cparam.tau_c = x0(2)
  end
  if setup.fit_param == 2
    param.mu    = x0(1)
    param.a     = x0(2)
  end
%
  [Y0,convergence_measure]       = remove_transient(param,Cparam,setup)
  setup.Y0                       = Y0(1:4);
  param.onoff                    = Y0(5);
  [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup)
%
  n_on           = length(sleep_on);
  n_off          = length(sleep_off);
%
  sleep_on_days  = sleep_on/(24*60*60);
  sleep_off_days = sleep_off/(24*60*60);
  [sleepparam] = find_sleep_param_from_on_off_times(sleep_on_days,sleep_off_days,setup)
%
  residual_SD = sqrt((sleepparam.mean_sleep_duration*24-setup.diary_SD).^2)
  residual_MS = sqrt((sleepparam.mean_midsleep*24-setup.diary_MS).^2) 
  residual = sqrt(residual_SD.^2 + residual_MS.^2) 
%
  if setup.plt_raster == 1
    figure(setup.ifig)
    plt_raster_line(sleepparam.sleep_onsets,sleep_param.sleep_offsets,0)
  end
%
end
