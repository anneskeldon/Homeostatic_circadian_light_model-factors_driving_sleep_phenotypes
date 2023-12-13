function F = objfunc_fsolve_meanSDmeanMS(x0,param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulation of the objective function for using fsolve to find parameters that
% match mean and observed mid-sleep and mean and observed sleep duration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if setup.fit_param == 1
    param.mu     = x0(1);
    Cparam.tau_c = x0(2);
  end
  if setup.fit_param == 2
    param.mu    = x0(1);
    param.a     = x0(2);
  end
%
  [Y0,convergence_measure] = remove_transient(param,Cparam,setup);
  setup.Y0    = Y0(1:4);
  param.onoff = Y0(5);
  [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup);
%
  if ~isempty(sleep_on) & ~isempty(sleep_off)
    sleep_on_days  = sleep_on/(24*60*60);
    sleep_off_days = sleep_off/(24*60*60);
    sleepparam = find_sleep_param_from_on_off_times(sleep_on_days,sleep_off_days,setup);
    if setup.plt_raster == 1
      figure(setup.ifig)
      plt_raster_line(sleepparam.sleep_onsets,sleepparam.sleep_offsets,0)
      xlim([-1 1])
      datetick('x')
    end
  else
    mean_SD = 2;
    mean_MS = 1;
  end
%  
  F(1) = sleepparam.mean_midsleep*24-setup.diary_MS;
  F(2) = sleepparam.mean_sleep_duration*24-setup.diary_SD;
%
end
