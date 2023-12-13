function [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrates the HCL equations from setup.T0 to setup.Tend with initial conditions given in setup.Y0
% and returns the sleep onset and offset times and
% parameters that specify the values of the dependent variables at the end of the
% time period.
%
% For SPONTANEOUS sleeping i.e. assumes transitions between sleep and wake occur at the thresholds.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sleep_on = [];
  sleep_off = [];
  T0        = setup.T0;
  Tend      = setup.Tend;
  Y0        = setup.Y0;
%
% Make sure that the initial conditions are such that if on wake, its below the upper threshold,
% and if on sleep it is above the lower threshold.
%
  if param.onoff == 1 &  (Y0(1) > param.H0p + circadian(Y0(2),Y0(3),param))
    param.onoff = 0;
  end
  if param.onoff == 0 &  (Y0(1) < param.H0m + circadian(Y0(2),Y0(3),param))
      param.onoff = 1;
  end
%
  while T0<Tend
     [T,Y,te,ye,ie] = integrate_equations(Cparam,param,setup,T0,Tend,Y0);
    if ie == 1 % go to sleep event
      sleep_on = [sleep_on,te(end)];
    end
    if ie == 2  % wake up event
      sleep_off = [sleep_off,te(end)];
    end
%
% Plot for checking
%   figure(1)
%   plot(T,Y(:,1),'k-')
%   hold on
%   if ie == 1
%     plot(te,ye(1),'ro')
%   end
%   if ie == 2
%     plot(te,ye(1),'bo')
%   end
% 
% Set   up initial conditions for next period
%
    H0 = Y(end,1);
    x0 = Y(end,2);
    y0 = Y(end,3);
    n0 = Y(end,4);
    Y0 = [H0 x0 y0 n0];
    T0 = T(end,1);
%
    if T0 < Tend
      param.onoff=abs(param.onoff-1);
    end
    onoff = param.onoff;
  end  % End while loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
