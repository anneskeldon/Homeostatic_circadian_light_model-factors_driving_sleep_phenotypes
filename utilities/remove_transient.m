function [Y0,convergence_measure] = remove_transient(param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeatedly integrates for the period from setup.T0 to setup.Tend until
% either there is little change between the predicted sleep onset and
% offset times or it has looped through setup.nloop times.
% Note, it implicitly assumes that setup.T0 to setup.Tend is an 
% integer number of days. 
%
% Returns the initial values for the integrator and the convergence_measure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  iloop               = 0;
  sleep_on            = [];
  convergence_measure = 1000; % seconds

  while convergence_measure > 10 & iloop < setup.nloops
    [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup);
    setup.Y0    = Y0;
    param.onoff = onoff;
    n_on        = length(sleep_on);
    n_off       = length(sleep_off);
    if iloop == 0
      sleep_on_old = sleep_on;
      n_on_old     = n_on;
    else
      nsleeps = min(n_on,n_on_old);
      convergence_measure = sqrt(sum( (sleep_on_old(1:nsleeps) - sleep_on(1:nsleeps) ).^2 ));
      sleep_on_old = sleep_on;
      n_on_old     = n_on;
    end
    iloop = iloop + 1;
  end
%
  Y0(5) = param.onoff;
%
end
