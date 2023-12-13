%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes a set of parameters, runs the HCL model in spontaneous mode, makes some plots and gives sleep onset, 
% offset, duration and midsleep and cmin. Here set up to loop over a grid of values varying tau_c 
% and the maximum level of the light (using a synthetic light profile). This generates data showing
% the entrainment regions (Arnold tongues).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  addpath('../utilities/')
  date = 20230923;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Parameter defaults
  [Cparam,param]    = default_HCL_parameters(30);
  [setup]           = default_setup_parameters();
  setup.maxstep     = 600;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
plt_check = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Need to make choices on light enviromnent and parameters
  Cparam.lightinput    = 0; % Means use the tanh profile in swlightinput.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  setup.nloops  = 40; % Maximum number of loops for convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  rundata_events = sprintf('temp/events.dat');
  setup.fid1 = fopen(rundata_events,'w');
  rundata_grid   = sprintf('temp/Arnold_tongue_%i.dat',date);
  setup.fid4     = fopen(rundata_grid,'w');
  save_header(param,Cparam,setup)
% 
  day_log   = linspace(1,4,129);
  lightvals = 10.^(day_log);
  tauvals   = 23.50:0.01:24.70;
  x0        = -0.81;
  y0        = -0.66;
  H0        = param.H0p+circadian(x0,y0,param);
  setup.Y0  = [H0 x0 y0 0];
  setup.T0  = 0;
  setup.Tend = 20*24*60*60;
  setup.Study_start_time = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  for ilight = 1:length(lightvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Cparam.daylight = lightvals(ilight) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    for itau = 1:length(tauvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
      Cparam.tau_c = tauvals(itau);
      setup.Tend = 20*24*60*60;
      [Y0,convergence_measure]     = remove_transient(param,Cparam,setup);
      param.onoff                  = 0;
      setup.Y0                     = Y0(1:4);
      setup.onoff                  = Y0(5);
      [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup);
      sleep_on_days  = sleep_on/(24*60*60);
      sleep_off_days = sleep_off/(24*60*60);
      sleepparam     = find_sleep_param_from_on_off_times(sleep_on_days,sleep_off_days,setup);
      setup.Y0       = Y0;
      param.onoff    = onoff;
%
%   Now run for two cycles and calculate Cmin
      setup.Tend = 2*24*60*60;
      [transitions] = default_protocol(setup);
      setup.transitions = transitions;
      [Htot,xtot,ytot,nreceptot,Ttot,param] = run_protocol(param,Cparam,setup);
      circadian_timeseries    = circadian(xtot,ytot,param);
      [TCmin,Cmin] = find_Cmin(Ttot/(24*60*60),circadian_timeseries);
      [TCmax,Cmax] = find_Cmin(Ttot/(24*60*60),-circadian_timeseries);
      [THmin,Hmin] = find_Cmin(Ttot/(24*60*60),Htot);
      [THmax,Hmax] = find_Cmin(Ttot/(24*60*60),-Htot);
%
% Plot to check
      if plt_check == 1
        figure(1)
        plot(Ttot/(24*60*60),Htot,'k')
        hold on
        plot(Ttot/(24*60*60),param.H0p+circadian_timeseries,'k')
        plot(Ttot/(24*60*60),param.H0m+circadian_timeseries,'k')
        plot(TCmin,param.H0m+Cmin,'ro')
        plot(TCmax,param.H0m-Cmax,'bo')
        plot(THmin,Hmin,'rx')
        plot(THmax,-Hmax,'bx')
        sleep_on  = sleepparam.mean_midsleep - 0.5*sleepparam.mean_sleep_duration;
        sleep_off = sleepparam.mean_midsleep + 0.5*sleepparam.mean_sleep_duration;
        plot([sleep_on sleep_off]/(24),[12 12],'b-')
        (THmin-TCmin)*24
        (THmax-TCmax)*24
        (TCmin-TCmax)*24+24
      end
%
% Add sleepparam in days, convert to hours for easier reading and save
%
      savesleepvals = [Cparam.daylight, Cparam.tau_c,sleepparam.mean_sleep_duration*24, sleepparam.std_sleep_duration*24, ...
      sleepparam.mean_midsleep*24, sleepparam.std_midsleep*24, TCmin(1)*24];
      fprintf(setup.fid4,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n',savesleepvals);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    end % End loop over tau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  end % End loop over light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
fclose(setup.fid1)
fclose(setup.fid4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
toc
