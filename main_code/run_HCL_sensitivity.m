%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the HCL model repeatedly in spontaneous mode with a synthetic light profile. Loops over
% parameter values - set up to to create the data for the sensitivity plots in the paper.
%
  addpath('../utilities')
%
  tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  date = 16092023;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Select which parameters to loop over 
%
  param_name = {'mu','H0bar','Delta','Camp','tau','G','exponent','chi','alpha0','b','gamma','k','beta','I0'};
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for iparam_name = 1:14  %Loop over the different parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set HCL model parameters to the defaults
    age                  = 65;
    [Cparam,param]       = default_HCL_parameters(age);
    [setup]              = default_setup_parameters();
    Cparam.lightinput    = 0; % Means use the tanh profile in swlightinput.m
    Cparam.daylight      = 700;
    Cparam.eveninglight  = 10; % Default value 40.
    H0p                  = param.H0p;
    H0m                  = param.H0m;
%
% Set up ranges that parameters take
%
    if iparam_name == 1
      paramvals = [13:0.2:24];
    end
    if iparam_name == 2
      paramvals = [-5:0.2:4];
    end
    if iparam_name == 3
      paramvals = [0.5:0.1:1.5]; 
    end
    if iparam_name == 4
      paramvals = [0.5:0.01:6];
    end
    if iparam_name == 5
      paramvals = [23.5:0.01:24.7];
    end
    if iparam_name == 6
      paramvals = [8:0.5:30];
    end
    if iparam_name == 7
      paramvals = [0.4:0.01:1]; 
    end
    if iparam_name == 8
      paramvals = [5:1:25 30:65]; 
    end
    if iparam_name == 9
      paramvals = [0.05:0.01:0.25]; 
    end
    if iparam_name == 10
      paramvals = [0:0.1:1]; 
    end
    if iparam_name == 11
      paramvals = [0.15:0.01:0.35]; 
    end
    if iparam_name == 12
      paramvals = [0.25:0.05:0.85]; 
    end
    if iparam_name == 13
      paramvals = [0.005:0.001:0.020]; 
    end
    if iparam_name == 14
      paramvals = [500:500:40000]; 
    end
%
% File where results will be saved
%
    rundata_events            = ('./temp/rundata_events.dat')
    setup.fid1                = fopen(rundata_events,'w');
    rundata_param_sensitivity = sprintf('./temp/param_sensitivity_age65_param_%s_%d_%i.dat',param_name{iparam_name},Cparam.eveninglight,date);
    setup.fid4                = fopen(rundata_param_sensitivity,'w');
    save_header(param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loop over parameter values and calculate metrics
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iparam = 1:length(paramvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      if iparam_name == 1
        param.mu  = paramvals(iparam);
      end
      if iparam_name == 2
        param.H0p  = H0p + paramvals(iparam);
        param.H0m  = H0m + paramvals(iparam);
      end
      if iparam_name == 3
        param.H0p  = (H0p+H0m)/2 + paramvals(iparam)/2
        param.H0m  = (H0p+H0m)/2 - paramvals(iparam)/2;
      end
      if iparam_name == 4
        param.a  = paramvals(iparam);
      end
      if iparam_name == 5
        Cparam.tau_c  = paramvals(iparam);
      end
      if iparam_name == 6
        Cparam.G  = paramvals(iparam);
      end
      if iparam_name == 7
        Cparam.p  = paramvals(iparam);
      end
      if iparam_name == 8
        param.Xs  = paramvals(iparam)*60*60;
        param.Xw  = paramvals(iparam)*60*60;
      end
      if iparam_name == 9
        Cparam.alpha0  = paramvals(iparam);
      end
      if iparam_name == 10
        Cparam.b  = paramvals(iparam);
      end
      if iparam_name == 11
        Cparam.gam  = paramvals(iparam);
      end
      if iparam_name == 12
        Cparam.k  = paramvals(iparam);
      end
      if iparam_name == 13
        Cparam.betaval  = paramvals(iparam);
      end
      if iparam_name == 14
        Cparam.I0  = paramvals(iparam);
      end
%
      saveparam = paramvals(iparam);
%
% Set initial conditions for integrator
%
      x0 = -0.81;
      y0 = -0.66;
      H0 = param.H0p + circadian(x0,y0,param);
      param.onoff = 0;  
      setup.Y0    = [H0 x0 y0 0];
      setup.T0   = 0;
      setup.Tend = 20*24*60*60; % Run for 20 days
%
% Set up spontaneous protocol 
%
      [transitions]          = default_protocol(setup); 
      setup.transitions      = transitions;
      setup.Study_start_time = setup.T0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Remove transient
      setup.nloops                 = 20;
      [Y0,convergence_measure]     = remove_transient(param,Cparam,setup);
      setup.Y0                     = Y0(1:4);
      param.onoff                  = Y0(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Now find sleep onset/offset times.
      [sleep_on,sleep_off,Y0,onoff]  = find_sleeps(param,Cparam,setup);
      sleep_on_days  = sleep_on/(24*60*60);
      sleep_off_days = sleep_off/(24*60*60);
      [sleepparam] = find_sleep_param_from_on_off_times(sleep_on_days,sleep_off_days,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Now calculate Cmin, Hmin and Hmax by integrating for 2 further cycles. Use run_protocol since
% this returns the timeseries for all variables.  
%
      setup.T0          = 0;
      setup.Tend        = 2*24*60*60;
      [transitions]     = default_protocol(setup);
      setup.transitions = transitions;
      [Htot,xtot,ytot,nreceptot,Ttot,param] = run_protocol(param,Cparam,setup);
      circadian_timeseries = circadian(xtot,ytot,param);
      [TCmin,Cmin]      = find_Cmin(Ttot/(24*60*60),circadian_timeseries);
      [TCmax,Cmax]      = find_Cmin(Ttot/(24*60*60),-circadian_timeseries);
      [THmin,Hmin]      = find_Cmin(Ttot/(24*60*60),Htot);
      [THmax,Hmax]      = find_Cmin(Ttot/(24*60*60),-Htot);
% 
% Plot to check
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
%
      savesleepvals = [saveparam, sleepparam.mean_sleep_duration*24, sleepparam.std_sleep_duration*24, ...
        sleepparam.mean_midsleep*24, sleepparam.std_midsleep*24, TCmin(1)*24]  
%
      fprintf(setup.fid4,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n',savesleepvals');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % End loop over parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fclose(setup.fid4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end % End loop over different parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
