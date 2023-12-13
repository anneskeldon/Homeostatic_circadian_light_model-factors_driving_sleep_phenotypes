%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in fitted parameters and light and sleep timing data for the cohort studied in the paper
% and constructs sleep on / off times and some graphs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  addpath('../data_raw')
  addpath('../data_generated/fitted_parameters/')
  addpath('../utilities')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Files containing fitted parameters
  filename_mu_ca  = 'fitted_values_HCL_SDMS_Actiwatch_imputed_mu_ca_20231028.dat';
  filename_mu_tau = 'fitted_values_HCL_SDMS_Actiwatch_imputed_mu_tau_20231118.dat';	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Read in parameter defaults
  [Cparam,param] = default_HCL_parameters(30);
  [setup]        = default_setup_parameters();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% And set up some others. 
  Cparam.gated      = 0;
  Cparam.lightinput = 1; % 0 Means use the tanh profile in swlightinput.m, 
                         % 1 means use light data in Cparam.timedata and Cparam.lightdata
  mu_tau            = 0; % Decide if looking at (mu,tau) fittings or (mu,ca) fittings. 1 means (mu,tau).
  idevice           = 1; % Using Actiwatch
  impute            = 1; % Using imputed data
%
% Read in fitted parameters 
%
  if mu_tau == 1
    param_fittings = readmatrix(filename_mu_tau);
    participants   = param_fittings(:,1);
  else
    param_fittings = readmatrix(filename_mu_ca);
    idx_vals       = find( (param_fittings(:,5) < 0.03) ...
      & param_fittings(:,2) > 0 ...
      & param_fittings(:,3) + 13 < param_fittings(:,2) );
    participants   = param_fittings(idx_vals,1)';
  end
  nparticipant   = length(participants);
  predicted_sleepon_sleepoff_tau = []; 
  predicted_sleepon_sleepoff_ca  = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Now loop over participants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  for jparticipant = 1:nparticipant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Select participant
    iparticipant = participants(jparticipant); 
% Find parameters for that participant
    idx = find(param_fittings(:,1) == iparticipant);
    param.mu     = param_fittings(idx,2);
    param.a      = param_fittings(idx,3);
    Cparam.tau_c = param_fittings(idx,4);
% Set up output files for the selected participant
    rundata_events     = sprintf('temp/events.dat');
    rundata_timeseries = sprintf('temp/timeseries_mu_tau_%i_participant_%i_from_20231108_fittings.dat',mu_tau,iparticipant);
    setup.fid1 = fopen(rundata_events,'w');
    setup.fid4 = fopen(rundata_timeseries,'w');
% Read in light data for the selected participant
    [time,light] = load_light_data(iparticipant,idevice,impute);
    setup.Study_start_time             = time(1);
    time                               = time-setup.Study_start_time; %Shift so timeseries starts at 0
% Set up initial conditions etc for the integrator
    setup.T0            = time(1);
    setup.Tend          = floor(time(end)/(24*60*60))*24*60*60;
    xc0                 = -0.81;
    yc0                 = -0.66;
    H0                  = param.H0p+circadian(xc0,yc0,param);
    setup.Y0            =  [H0 xc0 yc0 0];
    Cparam.timedata     = time;
    Cparam.lightdata    = light;
    [setup.transitions] = default_protocol(setup);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Plot raster plot
    figure(iparticipant)
    plt_raster_line(sleepparam.sleep_onsets,sleepparam.sleep_offsets,0)
    xlim([-1 1])
    datetick('x')
%
    if mu_tau == 1
      predicted_sleepon_sleepoff_tau = [predicted_sleepon_sleepoff_tau; ...
        repmat(iparticipant,length(sleepparam.sleep_onsets),1),sleepparam.sleep_onsets,sleepparam.sleep_offsets];
    else
      predicted_sleepon_sleepoff_ca  = [predicted_sleepon_sleepoff_ca; ...
        repmat(iparticipant,length(sleepparam.sleep_onsets),1),sleepparam.sleep_onsets,sleepparam.sleep_offsets];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculate timeseries
    [Htot,xtot,ytot,nreceptot,Ttot,param] = run_protocol(param,Cparam,setup);
% Plot two-process model
    figno = 100+iparticipant;
    plt_HCL_output(Htot,xtot,ytot,nreceptot,Ttot,figno,param,Cparam);
    Ttot = Ttot+setup.Study_start_time;
% Save timeseries
    save_timeseries(Htot,xtot,ytot,nreceptot,Ttot,param,Cparam,setup);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Cross check
    sleep_MS_SD_check(iparticipant,1) = sleepparam.mean_midsleep*24;
    sleep_MS_SD_check(iparticipant,2) = sleepparam.mean_sleep_duration*24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  end % End loop over participant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  if mu_tau == 1
    writematrix(predicted_sleepon_sleepoff_tau);
  else
    writematrix(predicted_sleepon_sleepoff_ca);
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
