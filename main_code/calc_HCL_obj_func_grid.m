%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in light and sleep timing data and calculates the objective function where the objective function
% is 
% sqrt( (mean_sleep_duration_simulated - mean_sleep_duration_observed)^2 
%    + (mean_mid-sleep_simulated - mean_mid-sleep_observed)^2 )
%
% Here set up to calculate the objective function for a grid of points (to create the data for Fig. 6b in the
% paper).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  addpath('../data_raw')
  addpath('../utilities')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Parameter defaults
  [Cparam,param] = default_HCL_parameters(65);
  [setup]        = default_setup_parameters();
  Cparam.gated   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Need to make choices on parameters
%
  Cparam.lightinput    = 1; % 0 Means use the tanh profile in swlightinput.m, 
                            % 1 means use light data in Cparam.timedata and Cparam.lightdata
  idevice          = 1;  % 1 is Actiwatch, 2 is HOBO
  impute           = 1;  % 1 is imputed, 0 is raw
  mu_tau           = 0;  % 1 for (mu,tau). any other value for (mu,ca)
  setup.nloops     = 20; % Maximum number of loops for removal of transient.
  setup.plt_raster = 0;
  setup.ifig       = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Read in diary data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  [data_mean_sleep_times,data_sleep_times] = load_diary_data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Choose participant
  participant_list = 32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  for jparticipant = 1:length(participant_list)
    iparticipant = participant_list(jparticipant);
    if mu_tau == 1
      fid5 = fopen(sprintf('./temp/obj_func_mu_tau_%i.dat',iparticipant),'w');
    else
      fid5 = fopen(sprintf('./temp/obj_func_mu_amp_%i.dat',iparticipant),'w');
    end
    iparticipant
    setup.diary_SD = data_mean_sleep_times.mean_sleep_duration(data_mean_sleep_times.participant_no==iparticipant);
    setup.diary_MS = data_mean_sleep_times.mean_midsleep(data_mean_sleep_times.participant_no==iparticipant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Read in light data
% Original excel files for light
    [time,light]                       = load_light_data(iparticipant,idevice,impute);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    setup.Study_start_time             = time(1);
    time                               = time-setup.Study_start_time; %Shift so timeseries starts at 0
% Set up initial conditions 
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
% Define grid of parameter values for calculation of the objective function
    tauvals     = 23.7:0.05:24.5;
    ampvals     = 0.5:0.5:10; 
    if mu_tau == 1
       muvals = 13.0:0.5:22.5; 
    else
       muvals = 13.0:0.25:22.5; 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Loop over mu
    for imu = 1:length(muvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Loop over tauvals / cavals
      param.mu   = muvals(imu);
      if mu_tau == 1
        param_vals = tauvals;
      else
        param_vals = ampvals;
      end
      for iparam = 1:length(param_vals)
        if mu_tau == 1
          itau            = iparam;
          Cparam.tau_c    = tauvals(itau);
          setup.fit_param = 1;
          x0              = [param.mu,Cparam.tau_c];
        else
          iamp           = iparam;
          param.a         = ampvals(iamp)
          setup.fit_param = 2;
          x0              = [param.mu,param.a];
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculate residual
        [errors] = objfunc_fsolve_meanSDmeanMS(x0,param,Cparam,setup)
        residual = sqrt(sum(errors.^2))
        residual_tot(imu,iparam) = residual;
% Save
        mu  = param.mu;
        amp = param.a;
        tau = Cparam.tau_c;
        fprintf(fid5,'%d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n',iparticipant,mu,amp,tau,residual,errors(1),errors(2))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
      end % End loop over inner parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    end % End loop over mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  end % End loop over participant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  fclose(fid5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure()
    x       = muvals;
    if mu_tau == 1    
      y       = tauvals;
    else
      y       = ampvals;
    end
  [XX,YY] = meshgrid(x,y);
  contourf(XX,YY,log10(residual_tot)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
toc
