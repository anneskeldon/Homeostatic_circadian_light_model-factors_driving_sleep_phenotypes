function run_HCL_parameter_fitting(idevice,idata_type,ifit_param,fit_method,participants,date)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in light and sleep timing data (provided in the folder data_raw) and fits parameters of the HCL 
% model.
% Designed specifically to run with data used in the paper, but could readily be adapted to use
% data from other sources.
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  device    = {'Actiwatch','HOBO'};  % i.e. light from the wrist (Actiwatch) or shoulder (HOBO)
  data_type = {'raw','imputed'};     % light data are raw or imputed
  fit_param = {'mu_tau','mu_ca'};    % choices of which parameters to fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up parameters 
%
 [setup] = default_setup_parameters();
 setup.fit_param  = ifit_param;
 setup.plt_raster = 0;
 impute           = idata_type-1;
%
  rundata_fittedvals_SDMS = sprintf('./temp/fitted_values_HCL_SDMS_%s_%s_%s_%i.dat',device{idevice},data_type{idata_type},fit_param{ifit_param},date)
  fid5                    = fopen(rundata_fittedvals_SDMS,'w');
%
% Parameter defaults
  [Cparam,param] = default_HCL_parameters(65);
  Cparam.gated   = 0;
%
%param.a      = 2.44; % high
%param.a      = 1.00; % low
%
  setup.nloops    = 20; % Maximum number of loops through data for removal of transients in calculation of objective functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Specify that light will be taken from a file
%
  Cparam.lightinput    = 1; % 0 Means use the tanh profile in swlightinput.m, 1 means use light data in Cparam.timedata and Cparam.lightdata
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Read in sleep data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
  [data_mean_sleep_times,data_sleep_times] = load_diary_data;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Loop over participants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  for iparticipant = participants
    iparticipant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Select participant sleep duration and midsleep
    setup.diary_SD = data_mean_sleep_times.mean_sleep_duration(data_mean_sleep_times.participant_no==iparticipant);
    setup.diary_MS = data_mean_sleep_times.mean_midsleep(data_mean_sleep_times.participant_no==iparticipant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Read in light data
    [time,light]                       = load_light_data(iparticipant,idevice,impute);
%
    setup.Study_start_time             = time(1);
    time                               = time-setup.Study_start_time; %Shift so timeseries starts at 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Initialise parameters for integration
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
% Fit iteratively, repeatedly looping through fitting for mu then tau until sufficient accuracy
% is reached (or been through the iterative loop 5 times).
%
    if ifit_param == 1 & fit_method == 'iterative'  
%
      fitparam.maxloops            = 10;
      mu_guess                     = interp1([4.5 12],[14.5 25],setup.diary_SD)
      tau_guess                    = 24.2;
      residual                     = 1;
      iloop                        = 1;
% Although remove transient at every stage in the calculation, do it here to try to get a better initial
% condition to start with and speed up later states.
      [Y0,convergence_measure]     = remove_transient(param,Cparam,setup);
      setup.Y0                     = Y0(1:4)
      param.onoff                  = Y0(5)
%
      while residual > 0.005 & iloop < 5;
%
% Find mu
        fitparam.fit_param_MS        = 0;
        fitparam.x0_bound1           = mu_guess;;
        fitparam.delta               = 1;
        [mu_new,fval_SD,exitflag_SD] = find_fit_param_SD_or_MS(fitparam,param,Cparam,setup)
        param.mu                     = mu_new
        mu_guess                     = mu_new;
%
% Now that have a good guess for mu, again remove transient, to get a better initial
% condition for finding tau
%
      [Y0,convergence_measure]     = remove_transient(param,Cparam,setup);
      setup.Y0                     = Y0(1:4);
      param.onoff                  = Y0(5);
%     
% Find tau
        fitparam.fit_param_MS        = 1;
        fitparam.x0_bound1           = tau_guess;
        fitparam.delta               = 0.1;
        [tau_new,fval_MS,exitflag]   = find_fit_param_SD_or_MS(fitparam,param,Cparam,setup)
        Cparam.tau_c                 = tau_new
        tau_guess                    = tau_new;
% Find overall residual
        setup.ifig       = iparticipant;
        x0               = [param.mu,Cparam.tau_c];
        errors           = objfunc_fsolve_meanSDmeanMS(x0,param,Cparam,setup)
        residual         = sqrt(sum(errors.^2))
%
        iloop                        = iloop + 1;
      end
       amp_guess = param.a;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Fit for two parameters simultaneously. 
%
    if fit_method == 'twoparams'  
%
      setup.plt_raster  = 0
      mu_guess          = param.mu;
      amp_guess         = param.a;
      tau_guess         = Cparam.tau_c;

      if setup.fit_param == 1  % For fitting mu and tau
        x0_guess          = [mu_guess tau_guess];
      end
%
      if setup.fit_param == 2  % For fitting mu and amp
        x0_guess          = [mu_guess amp_guess];
      end
%
      options           = optimset('TolX',1e-3,'TolFun',1e-3);
      [x,residual,exitflag,output] = fminsearch(@(x0) objfunc_meanSDmeanMS(x0,param,Cparam,setup),x0_guess,options)
%      [x,residual,exitflag,output] = fminunc(@(x0) objfunc_meanSDmeanMS(x0,param,Cparam,setup),x0_guess,options)
%       
      if setup.fit_param == 1
        mu_guess  = x(1);
        tau_guess = x(2);
      end
%
      if setup.fit_param == 2
        mu_guess  = x(1);
        amp_guess = x(2);
      end
%
% check
%      setup.plt_raster = 1;
      setup.ifig       = 100+iparticipant;
      errors           = objfunc_fsolve_meanSDmeanMS(x,param,Cparam,setup)
      residual         = sqrt(sum(errors.^2))
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Save output
%
    fprintf(fid5,'%d %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n',iparticipant,mu_guess,amp_guess,tau_guess,residual,setup.diary_MS,setup.diary_MS+errors(1),setup.diary_SD,setup.diary_SD+errors(2))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  end % End loop over participant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  fclose(fid5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
end
