%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Code to calculate light metrics 
% idevice = 1 Actiwatch (i.e. light data collected from the wrist)
% idevice = 2 HOBO (i.e. light data collected from the shoulder)
%
  addpath('../data_raw/')
  addpath('../utilities/')
%
tic
%
% Set up some choices
  write_table      = 1;  % Specifies whether or not to write out the table that is generated
  iphase           = 1;  % Calculate the phase? 1 means yes.
  iphase_response  = 1;  % Calculate the phse response? ("biological effect of light" in the paper). 1 means yes.
  impute           = 1;  % 0 for raw light data; 1 for imputed
%
  iloop  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose device
  for idevice = 1%:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose participants
    if idevice == 1
      participants = [1:21 23:35];
    end 
    if idevice == 2
      participants = [1 2 4 7 11 12 13 17 20 24 25 29 30 31 32 34];
    end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each participant calculate light metrics
    for iparticipant = participants
      tic
      iparticipant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load participant light data
      [time,light]     = load_light_data(iparticipant,idevice,impute);
%
      time_days_long   = time/(24*60*60);
      light_long       = light;
%
% Truncate so that the timeseries starts at midnight; 
      start_day        = ceil(time(1)/(24*60*60));
      time_days_short  = time_days_long(time_days_long>=start_day);
      light_short      = light_long(time_days_long>=start_day);
%
% Work out light metrics - note light metrics will only be calculated for each complete 24 hour period starting
% at midnight.
      [start_day,total_light,total_loglight,time_half_light,time_half_loglight,hours_brightlight,Cmax,TCmax_hr,Cmin,TCmin_hr,convergence_measure,xc0,x0,xcfinal,xfinal,xcchange,xchange,phase_initial,phase_final,phase_change] = find_daily_light_metrics(time_days_short,light_short,iphase,iphase_response);
%
% Save metrics in a table 
      device       = repmat(idevice,length(start_day),1);
      participant  = repmat(iparticipant,length(start_day),1);
        if iloop == 1
          light_table      = table(device,participant,start_day',total_light',total_loglight',time_half_light',time_half_loglight',hours_brightlight',Cmax',TCmax_hr',Cmin',TCmin_hr',convergence_measure',xc0',x0',xcfinal',xfinal',xcchange',xchange',phase_initial',phase_final',phase_change');
        else
          sub_light_table  = table(device,participant,start_day',total_light',total_loglight',time_half_light',time_half_loglight',hours_brightlight',Cmax',TCmax_hr',Cmin',TCmin_hr',convergence_measure',xc0',x0',xcfinal',xfinal',xcchange',xchange',phase_initial',phase_final',phase_change');
          light_table = [light_table;sub_light_table];
        end
      iloop = iloop+1;
      toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % End loop over participant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end % End loop over device
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out table
  if write_table == 1
    light_table.Properties.VariableNames = ["device","participant","start_day","total_light","total_loglight","time_half_light","time_half_loglight",...
       "hours_brightlight","Cmax","TCmax_hr","Cmin","TCmin_hr","convergence_measure",...
       "x0","y0","xfinal","yfinal","xchange","ychange","phase_initial","phase_final","phase_change"];
    if impute == 0
      light_table_raw = light_table;
      writetable(light_table_raw)
    end
    if impute == 1
      light_table_imputed = light_table;
      writetable(light_table_imputed)
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
