function [start_day,total_light,total_loglight,time_half_light,time_half_loglight,hours_brightlight,Cmax,TCmax_hr,Cmin,TCmin_hr,convergence_measure,xc0,x0,xcfinal,xfinal,xcchange,xchange,phase_initial,phase_final,phase_change] = find_daily_light_metrics(time_days_long,light_long,iphase,iphase_response)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the light metrics reported in the paper.
%
% ASSUMES TIMESERIES STARTS AT MIDNIGHT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate number of whole days
%
  ndays = floor(time_days_long(end)-time_days_long(1));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Works though each day and calculates metrics for each day.
  for iday = 1:ndays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select light and time window. 
    start_day(iday)      = time_days_long(1)+(iday-1);
    end_day              = time_days_long(1)+iday;
    time_iday            = time_days_long(time_days_long>=start_day(iday)&time_days_long<=end_day);
    light_iday           = light_long(time_days_long>=start_day(iday)&time_days_long<=end_day);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate mean of light and loglight
    total_light(iday)    = sum(light_iday)/length(light_iday);
    total_loglight(iday) = sum(log10(light_iday+1))/length(light_iday);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If at least some light is received, calculate more metrics 
    if total_loglight(iday)>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate time of half light
      [sorted_time,sorted_light,cumulative_light,cumulative_loglight,time_half_light(iday),time_half_loglight(iday)] ...
         = find_time_half_light(light_iday,time_iday);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate phase (converged phase assuming that the same light was received every day)
      if iphase == 1
        [Cmax(iday),TCmax_hr(iday),Cmin(iday),TCmin_hr(iday),convergence_measure(iday)] = find_phase(time_iday,light_iday);
      else
        Cmax(iday)                = 0;
        TCmax_hr(iday)            = 0;
        Cmin(iday)                = 0;
        TCmin_hr(iday)            = 0;
        convergence_measure(iday) = 0;
      end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the new light metric - called the "biological effect of light" in the paper
      if iphase_response == 1
        [xc0(iday),x0(iday),xcfinal(iday),xfinal(iday),xcchange(iday),xchange(iday),phase_initial(iday),phase_final(iday),phase_change(iday),x,xc,Tday] = find_phase_response(time_iday,light_iday);
      else
        xc0(iday)            = 0;
        x0(iday)             = 0;
        xcfinal(iday)        = 0;
        xfinal(iday)         = 0;
        xcchange(iday)       = 0;
        xchange(iday)        = 0;
        phase_initial(iday)  = 0;
        phase_final(iday)    = 0;
        phase_change(iday)   = 0;
      end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
      time_half_light(iday)    = 0;
      time_half_loglight(iday) = 0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the hours of bright light 
    hours_brightlight(iday)    = find_hours_bright_light(light_iday,500,time_iday*24*60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end % End loop over days
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
