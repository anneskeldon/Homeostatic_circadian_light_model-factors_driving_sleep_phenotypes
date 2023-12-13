function [sorted_time,sorted_light,cumulative_light,cumulative_loglight,time_half_light,time_half_loglight]=find_time_half_light(light,time_days)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given one day of data, calculate the time of half light. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  timeofday = mod(time_days,1);
% Now sort times and light values based on time of day
  [sorted_time,sort_ind] = sort(timeofday,'ascend');
  sorted_light    = light(sort_ind);
  total_light     = sum(sorted_light);
  total_log_light = sum(log10(sorted_light+1));
%
  for i=1:length(sorted_light)
    cumulative_light(i)    = sum(sorted_light(1:i));
    cumulative_loglight(i) = sum(log10(sorted_light(1:i)+1));
  end
  ind_light    = find(cumulative_light>total_light/2);  
  ind_loglight = find(cumulative_loglight>total_log_light/2);  
%
  time_half_light    = sorted_time(ind_light(1));
  time_half_loglight = sorted_time(ind_loglight(1));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
