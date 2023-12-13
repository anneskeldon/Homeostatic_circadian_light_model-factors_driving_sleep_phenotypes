function [sleepparam] = find_sleep_param_from_on_off_times(sleep_on,sleep_off,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From a list of sleep onset and offset times (given in days) calculate a variety of metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n_on           = length(sleep_on);
  n_off          = length(sleep_off);
%
  if sleep_on(1) < sleep_off(1)
    if sleep_on(end) < sleep_off(end)
      sleep_onsets  = sleep_on(1:n_on)' + setup.Study_start_time/(24*60*60);
      sleep_offsets = sleep_off(1:n_on)' + setup.Study_start_time/(24*60*60);
    else
      sleep_onsets  = sleep_on(1:n_off)' + setup.Study_start_time/(24*60*60);
      sleep_offsets = sleep_off(1:n_off)' + setup.Study_start_time/(24*60*60);
    end
  else
    if sleep_on(end) < sleep_off(end)
      sleep_onsets  = sleep_on(1:n_on-1)' + setup.Study_start_time/(24*60*60);
      sleep_offsets = sleep_off(2:n_on)' + setup.Study_start_time/(24*60*60);
    else
      sleep_onsets  = sleep_on(1:n_off-1)' + setup.Study_start_time/(24*60*60);
      sleep_offsets = sleep_off(2:n_off)' + setup.Study_start_time/(24*60*60);
    end
  end
%
  if ~isempty(sleep_onsets) & ~isempty(sleep_offsets)
    sleepparam.sleep_onsets        = sleep_onsets;
    sleepparam.sleep_offsets       = sleep_offsets;
    sleepparam.sleep_durations     = sleep_offsets-sleep_onsets;
    sleepparam.midsleeps           = mod(sleep_offsets,1) - 0.5*sleepparam.sleep_durations;
    sleepparam.last_sleep_onset    = sleep_onsets(end);
    sleepparam.last_sleep_offset   = sleep_offsets(end);
    sleepparam.last_sleep_duration = sleep_offsets(end) - sleep_onsets(end);
    sleepparam.last_midsleep       = mod(sleep_offsets(end),1) - 0.5*sleepparam.sleep_durations(end);
    sleepparam.mean_sleep_duration = mean(sleepparam.sleep_durations);
    sleepparam.std_sleep_duration  = std(sleepparam.sleep_durations);
    [mean_MS,std_midsleep_csd]     = circ_mean_std(mod(sleepparam.midsleeps,1)*24);
    sleepparam.mean_midsleep       = mean_MS/24;     % in days
      if sleepparam.mean_midsleep > 0.5
        sleepparam.mean_midsleep   = sleepparam.mean_midsleep - 1;
      end
    sleepparam.std_midsleep        = std_midsleep_csd/24; % in days
  else
    sleepparam.sleep_onsets        = -100;
    sleepparam.sleep_offsets       = -100;
    sleepparam.sleep_durations     = -100;
    sleepparam.midsleeps           = -100;
    sleepparam.last_sleep_onset    = -100;
    sleepparam.last_sleep_offset   = -100;
    sleepparam.last_sleep_duration = -100;
    sleepparam.last_midsleep       = -100;
    sleepparam.mean_sleep_duration = -100;
    sleepparam.std_sleep_duration  = -100;
    sleepparam.mean_midsleep       = -100;
    sleepparam.std_midsleep        = -100;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
