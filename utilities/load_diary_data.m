function [data_mean_sleep_times,data_sleep_times] = load_diary_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Designed to read in diary data in folder data_raw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  data_raw = readtable('sleep_diary_data.xlsx');
%
  data_raw.date_finalawake      = datenum(data_raw.date_on_diary);
  data_raw.date_gotosleep       = data_raw.date_finalawake;
%
% Now correct dates since diary records date of day on which participant wakes up
%
  gotosleep                     = data_raw.gotosleep;
  finalawake                    = data_raw.finalawake;
%
  data_raw.date_gotosleep((finalawake-gotosleep)<0)  = data_raw.date_gotosleep((finalawake-gotosleep)<0) - 1;
%
  data_sleep_times.participant_no = data_raw.participant_no;
  data_sleep_times.sontime        = data_raw.date_gotosleep   + data_raw.gotosleep + data_raw.fallasleep/(60*24);
  data_sleep_times.sofftime       = data_raw.date_finalawake  + data_raw.finalawake;  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% And calculate mean sleep and sleep duration
  participants  = unique(data_sleep_times.participant_no);
  nparticipants = length(participants);
%
  for inum = 1:nparticipants
    ipart = participants(inum);
%
    idx_ipart      = find(data_sleep_times.participant_no == ipart & ~isnan(data_sleep_times.sontime) ...
      & ~isnan(data_sleep_times.sofftime));
    sontime  = data_sleep_times.sontime(idx_ipart);
    sofftime = data_sleep_times.sofftime(idx_ipart);
    diary_SD = (sofftime-sontime)*24;
    diary_MS = mod(sofftime,1)*24-0.5*diary_SD;
%
% Note the means only include the field days. The last diary day corresponds to the laboratory day.
%
    data_mean_sleep_times.participant_no(inum)      = ipart;
    data_mean_sleep_times.mean_sleep_duration(inum) = mean(diary_SD(1:end-1));
    data_mean_sleep_times.std_sleep_duration(inum)  = std(diary_SD(1:end-1));
    data_mean_sleep_times.mean_midsleep(inum)       = mean(diary_MS(1:end-1));
    data_mean_sleep_times.std_midsleep(inum)        = std(diary_MS(1:end-1));
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
