function [time_secs,light] = loadlightdata(iparticipant,idevice,impute)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifically set up to read data in the folder data_raw
% Device 1 = Actiwatch
% Device 2 = HOBO
% impute 1 = imputed light
% impute 0 = raw light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  filename_light = sprintf('light_data_%i.xlsx',iparticipant);
  light_data = readtable(filename_light);
% 
  time_days  = datenum(light_data.DateTime);
%
  if idevice == 1
    if impute == 0
      light = light_data.WhiteLight;
    end
    if impute == 1
      light = light_data.ImputedWhiteLight;
    end
  end
  if idevice == 2
    if impute == 0
      light = light_data.LuxHobo2;
    end
    if impute == 1
      light = light_data.ImputedLuxHobo2;
    end
  end
%
  time_secs = time_days*(24*60*60);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
