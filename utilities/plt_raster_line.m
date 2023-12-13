function plt_raster_line(sleep_times,wake_times,vshift)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes sleep and wake times (in days) and plots them as a raster plot.
% vshift shifts the position of the bar vertically.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
yyaxis left
%
for iday = 1:length(sleep_times)
  day = floor(wake_times(iday));
  plot([sleep_times(iday) wake_times(iday)]-day,[day day]+vshift,'k-','LineWidth',4)
  hold on
end
%
set(gca,'YDir','Reverse','TickDir','out','Xtick',-0.5:0.25:0.5,'XtickLabels',{'12:00','18:00','0:00','6:00'});
%
yyaxis right
xlim([-13.5 8.5])
%
end
