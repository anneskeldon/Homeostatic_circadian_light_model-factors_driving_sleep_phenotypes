function plt_HCL_output(H,xpace,ypace,nrecep,t,figno,param,Cparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots out the HCL model and the various associated variables. Uses figures
% figno and figno+1.
%
% H is the homeostatic sleep pressure
% xpace and ypace are the components of the van der Pol oscillator model
% nrecep is the fraction of activated photoreceptors
% t is time
% param and Cparam are structures which contain the parameter values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
  figure(figno+1)
%
  lower_threshold = param.H0m + circadian(xpace,ypace,param);
  upper_threshold = param.H0p + circadian(xpace,ypace,param);
  Tplot           = t/(60*60); %Plot time in hours
%
  subplot(6,1,1)
  plot(Tplot,H,'k')
  hold on
 plot(Tplot,lower_threshold,'b')
 plot(Tplot,upper_threshold,'b')
%
  axis([0 Tplot(end) 5 25]);
  subplot(6,1,2)
  plot(Tplot,xpace,'b')
  hold on
  subplot(6,1,3)
  plot(Tplot,ypace,'r')
 
  subplot(6,1,4)
  plot(Tplot,circadian(xpace,ypace,param)/param.a,'k')
%
  subplot(6,1,5)
  light = swlightinput(t,Cparam);
  plot(Tplot,nrecep,'k')
  hold on
  subplot(6,1,6)
  hold on
  plot(Tplot,light,'k')
%
  figure(figno)
  lw = 2;
  subplot(3,1,1)
  plot(Tplot/24-14,light,'-','Color',[1 0.9 0.1],'LineWidth',lw)
  set(gca,'TickDir','out','XTickLabels',{})
  subplot(3,1,2)
  plot(Tplot/24-14,xpace,'b','LineWidth',2)
  hold on
  plot(Tplot/24-14,ypace,'r','LineWidth',2)
  set(gca,'TickDir','out','XTickLabels',{})
  subplot(3,1,3)
  plot(Tplot/24-14,H,'k','LineWidth',2)
  set(gca,'TickDir','out')
%
end
