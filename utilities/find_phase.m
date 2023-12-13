function [Cmax,TCmax_hr,Cmin,TCmin_hr,convergence_measure]=find_phase(time_days,light)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Designed to find the phase as a result of a single (repeated) day of light data. 
% Assumes data is from t=0 for 24 hours and uses the Forger et al 1999 light-circadian model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Lparam.time_secs = (time_days-time_days(1))*24*60*60;
  Lparam.light     = light;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify up model parameters
  Cparam.b=0.4; Cparam.G=19.9; Cparam.alpha0=0.16; Cparam.p=0.6; Cparam.I0=9500; 
  Cparam.kappa=12/pi*60*60; Cparam.lambda=60/(60*60); Cparam.mu=0.23; Cparam.tau_c=24.20; 
  Cparam.betaval=0.013; Cparam.f=0.99669; Cparam.k=0.55;
%
  nloops     = 200;  % Maximum number of loops of integration for removal of transient. 
  ndays      = 1;  % Number of days in each loop.
%
% Initial conditions
  x0=0.8805; xc0=0.5574; n0=0.9270; 
%
% model
  odefun=@(t,Y)forger1999model(t,Y,Cparam,Lparam);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First integrate to remove transient
%
  resnorm = 10;
  iloop   = 1;
  while (iloop<nloops) & (resnorm>0.001) 
    [T,Y] = ode15s(odefun,[0 ndays*24*60*60],[x0 xc0 n0]);
    resnorm(iloop) = sqrt(sum((([x0 xc0 n0] - Y(end,:))./Y(end,:)).^2));
%
% look at convergence
%    figure(10)
%    subplot(1,2,1)
%    plot(iloop,log(resnorm(iloop)),'x')
%    hold on
%    subplot(1,2,2)
%    hold on
%    plot(iloop,x0,'bo')
%    plot(iloop,xc0,'ko')
%    plot(iloop,n0,'go')
%
    x0=Y(end,1); xc0=Y(end,2); n0=Y(end,3);
    iloop=iloop+1;
  end
  nresnorm = length(resnorm);
  convergence_measure = sum(resnorm(end:end-9)); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now integrate for ndays and then calculate the maximum and minimum
  for i=1:ndays
    options=odeset('Maxstep',60,'RelTol',1e-4,'BDF','on','AbsTol',1e-4);
%
    [T,Y]=ode15s(odefun,[0 24*60*60],[x0 xc0 n0],options);
    x=Y(:,1); xc=Y(:,2); nrecep=Y(:,3); Thour=T/(60*60); Tday=T/(24*60*60);
%
% Find maximum and minimum
    [Cmax,TCmax_ind]=max(x);
    [Cmin,TCmin_ind]=min(x);
%
    Tindmax = floor(mean(TCmax_ind));
    Tindmin = floor(mean(TCmin_ind));
    TCmax   = T(Tindmax);
    TCmin   = T(Tindmin);
%
    TCmax_hr = mod(TCmax/(60*60) + time_days(1)*24,24);
    TCmin_hr = mod(TCmin/(60*60) + time_days(1)*24,24);
%
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dy=forger1999model(t,y,Cparam,Lparam)
%
% Forger et al 1999 light-circadian model https://doi.org/10.1177/074873099129000867
%
  x      = y(1);
  xc     = y(2);
  nrecep = y(3);
%
% Process L
  alpha = Cparam.alpha0*(lightinput(t,Lparam)/Cparam.I0)^Cparam.p;
  Bhat  = Cparam.G*alpha*(1-nrecep);
  B     = (1-Cparam.b*x)*(1-Cparam.b*xc)*Bhat;
%
% Differential equations
  dxdt  = (xc + B )/Cparam.kappa ; 
  dxcdt = (Cparam.mu*(xc-4*xc^3/3) ...
    -((24/(Cparam.f*Cparam.tau_c))^2+Cparam.k*B)*x)/Cparam.kappa; 
%
  dnrecepdt = Cparam.lambda*(alpha*(1-nrecep)-Cparam.betaval*nrecep);
%
  dy = zeros(3,1);
  dy(1) = dxdt;
  dy(2) = dxcdt;
  dy(3) = dnrecepdt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lightinput=lightinput(t,Lparam)
  lightinput=interp1(Lparam.time_secs,Lparam.light,t,'linear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
