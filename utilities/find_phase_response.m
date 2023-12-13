function [xc0,x0,xcfinal,xfinal,xcchange,xchange,phase_initial,phase_final,phase_change,x,xc,Tday] = find_phase_response(time_days,light)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the phase response of the inputted light. When the light is for one day, this is the
% metric called the "biological effect of light" in the paper 
% Essentially the same as find_phase.m, but there is no transient removal (could code more neatly
% as a single function with an option for whether or not a transient is removed).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Lparam.time_secs = (time_days-time_days(1))*24*60*60;
  Lparam.light     = light;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up parameters
%
  Cparam.b=0.4; Cparam.G=19.9; Cparam.alpha0=0.16; Cparam.p=0.6; Cparam.I0=9500; 
  Cparam.kappa=12/pi*60*60; Cparam.lambda=60/(60*60); Cparam.mu=0.23; Cparam.tau_c=24.20; 
  Cparam.betaval=0.013; Cparam.f=0.99669; Cparam.k=0.55;
%
  ndays      = 1;  % Number of days in each loop.
%
% Initial conditions (note that in HCL, xc is equiv of x, x is equiv of y).
 xc0 = -0.8738; x0 = -0.5534; n0 = 0.3181; % For midnight
%
% Set up model
  odefun=@(t,Y)forger1999model(t,Y,Cparam,Lparam);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate circadian model
%
% Integrate for ndays where ndays = 1
  for i = ndays
    options=odeset('Maxstep',60,'RelTol',1e-4,'BDF','on','AbsTol',1e-4);
%
    [T,Y]=ode15s(odefun,[0 24*60*60],[x0 xc0 n0],options);
    x=Y(:,1); xc=Y(:,2); nrecep=Y(:,3); Thour=T/(60*60); Tday=T/(24*60*60);
%
    xfinal   = x(end);
    xcfinal  = xc(end);
    xchange  = x(end) - x0;
    xcchange = xc(end) - xc0;
    phase_initial = atan2(x0,xc0);
    phase_final   = atan2(xfinal,xcfinal);
    phase_change  = phase_final-phase_initial;
%
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % End of find_phase_response function
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
  dnrecepdt = Cparam.lambda*(alpha*(1-nrecep)-Cparam.betaval*nrecep);
%
  dy = zeros(3,1);
  dy(1) = dxdt;
  dy(2) = dxcdt;
  dy(3) = dnrecepdt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lightinput=lightinput(t,Lparam)
  lightinput=interp1(Lparam.time_secs,Lparam.light,t,'linear');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
