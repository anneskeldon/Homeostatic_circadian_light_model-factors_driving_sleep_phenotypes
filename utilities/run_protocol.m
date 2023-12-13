function [Htot,xtot,ytot,nreceptot,Ttot,param] = runprotocol(param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run protocol
%
% Loops through each of the time windows specified in the protocol file either running
% in spontaneous mode if setup.transitions(transition_number,3) is equal to 0 and
% in overriding thresholds if setup.transitions(transition_number,3) is equal to 1. 
%
% This version only runs in SPONTANEOUS MODE. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initital conditions for the integrator
  Y0 = setup.Y0; %[H0 x0 y0 0];
%
  Htot  = []; xtot  = []; ytot  = []; nreceptot = []; Ttot  = [];
  te    = []; ie    = []; ye    = [];
  tetot = []; ietot = []; yetot = [];
%
% Number of time windows. 
  ntimes = length(setup.transitions(:,1));
%
  saveparam = [param.a param.mu];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for transition_number = 1:ntimes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T0                  = setup.transitions(transition_number,1);
    Tend                = setup.transitions(transition_number,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If unconstrained, integrate as normal, recording sleep/wake transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if setup.transitions(transition_number,3) == 0  % Spontaneous transitions
%
      wokewithalarm       = 0;
      selfselectedbedtime = 1;
      scheduledwake       = 0;
%
% If above the upper threshold at the start, switch to sleep
      H_threshold = param.H0p + circadian(Y0(2),Y0(3),param);
      if Y0(1) >= H_threshold
        param.onoff = 0; % Start on sleep
      end
%
      while T0<Tend 
        [T,Y,te,ye,ie] = integrate_equations(Cparam,param,setup,T0,Tend,Y0);
%
%  Append data for latest period with prior periods. Don't include the last data point
        if transition_number == ntimes
          Htot      = [Htot;Y(:,1)];
          xtot      = [xtot;Y(:,2)];
          ytot      = [ytot;Y(:,3)];
          nreceptot = [nreceptot;Y(:,4)];
          Ttot      = [Ttot;T];
        else
          Htot      = [Htot;Y(1:end-1,1)];
          xtot      = [xtot;Y(1:end-1,2)];
          ytot      = [ytot;Y(1:end-1,3)];
          nreceptot = [nreceptot;Y(1:end-1,4)];
          Ttot      = [Ttot;T(1:end-1)];
        end
%
        if(~isempty(ie) & T<Tend)
          ie(end)=param.onoff;
          tetot = [tetot;te]; yetot = [yetot;ye]; ietot = [ietot;ie];
        end
%
        writeoutinfo(setup.fid1,ie,te,ye,saveparam,scheduledwake,wokewithalarm,selfselectedbedtime);
        te = []; ie = []; ye = [];
%
% Set   up initial conditions for next period
%
        H0 = Y(end,1);
        x0 = Y(end,2);
        y0 = Y(end,3);
        n0 = Y(end,4);
        Y0 = [H0 x0 y0 n0];
        T0 = T(end,1);
        if T0 < Tend
           param.onoff=abs(param.onoff-1);
        end
      end  % End while loop
    end % End spontaneous transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If constrained to be awake, integrate H and C but do not allow switches to sleep.
    if setup.transitions(transition_number,3) == 1  % Forced wake
      wokewithalarm       = 0;
      selfselectedbedtime = 0;
      scheduledwake       = 1;
    end % End forced wake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  end % End loop through protocol segments (ntransitions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
end % End function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function writeoutinfo(fid1,ie,te,ye,saveparam,scheduledwake,wokewithalarm,selfselectedbedtime)
  for i = 1:length(te)
    savevector = [saveparam(1) saveparam(2) scheduledwake ie(i) te(i)/(60*60) ye(i,:) wokewithalarm selfselectedbedtime 0 0 0];
    fprintf(fid1,'%7.4f %7.4f %3.1f %3.1f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %3.1f %3.1f %3.1f \n',...
      savevector);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
