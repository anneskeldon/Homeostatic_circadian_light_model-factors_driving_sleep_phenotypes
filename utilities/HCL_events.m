function [value,isterminal,direction] = HCL_events(t,Y,Cparam,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specifies the conditions for switching between wake and sleep and sleep and wake for the HCL model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  H  = Y(1);
  dy = HCL_equations(t,Y,Cparam,param);
  C  = circadian(Y(2),Y(3),param);
%
  if param.onoff == 1              %if awake, stop at upper threshold
    H_threshold   = param.H0p + C;
    value(1)      = H-H_threshold;
    isterminal(1) = 1; % 1 stop the integration, 0 continue.
    direction(1)  = 0; % 0, any direction, 1 positive direction,  -1 only negative direction.
  else                           %if asleep, stop at lower threshold
    H_threshold   = param.H0m + C;  %if asleep, stop at lower threshold
    value(2)      = H-H_threshold;
    isterminal(2) = 1; % 1 stop the integration, 0 continue.
    direction(2)  = 0; % 0, any direction, 1 positive direction,  -1 only negative direction.
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
