function dy = HCL_equations(t,Y,Cparam,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differential equations for the HCL equations. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  H(1)      = Y(1); % Homeostatic sleep pressure
  xpace(1)  = Y(2); % x variable of the van der Pol oscillator model for the pacemaker
  ypace(1)  = Y(3); % y variable of the van der Pol oscillator model for the pacemaker
  nrecep(1) = Y(4); % fraction of activated receptors
%
  if param.onoff==0  %sleep
    chi = param.Xs;
  end
%
  if param.onoff==1  %wake
    chi=param.Xw;
  end
%
% Circadian functions
  if Cparam.gated == 1
    Itilde = param.onoff*swlightinput(t,Cparam,param);
  else
    Itilde = swlightinput(t,Cparam,param);
  end
%
  alpha  = Cparam.alpha0*(Itilde/Cparam.I0)^Cparam.p;
  Bhat   = Cparam.G*alpha*(1-nrecep);
  B      = (1-Cparam.b*xpace)*(1-Cparam.b*ypace)*Bhat;
%
% Differential equations
%
  dHdt      = (-H + param.onoff*param.mu)/chi;
%
  dxpacedt  = ( Cparam.gam*(xpace-4*xpace^3/3)-...
    ypace*((24/(Cparam.f*Cparam.tau_c))^2+Cparam.k*B) )/Cparam.kappa ;
  dypacedt  = ( xpace + B )/Cparam.kappa;
  dnrecepdt = Cparam.lambda*(alpha*(1-nrecep)-Cparam.betaval*nrecep);
%
  dy    = zeros(4,1);
  dy(1) = dHdt;
  dy(2) = dxpacedt;
  dy(3) = dypacedt;
  dy(4) = dnrecepdt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
