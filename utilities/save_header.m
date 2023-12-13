function save_header(param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Writes out information at the top of a file specified in fid4 listing values of all parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = setup.fid4;
%
  fprintf(fid,'%7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n', ...
    [param.a, param.mu, param.H0p, param.H0m, param.Xs, param.Xw; ...
     Cparam.alpha0, Cparam.I0, Cparam.p, Cparam.G, Cparam.b, Cparam.gam;...
     Cparam.kappa, Cparam.f, Cparam.lambda, Cparam.betaval, Cparam.tau_c, 0; ...
     Cparam.daylight, Cparam.eveninglight, Cparam.dawn, Cparam.dusk, Cparam.gated, 0; ...
     0, 0, 0, 0, 0, 0; ...
    ]')
%
end
