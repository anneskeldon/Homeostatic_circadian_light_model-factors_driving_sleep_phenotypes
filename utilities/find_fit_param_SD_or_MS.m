function [x,fval,exitflag] = find_param(fitparam,param,Cparam,setup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses fzero to find either mu by fitting for sleep duration or tau by fitting for midsleep
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First finds a window in which there is a zero. (Uses the fact that mid-sleep and sleep duration
% vary monotonically with the fit parameters.
%
  if fitparam.fit_param_MS == 1
   objfun=@(param,Cparam,setup,x0)objfunc_meanMS(param,Cparam,setup,x0)
  end
  if fitparam.fit_param_MS == 0
   objfun=@(param,Cparam,setup,x0)objfunc_meanSD(param,Cparam,setup,x0)
  end
%
  maxloops = fitparam.maxloops;
  iloop    = 0;
%
  x0_bound1   = fitparam.x0_bound1;
  x0          = x0_bound1;
  [residual1] = objfun(param,Cparam,setup,x0);
  residual2   = residual1;
    if residual1 >0
      while (residual1*residual2 > 0) & iloop < maxloops
        x0_bound2   = x0_bound1 - fitparam.delta;
        x0          = x0_bound2;
        [residual2] = objfun(param,Cparam,setup,x0);
          if residual1*residual2 > 0
            x0_bound1 = x0_bound2;
          end
        iloop = iloop + 1
      end
    else
      while (residual1*residual2 > 0) & iloop < maxloops
        x0_bound2  = x0_bound1 + fitparam.delta;
        x0         = x0_bound2;
        [residual2] = objfun(param,Cparam,setup,x0);
            if residual1*residual2 > 0
              x0_bound1 = x0_bound2;
            end
        iloop = iloop + 1;
      end
    end
%
% If an interval has been found containing a root, then use fzero to get an accurate
% solution.
%
  if iloop < maxloops
    x0_guess = [x0_bound1 x0_bound2];
    options                  = optimset('TolX',1e-5,'TolFun',1e-5);
    [x,fval,exitflag,output] = fzero(@(x0) objfun(param,Cparam,setup,x0),x0_guess,options)
  else
    x        = x0_bound1;
    fval     = -1;
    exitflag = -3;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
