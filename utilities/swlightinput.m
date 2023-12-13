function swlightinput = swlightinput(t,Cparam,SWparam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian 
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%       
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluates value of light for time t.  Expecting t in seconds
%
% Cparam.lightinput = 0 % Uses mathematical function consisting sum of tanh's as introduced in
%                       % Skeldon et al Sci Rep 7:45158 (2017) https://doi.org/10.1038/srep45158
% Cparam.lightinput = 1 % Uses data in Cparam.timedata, Cparam.lightdata, uses linear interpolation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Light data from a function. 
%
 if Cparam.lightinput == 0
    eveninglight = Cparam.eveninglight;
    daylight     = Cparam.daylight;
    dawn         = Cparam.dawn;
    dusk         = Cparam.dusk;
    steepness    = 1/6000;
    a            = eveninglight;
    b            = (daylight-eveninglight)/2;
    lightf       = @(a,b,dawn,dusk,t)a+b*(tanh(steepness*(t-dawn)) ...
        -tanh(steepness*(t-dusk)));
    swlightinput = lightf(a,b,dawn,dusk,mod(t,24*60*60));
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Light data read interpolated from Cparam.lightdata
  if Cparam.lightinput == 1
    swlightinput=interp1(Cparam.timedata,Cparam.lightdata,t,'linear');  
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
