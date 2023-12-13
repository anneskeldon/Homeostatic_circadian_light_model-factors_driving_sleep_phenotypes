function [TCmin,Cmin] = find_Cmin(T,timeseries)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Released as part of the codebase to replicate results reported in
% "Method to determine whether sleep phenotypes are driven by endogenous circadian
% rhythms or environmental light by combining longitudinal data and personalised mathematical models"
% Skeldon et al, PLoS Comput Biol, provisionally accepted Dec 2023.
%
% Author: A.C. Skeldon, a.skeldon@surrey.ac.uk, University of Surrey, 2023.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the maximum and minimum of values in the timeseries, one value for each day. 
% Assumes time T is in days and starts from 0. 
% If there are multiple points in the timeseries which have the same min or max value then the
% value at the mean of the indices is taken. Implicit in this is that identical min (or max) values
% occur sequentially. This is not foolproof. If results appear inconsistent / odd it is always a 
% good idea to cross-check by plotting the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ndays = floor(T(end));
%
  for n = 1:ndays
    start_time = n-1 - 2/24; % 
    end_time   = n   - 2/24; %
    Cvals      = timeseries(T>start_time & T<end_time);
    Tvals      = T(T>start_time & T<end_time);
    [max_Cvals,indmax] = max(Cvals);
    [min_Cvals,indmin] = min(Cvals);
    Tindmax = floor(mean(indmax));
    Tindmin = floor(mean(indmin));
    TCmax(n) = Tvals(Tindmax);
    TCmin(n) = Tvals(Tindmin);
    Cmax(n)  = max_Cvals;
    Cmin(n)  = min_Cvals;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
