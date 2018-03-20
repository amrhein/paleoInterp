% function [ym,yms] = get_ym(y,t,N)
% Least-squares estimate of the mean of a potentially unevenly-spaced time
% series computed using an estimate of the signal structure function
%
% D Amrhein, July 2013
% Edited for readability, September 2015
%
% % INPUTS
%
% y     vector series of values to be interpolated
% t     vector of times (possibly unevenly spaced) associated with y
% N     observational uncertainty covariance matrix having units of
%       sigma^2. In the case where observational uncertainty does not covary in
%       time and has a constant value sigma, N is eye(length(y))*sigma^2
%
% % OUTPUTS
%
% ym    least-squares estimate of the record mean
% yms   uncertainty of the mean estimate


function [ym,yms] = get_ym(y,t,N)

y = y(:);
t = t(:);

msy = mean(y.^2);
ly = length(y);
ynm = y-nanmean(y);

% compute the structure function estimate
[blag,brms,bvar,lagv,rmsv,lagm] = strufun(t,ynm,N,20);
bg2 = (blag>=0 & ~isnan(brms));
p = polyfit(log10(blag(bg2)),log10(brms(bg2)),1);
a = 10^p(2);
b = p(1);

strf = @(tau) a*tau.^(b);

% lagm is a lower triangular matrix. make it full with a diagonal of zeros:
lagmf = lagm+lagm';
lagmf(~~eye(size(lagmf))) = 0;

% generate the signal covariance function
S = mean(ynm.^2) - 0.5*strf(abs(lagmf));
S(abs(S)==inf)=0; % can happen for blue spectra
%S(S(:)<0)=0;

% compute an oft-used inverse
iSN = inv(S+N);

% mean estimator. See Rybicki and Press 1992 or Wunsch 2006 2.413
ee = ones(ly,1);
ym = ee'*(iSN*y)/(ee'*iSN*ee);
yms = 1/(ee'*iSN*ee);


