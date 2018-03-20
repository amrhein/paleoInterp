function [blag,brms,bvar,lagv,rmsv,lagm] = strufun(t,s,N,nb)
% Computes the structure function for a time series following
% Press et al.1992, 'Time Delay of Gravitational Lens 0957 + 561 I., Eq. 22
% D Amrhein, September 2013
% Updated for readability September 2015
%
% % INPUTS
%
% s     vector series of values to be interpolated
% t     vector of times (possibly unevenly spaced) associated with y
% N     observational uncertainty covariance matrix having units of
%       sigma^2. In the case where observational uncertainty does not covary in
%       time and has a constant value sigma, N is eye(length(y))*sigma^2
% nb    number of bins used for binned estimates of lag, variance, and rms
%
% % OUTPUTS
%

% Compute lags and squared differences (less noise) for all pairs of
% measurements. Store these in matrices. If computed between all possible
% pairs, redundant calculations will be performed, so only compute the
% upper triangular part.
l = length(t);
rmsm = [];
lagm = [];
for ii = 2:l
    for jj = 1:ii
        lagm(ii,jj) = abs(t(ii)-t(jj));
        rmsm(ii,jj) = (s(ii)-s(jj))^2 - N(ii,ii) - N(jj,jj);
    end
end

% Set negative squared differences to NaN to eliminate them from subsequent
% calculations.
rmsm(rmsm<0) = nan;

% Vectorize nonzero values of rmsm and lagm
li = ~triu(ones(size(rmsm)));
rmsv = rmsm(li(:));
lagv = lagm(li(:));

% Values of the structure function will be binned in log-log space for
% computation of a power-law in lsinterp.m.
% compute bin edges
b_e = logspace(log10(min(lagv)),log10(max(lagv)),nb+1);
% initialize binned avg rms
brms = nan(nb,1); 
% initialize intra-bin variance of binned rms
bvar = brms;
% initialize bin centers
blag = b_e(1:end-1)+diff(b_e)/2;

% Populate bin means and variances
for ii = 1:nb
    brms(ii) = nanmean(rmsv(lagv>=b_e(ii) & lagv<b_e(ii+1)));
    bvar(ii) = nanvar(rmsv(lagv>b_e(ii) & lagv<b_e(ii+1)));
end
    
% Make these column vectors
brms = brms(:);
blag = blag(:);
bvar = bvar(:);
    
    
    
    

