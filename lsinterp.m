% function [yi,P,S,ym,yms] = lsinterp(y,t,tp,N)
%
% Objectively interpolates unevenly sampled data using algebra from 
% "Press et al. 1992, 'Time Delay of Gravitational Lens 0957 + 561 I." and
% Wunsch 2006, "Discrete Inverse and State Estimation Problems"
%
% D Amrhein, September 2013
% Updated for readability September 2015
%
% % INPUTS
%
% y     vector series of values to be interpolated
% t     vector of times (possibly unevenly spaced) associated with y
% tp    vector of times to which the data will be objectively interp'd
% N     observational uncertainty covariance matrix having units of
%       sigma^2. In the case where observational uncertainty does not covary in
%       time and has a constant value sigma, N is eye(length(y))*sigma^2
%
% % OUTPUTS
%
% yi    values at the interpolant times tp
% P     minimum uncertainty of the interpolated record
% S     estimate of the signal covariance computed by
%       fiting an exponent to the structure function
% ym    least-squares estimate of the record mean
% yms   uncertainty of the mean estimate

function [yi,P,S,ym,yms] = lsinterp(y,t,tp,N)

% set to be column vectors
y = y(:); t = t(:);

% compute least-squares means using get_ym
[ym,yms] = get_ym(y,t,N);

ly = length(y);

ymm = y - ym;
msy = mean(ymm.^2);

% compute the structure function estimate and fit a line to it in log-log
% space
[blag,brms,bvar,lagv,rmsv,lagm] = strufun(t,ymm,N,20);
bg2 = (blag>0 & ~isnan(brms));
p = polyfit(log10(blag(bg2)),log10(brms(bg2)),1);
a = 10^p(2);
b = p(1);

disp(['a = ' num2str(a) ', b = ' num2str(b)])
% NB: low values of b sometimes assoc'd with aesthetically unpleasing results

% The estimated power law structure function
strf = @(tau) a*tau.^(b);

% lagm is a lower triangular matrix. make it full with a diagonal of zeros:
lagmf = lagm+lagm';
lagmf(~~eye(size(lagmf))) = 0;

% generate the signal covariance function
S = msy - 0.5*strf(abs(lagmf));
%S(S(:)<0)=0;

SN = (S+N);

% Generate the correlation matrix Ss (S* in Rybicki and Press 1992 eqn 5)
Ss = [];

for ii = 1:ly
    Ss(ii,:) = msy - 0.5*strf(abs(t(ii) - tp));
end

% Compute the interpolated time series. Rename it to be a column vector.
yi = Ss'*(SN\ymm);
yi = yi(:);

% Generate an estimate of the solution covariance matrix at the
% interpolated points
lyp = length(tp);
Rxx = [];
for ii = 1:lyp
    Rxx(ii,:) = msy - 0.5*strf(abs(tp(ii) - tp));
end

P = Rxx - Ss'*(SN\Ss); % minimum uncertainty
