% test_lsinterp.m
% A script that generates synthetic time series to 
% demonstrate functionality of lsinterp.
% D Amrhein September 2015

clear

%% Generate a synthetic time series with variable observational spacing.
% Length of synthetic time series
L = 50;

% Examples:

% % 1. White noise with random uniformly distributed [0 1] time steps
% y = (randn(L,1));
% t = cumsum(rand(L,1))

% % 2. AR1 with random chi-squared time steps
% y = cumsum(randn(L,1));
% t = cumsum((randn(L,1)).^2);

% 3. AR1 with random uniformly distributed [0 1] time steps
y = cumsum(randn(L,1));
t = cumsum(rand(L,1));

% Observational noise covariance matrix
N = 0.2^2*eye(length(y)); % e.g., no error covariance between observations

%%
% Compute and plot interpolation

% the grid we're interpolating to
TSTEP = 0.1;
tp = min(t):TSTEP:max(t);

close all

[yi,P,S,ym,yms] = lsinterp(ymm,t,tp,N);
se = sqrt(diag(P));
hold on
ciplot(yi-se,yi+se,tp)
plot(tp,yi,'k')

% Overlay the mean-removed raw time series
ym = get_ym(y,t,N);
ymm = y-ym;
plot(t,ymm,'rd','markerfacecolor','r')
