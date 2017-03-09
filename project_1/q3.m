a = 1; % Set a
b = 2; % Set b
D = 0.01; % Diffusivity constant
Q = 0.01;
tmax = 4;

% Set h
N = 100;
h = (b - a) / N; % h = 0.01

% Set k
maxdt = 800; % Works for this value
% maxdt = 780; % Breaks for this value
k = tmax / maxdt; % k = 0.005 or just over 0.005

[u, r, t] = solveq2(a,b,Q,D,N,maxdt,k); % Solve the equation

% Set point
pr = N/2 + 1;
r0 = r(pr);
t0 = t(end);
u0 = (b-r0)/(b-a) + (r0-b)*(r0-a)*exp(-t0);

E1 = abs(u0 - u(pr,end)) % Error is printed to compare