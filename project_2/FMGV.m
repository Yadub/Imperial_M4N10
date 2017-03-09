function [u, error] = FMGV(u0, NV, S, b)
% NV multigrid cycles for u0 initial conditions with 
% u(1,theta) = f(theta) boundary conditions

[Mp2,Np1] = size(u0);
M = Mp2 - 2;
N = Np1 - 1;

dr = (b - 1) / M;     % Set delta
dtheta = 2 * pi / N;  % Set delta theta
r = 1:dr:b;           % Initialize r array

if nargin < 4
    S = zeros(Mp2,Np1);
end

%   determine maximum number of grid levels
klevel = round(log(min(M,N)) / log(2))

% initialization for when full MG is implemeted
% u0 = FullMG(f);

% Now do Multigrid V-cycles
for k=1:NV
    u0 = MultiGridV( u0, S, dr, dtheta, r);
end

% and tell us the answer
error = max(max(residual(u0, S, dr, dtheta, r)));
u = u0;
  
