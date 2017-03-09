function [u, largest_residual] = FMGV(u0, f, dr, dtheta, r, NV)
% Performs NV multigrid V cycle for the Poisson equation in an annulus,
% D^2 u = f, for u0 initial conditions with u(1,theta) = f(theta)
% and given delta r, delta theta, and a vector of r values
% Editted code taken from Blackboard done by Dr Mestel

if nargin < 6   % Incase number of V cycles hasn't been specified
    NV = 1;
end

[Mp2,Np1] = size(u0);   % Extract M and N
M = Mp2 - 2;
N = Np1 - 1;

% Do Multigrid V-cycles
for k=1:NV
    u0 = MultiGridV( u0, f, dr, dtheta, r);
end

% and tell us the answer
largest_residual = max(max(residual(u0, f, dr, dtheta, r)));
u = u0;

end