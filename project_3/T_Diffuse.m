function [ T ] = T_Diffuse( T, dt, dr, dtheta, r, k )
%ANNULUSDIFFUSET: Pure diffusion of the temperature by taking a 
%                   Crank Nicolson step in the annulus

if nargin < 6, k = 1; end

rhs = T + .5* dt * k * T_delsqr(T, dr, dtheta, r);

tol = 1e-10;
maxres = 1;
while maxres > tol
    T = T_MultiGridV( T, rhs, dt, dr, dtheta, r, k);
    res = T_residual( T, rhs, dt, dr, dtheta, r, k);
    maxres = max(max(abs(res)));
end

end

