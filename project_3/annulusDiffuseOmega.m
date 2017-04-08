function [ omega ] = annulusDiffuseOmega( omega, omegaRhs, psi, dt, dr , dtheta, r, Pr )
%ANNULUSDIFFUSEOMEGA: Pure diffusion of the voriticity by taking a 
%                       Crank Nicolson step in the annulus

% Extract M and N
[Mp1,~] = size(omega);  
M = Mp1 - 1;
% Add buoyancy term (omrhs)
rhs = omega + dt * omegaRhs + .5 * dt * Pr * T_delsqr( omega, dr, dtheta, r);  

tol = 1e-10;
maxres = 1;
while maxres > tol
    omega = omega_MultiGridV( omega, rhs, psi, dt, dr, dtheta, r, Pr, M);
    res = omega_residual( omega, rhs, dt, dr, dtheta, r, Pr);
    maxres = max(max(abs(res)));
end

end

