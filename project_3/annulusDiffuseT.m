function [ T ] = annulusDiffuseT( T, dt, dr, dtheta, r )
%ANNULUSDIFFUSET: Pure diffusion of the temperature by taking a 
%                   Crank Nicolson step in the annulus

rhs = T + dt * T_delsqr(T, dr, dtheta, r) / 2;

tol = 1e-10;
maxres = 1;
while maxres > tol
    T = T_MultiGridV( T, rhs, dt, dr, dtheta, r);
    res = T_residual( T, rhs, dt, dr, dtheta, r);
    maxres = max(max(abs(res)));
end

end

