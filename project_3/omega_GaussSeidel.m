function omega = omega_GaussSeidel( omega, f, psi, dt, dr, dtheta, r, iters, Pr, bc)
% Computes an estimate u given initial guess u0 for the Poisson Equation
% in an annulus defined by delta r, delta theta, and the vector of r values
% Smooths for #iters

% Extract M and N
[Mp1,N] = size(omega);  
M = Mp1 - 1;
% Coefficients for diffusion equation %%%%%%% Need Kappa???? %%%%%%%%%%
dtdr = Pr * dt / 2 / dr;
dr2 = dtdr / dr;
rtheta = Pr * dt / 2 / dtheta / dtheta;
% Array of interest
a0toNm1 = [N 1:N-1];
a2toNp1 = mod(1:N,N) + 1;
% Implement boundary condition if on the coarsest level
if bc
    omega = omega_BC(omega, psi, dr);
end
% Iterate
for i = 1:iters
    for m = 2:M % <- Not changing r=1 and r=b values
        dr2dt2 = rtheta / r(m)^2;
        rdr = dtdr / r(m) / 2;
        sumdr2dt2 = 2 * (dr2 + dr2dt2) + 1;
        % Compute values periodic in theta
        omega(m, :) = ((dr2 + rdr) * omega(m+1, :) + (dr2 - rdr) * omega(m-1, :) + ...
                dr2dt2 * (omega(m, a2toNp1) + omega(m, a0toNm1)) + f(m,:)) / sumdr2dt2;

    end
    % Implement boundary condition if on the coarsest level
    if bc
        omega = omega_BC(omega, psi, dr);
    end
end

end