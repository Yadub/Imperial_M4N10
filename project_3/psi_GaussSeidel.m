function psi = psi_GaussSeidel( psi, f, dr, dtheta, r, iters)
% Computes an estimate u given initial guess u0 for the Poisson Equation
% in an annulus defined by delta r, delta theta, and the vector of r values
% Smooths for #iters

% Extract M and N
[Mp1,N] = size(psi);  
M = Mp1 - 1;
% Coefficients for diffusion equation %%%%%%% Need Kappa???? %%%%%%%%%%
dtdr = 1 / dr;
dr2 = dtdr / dr;
rtheta = 1 / dtheta / dtheta;
% Array of interest
a0toNm1 = [N 1:N-1];
a2toNp1 = mod(1:N,N) + 1;

for i = 1:iters
    for m = 2:M % <- Not changing r=1 and r=b values
        dr2dt2 = rtheta / r(m)^2;
        rdr = dtdr / r(m) / 2;
        sumdr2dt2 = 2 * (dr2 + dr2dt2) + 1;
        % Compute values periodic in theta
        psi(m, :) = ((dr2 + rdr) * psi(m+1, :) + (dr2 - rdr) * psi(m-1, :) + ...
                dr2dt2 * (psi(m, a2toNp1) + psi(m, a0toNm1)) + f(m,:)) / sumdr2dt2;

    end
end

end