function res = psi_residual(psi, f, dr, dtheta, r)
% Computes the residuals for the Poisson equation in an annulus 
% D^2 psi = f given delta r, delta theta, and vector of r values
% Editted code taken from Blackboard done by Dr Mestel

% Extract M and N
[Mp1,N] = size(psi);  
M = Mp1 - 1;
% Coefficients for diffusion equation
rr = 1 / dr;
rr2 = rr / dr;
rtheta = 1 / dtheta / dtheta;
% Array of interest
a0toNm1 = [N 1:N-1];
a2toNp1 = mod(1:N,N) + 1;

res = zeros(size(psi));   % Initalize matrix to store residuals
for m = 2:M
    r2dt2 = rtheta / r(m)^2;
    rdr = rr / r(m) / 2;
    sumdr2dt2 = 2 * (rr2 + r2dt2) + 1;
    % Compute values periodic in theta
    res(m, :) = f(m,:) - sumdr2dt2 * psi(m,:) + (rr2 + rdr) * psi(m+1, :) ...
            + (rr2 - rdr) * psi(m-1, :) + r2dt2 * ( psi(m, a2toNp1) + psi(m, a0toNm1) );
end

end
