function res = T_residual(T, f, dt, dr, dtheta, r, k)
% Computes the residuals for the Poisson equation in an annulus 
% D^2 T = f given delta r, delta theta, and vector of r values
% Editted code taken from Blackboard done by Dr Mestel

% Extract M and N
[Mp1,N] = size(T);  
M = Mp1 - 1;
% Coefficients for diffusion equation %%%%%%% Need Kappa???? %%%%%%%%%%
rr = k * dt / 2 / dr;
rr2 = rr / dr;
rtheta = k * dt / 2 / dtheta / dtheta;
% Array of interest
a0toNm1 = [N 1:N-1];
a2toNp1 = mod(1:N,N) + 1;

res = zeros(size(T));   % Initalize matrix to store residuals
for m = 2:M
    r2dt2 = rtheta / r(m)^2;
    rdr = rr / r(m) / 2;
    sumdr2dt2 = 2 * (rr2 + r2dt2) + 1;
    % Compute values periodic in theta
    res(m, :) = f(m,:) - sumdr2dt2 * T(m,:) + (rr2 + rdr) * T(m+1, :) ...
            + (rr2 - rdr) * T(m-1, :) + r2dt2 * ( T(m, a2toNp1) + T(m, a0toNm1) );
end

end
