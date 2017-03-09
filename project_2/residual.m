function res = residual(u, f, dr, dtheta, r)
% Computes the residuals for the Poisson equation in an annulus 
% D^2 u = f given delta r, delta theta, and vector of r values
% Editted code taken from Blackboard done by Dr Mestel

[Mp2,Np1] = size(u);    % Get M and N
M = Mp2 - 2;
N = Np1 - 1;

res = zeros(size(u));   % Initalize matrix to store residuals
for m = 2:M+1
    for n = 1:N
        dr2 = 1 / dr^2;
        r2dt2 = 1 / r(m)^2 * 1 / dtheta^2;
        rdr = 1 / r(m) * 1 / (2 * dr);
        sumdr2dt2 = 2 * (dr2 + r2dt2);

        res(m, n) = f(m,n) + sumdr2dt2 * u(m,n) - ((dr2 + rdr) * u(m+1, n) ...
            + (dr2 - rdr) * u(m-1, n) + r2dt2 * (u(m, n+1) ...
            + u(m, mod(n-2, N)+1)));
    end
    res(m, N+1) = res(m, 1); % BC: Periodic (Dirichlet)
end
res(M+2, :) = res(M, :); % BC: Insulating (Neumann)

end
