function res = residual(u, f, dr, dtheta, r) 
% Residuals over the grid (1/m, 1/n)?

[Mp2,Np1] = size(u);
M = Mp2 - 2;
N = Np1 - 1;

res = zeros(size(u)); 
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
    res(m, N+1) = res(m, 1); % Periodic
end
res(M+2, :) = res(M, :); % Insulating

end