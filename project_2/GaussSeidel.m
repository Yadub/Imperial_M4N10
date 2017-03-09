function u = GaussSeidel( u0, f, dr, dtheta, r, iters  )
% Computes an estimate u given initial guess u0 for the Poisson Equation
% in an annulus defined by delta r, delta theta, and the vector of r values
% Smooths for #iters

[Mp2, Np1] = size(u0);  % Extract size of the problem to be solved
M = Mp2 - 2;            % Set size of M
N = Np1 - 1;            % Set size of N

u = u0;
for i = 1:iters
    for m = 2:M+1
        for n = 1:N
            dr2 = 1 / dr^2;
            dr2dt2 = 1/r(m)^2 * 1/dtheta^2;
            rdr = 1/r(m) * 1/(2 * dr);
            sumdr2dt2 = 2 * (dr2 + dr2dt2);
            
            u(m, n) = ((dr2 + rdr) * u(m+1, n) + (dr2 - rdr) * u(m-1, n) + ...
                dr2dt2 * (u(m, n+1) + u(m, mod(n-2, N)+1)) - f(m,n)) / sumdr2dt2;
        end
        u(m, N+1) = u(m, 1); % BC: Periodic (Dirichlet)
    end
    u(M+2, :) = u(M, :); % BC: Insulating (Neumann)
end

end