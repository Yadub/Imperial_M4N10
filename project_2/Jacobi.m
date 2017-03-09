function [u, iters] = Jacobi(u0, r, dr, dtheta, S, tol)
% This function takes the intial u0 grid for an annulus and computes the 
% solution using Jacobi iteration with source function S

[Mp2, Np1] = size(u0);  % Extract size of the problem to be solved
M = Mp2 - 2;            % Set size of M
N = Np1 - 1;            % Set size of N

u = u0;                 % Initalize another u to store updated values

iters = 0;
residual = 1;
while residual > tol
    for m = 2:M+1
        for n = 1:N
            dr2 = 1 / dr^2;
            r2dt2 = 1 / r(m)^2 * 1 / dtheta^2;
            rdr = 1 / r(m) * 1 / ( 2 * dr );
            sumdr2dt2 = 2 * (dr2 + r2dt2);
            
            u(m, n) = ( (dr2 + rdr) * u0( m + 1, n) + (dr2 - rdr) ...
                          * u0( m - 1, n) + r2dt2 * ( u0( m, n + 1 ) + ...
                          u0(m, mod( n - 2 , N ) + 1 )) + S( m, n)) / sumdr2dt2;
        end
        u(m, N+1) = u(m, 1);          % BC: Periodic (Dirichlet)
    end
    residual = max(max(abs(u - u0))); % Compute the residual
    u(M+2, :) = u(M, :);              % BC: Insulating (Neumann)
    u0 = u;                           % Update u
    iters = iters + 1;                % Increase iteration count
end

end
