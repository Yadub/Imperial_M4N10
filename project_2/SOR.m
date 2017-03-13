function [u, iters, res] = SOR(u0, f, w, r, dr, dtheta, tol)
% This function takes the intial u0 grid for an annulus and computes the 
% solution using Jacobi iteration with source function S

[Mp2, Np1] = size(u0);  % Extract size of the problem to be solved
M = Mp2 - 2;            % Set size of M
N = Np1 - 1;            % Set size of N

u = u0;                 % Initalize another u to store updated values

iters = 0;
res = 1;
while res > tol
    % Calculate new iterates of u
    for m = 2:M+1
        dr2 = 1 / dr^2;
        r2dt2 = 1 / r(m)^2 * 1 / dtheta^2;
        rdr = 1 / r(m) * 1 / ( 2 * dr );
        sumdr2dt2 = 2 * (dr2 + r2dt2);
        for n = 1:N
            u(m, n) = (1 - w) * u(m, n) + ((dr2 + rdr) * u(m+1, n) ...
                        + (dr2 - rdr) * u(m-1, n) + r2dt2 * (u(m, n+1) ...
                        + u(m, mod(n-2, N)+1)) - f(m, n)) * w / sumdr2dt2;
        end
        u(m,N+1) = u(m, 1);           % BC: Periodic
    end
    u(M+2,:) = u(M,:);                % BC: Insulating (Neumann)
    res = max(max(abs(u - u0))); % Compute the residual
%     res = max(max(abs(residual( u, f, dr, dtheta, r)))); % For Q2
    u0 = u;                           % Update u
    iters = iters + 1;                % Increase iteration count
end

end
