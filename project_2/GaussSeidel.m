function u = GaussSeidel( u0, f, dr, dtheta, r, iters  )


[Mp2, Np1] = size(u0);  % Extract size of the problem to be solved
M = Mp2 - 2;            % Set size of M
N = Np1 - 1;            % Set size of N

u = u0;

for i = 1:iters
    for m = 2:M+1
        for n = 1:N
            t1 = 1/dr^2;
            t2 = 1/r(m)^2 * 1/dtheta^2;
            t3 = 1/r(m) * 1/(2 * dr);
            t4 = 2 * (t1 + t2);
            u(m, n) = ((t1 + t3) * u(m+1, n) + (t1 - t3) * u(m-1, n) + ...
                t2 * (u(m, n+1) + u(m, mod(n-2, N)+1)) - f(m,n)) / t4;
        end
        u(m, N+1) = u(m, 1); % BC: Periodic
    end
    u(M+2, :) = u(M, :); % BC: Insulating
end

end