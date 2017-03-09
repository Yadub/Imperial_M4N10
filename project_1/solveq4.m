function [ u, r , t] = solveq4( Q,D )
%SOLVEQ4 find u_inf(r) for various values of Q and D

% Fix a and b
a = 1;
b = 17;

N = 2^9; % Set grid size

dr = (b - a) / N;       % Steps in r (h = 0.001)

dt = 0.01;                 % Steps in time (k)

u = zeros( N + 1, 1); % Array for u(r,t) initialised to zero
unew = zeros( N + 1, 1);
r = zeros(1,N+1);
V = zeros(1,N+1);

% Initalize u_0 and compute S(r,t)
for n = 1:N+1
    r(n) = a + (n-1) * dr;
    V(n) = (Q+D) / r(n);
end
u(1,1) = 1; % Inital boundary condition

t = 0;
diff = 1;
while diff > 1e-10  % Step in time until difference is small
    t = t + dt;
    unew(1) = 1;   % Boundary Condition
    for n = 2:N
        u_r = (u(n+1) - u(n-1) ) / (2*dr);
        u_rr = ( u(n+1)-2*u(n)+u(n-1) ) / dr^2;
        Vterm = V(n) * u_r;
        Dterm = D * u_rr;
        unew(n) = u(n) + dt * (Vterm + Dterm);
    end
    unew(N+1) = 0; % Boundary Condition
    diff = norm(u - unew); % Compute difference
    u = unew;   % Update new values
end

end

