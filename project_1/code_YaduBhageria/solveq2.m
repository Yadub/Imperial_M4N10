function [ u, r, time ] = solveq2(a, b, Q, D, N, Nt, k, display )
%SOLVEQ2 Function to solve the function u(r,t) defined in the document.
% a < r < b
% Nt time steps of size k
% grid of size N+1

if nargin < 8
    display = false;
else
    display = true;
end

dr = (b - a) / N;       % Steps in r (h)

dt = k;                 % Steps in time (k)
time = zeros(Nt+1,1);   % Set up the time vector.
for j=1:Nt
    time(j+1)=time(j)+dt;
end

% u = (b-r)/(b-a) + (r-b)*(r-a)*exp(-t)
% u_t = - (r-b)*(r-a)*exp(-t)
% u_r = - 1 / (b-a) + ( 2*r - a - b ) * exp(-t)
% u_rr = 2*exp(-t)
% S(r,t) = u_t - V(r) * u_r - D (u_rr + u_r / r);

u = zeros( N + 1, Nt + 1); % Array for u(r,t) initialised to zero
S = zeros( N + 1, Nt + 1); % Array for S(r,t) initialised to zero and set up
r = zeros(1,N+1);
V = zeros(1,N+1);

% Initalize u_0 and compute S(r,t)
for n = 1:N+1
    r(n) = a + (n-1) * dr;
    V(n) = (Q+D) / r(n);
    u(n,1) = (b - r(n) ) / (b - a) + (r(n)-b) * (r(n)-a);
    for j = 1:Nt+1
        t = time(j);
        u_t = - (r(n)-b) * (r(n)-a) * exp(-t) ;
        u_r = - 1 / (b-a) + ( 2*r(n) - a - b ) * exp(-t);
        u_rr = 2*exp(-t);
        S(n,j) = u_t - V(n) * u_r - D * u_rr;
    end
end

for j=1:Nt    % Step in time
    u(1,j+1) = 1;   % Boundary Condition
    for n = 2:N
        u_r = (u(n+1,j) - u(n-1,j) ) / (2*dr);
        u_rr = ( u(n+1,j)-2*u(n,j)+u(n-1,j) ) / dr^2;
        Vterm = V(n) * u_r;
        Dterm = D * u_rr;
        u(n,j+1) = u(n,j) + dt * (Vterm + Dterm + S(n,j));
    end;
    u(N+1,j+1) = 0; % Boundary Condition
end;

if display == true
    figure();    % Contour u(x,t)
    contour(time,r,u,30);
    xlabel('time'); ylabel('r'); colorbar;

    figure();    % u(x) for various t in frame moving with Vmax
    hold on
    for j=1:uint32(Nt/10):Nt
        plot(r , u(:,j), 'DisplayName',['t = ', num2str( time(j) ) ]);
        xlim([a,b])
    end;
    plot(r , u(:,end), 'DisplayName',['t = ', num2str( time(end) ) ]);
    xlabel('r'); ylabel('u(r,t)'); legend('show');
    hold off
end

end

