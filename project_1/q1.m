% Explicit method to solve u_t = V(r) u_r + D (u_xx + u_r / r) + S(r,t) 
% Parameters: a < r < b, N : number of points
% D > 0 : diffusivity, V(r) = Q/r : radial velocity
% Q > 0, t > 0

% For Q1: S = 0

a = 10; % Set a
b = 20; % Set b
D = 0.01; % Diffusivity constant
Q = 0.01;
N = 2^10; % Number of points to split the grid in

maxdt = 2^12; % Maximum number of timesteps
tmax = 10;
Nt = maxdt; 
dt = tmax / Nt;

dr = (b - a) / N; % sice of steps in r
time = 0;

u = zeros( N + 1, maxdt + 1); % Array for u(r,t) initialised to zero
S = zeros( N + 1, maxdt + 1); % Array for S(r,t) initialised to zero
r = zeros(1,N+1);
V = zeros(1,N+1);

for n = 1:(N+1)        % Set up grid r(n) and initialise u(n,0)
    r(n) = a + (n-1) * dr;
    V(n) = (Q+D) / r(n);
    if ( abs( r(n) ) <= (b+a)/2 )
        u(n,1) = 1;
%     else
%         u(n,1) = 0;
    end
end

time(1)=0;
for j=1:maxdt    % Step in time
    time(j+1) = time(j) + dt;
    u(1,j+1) = 1;   % Boundary Condition
    for n = 2:N
        u_r = (u(n+1,j) - u(n-1,j) ) / (2*dr);
        u_rr = ( u(n+1,j)-2*u(n,j)+u(n-1,j) ) / dr^2;
        Vterm = V(n) * u_r;
        Dterm = D * u_rr ;
        u(n,j+1) = u(n,j) + dt * (Vterm + Dterm + S(n,j+1) );
    end;
    u(N+1,j+1) = 0; % Boundary Condition
end;


figure();    % Contour u(x,t)
contour(time,r,u,30);
xlabel('time'); ylabel('r'); colorbar;

figure();    % u(x) for various t in frame moving with Vmax
hold on;
for j=1:uint32(maxdt/10):maxdt
    plot(r , u(:,j),'DisplayName',['t = ', num2str( time(j) ) ]);
    xlim([a,b])
end;
plot(r , u(:,end), 'DisplayName',['t = ', num2str( time(end) ) ]);
xlabel('r'); ylabel('u(r,t)'); legend('show');
hold off;