% Yadu Bhageria
% 00733164
b = 5;                 % Chosen value of b
tol = 1e-8;             % Set tolerance

M = 2^6;                % Set M
N = 2^7;                % Set N

dr = (b - 1) / M;       % Set delta
dtheta = 2 * pi / N;    % Set delta theta

r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array

p = 2;                  % Set value of p
u0 = zeros(M+2,N+1);    % Initialize initial u
u0(1,1:N) = (1 - b)^2 * sin( p * theta(1:N) ); % BC: u(1,theta)
u0(M+2,:) = u0(M,:);    % BC: Insulation (Neumann)
u0(:,N+1) = u0(:,1);    % BC: Periodic

% The exact solution
S = zeros(M+1,N+1);     % Initalize the forcing term S
U = zeros(M+1,N+1);     % Initialize exact solution of u
for m = 1:M+1
  for n = 1:N+1
    U(m,n) = ( (r(m) - b)^2 ) / r(m) * sin( p * theta(n) );
    S(m,n) = ( 2 * b * p^2 / r(m) + (1 - p^2) * (1 + b^2 / r(m)^2) ) ...
        * sin(p * theta(n)) / r(m);
  end
end

% Compute the solution
[u, iters] = Jacobi( u0, S, r, dr, dtheta, tol );

% Get x-y grid
[R, Theta] = meshgrid(r, theta);
x = R.*cos(Theta);
y = R.*sin(Theta);

% Plot solution
fig_main = figure();
subplot(1,2,1);
contour(x,y,U(1:M+1, 1:N+1)', 50);
title('Exact Solution'); xlabel('x'); ylabel('y');
subplot(1,2,2);
contour(x,y,u(1:M+1, 1:N+1)', 50);
title('Estimated Solution'); xlabel('x'); ylabel('y');

Error = (U(1:M+1, 1:N)-u(1:M+1, 1:N));
avgError = mean(mean(abs(Error)));

