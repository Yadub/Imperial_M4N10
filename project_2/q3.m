% Yadu Bhageria
% 00733164
b = 10;                 % Chosen value of b
M = 2^7;                % Set M
N = 2^7;                % Set N

dr = (b - 1) / M;       % Set delta
dtheta = 2 * pi / N;    % Set delta theta

r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array

p = 3;                  % Set value of p for f(theta)
u0 = zeros(M+2,N+1);    % Initialize initial u
u0(1,1:N) = (1 - b)^2 * sin( p * theta(1:N) ); % BC: u(1,theta)
u0(M+2,:) = u0(M,:);    % BC: Insulation (Neumann)
u0(:,N+1) = u0(:,1);    % BC: Periodic (Dirichlet)

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

tic; % Compute the solution and time taken
[u, largest_residual] = FMGV( u0, S, dr, dtheta, r, 2); % Source term is f
time_taken = toc

% Get x-y grid
[R, Theta] = meshgrid(r, theta);
x = R.*cos(Theta);
y = R.*sin(Theta);

% Plot solution
fig_main = figure();
subplot(1,2,1);
contour(x,y,U(1:M+1, 1:N+1)', 50);
title('Exact Solution');

subplot(1,2,2);
contour(x,y,u(1:M+1, 1:N+1)', 50);
title('Estimated Solution');

% Method accuracy
Errors = abs(U-u(1:M+1, 1:N+1));
display(largest_residual)
largest_difference = max(max(Errors))
