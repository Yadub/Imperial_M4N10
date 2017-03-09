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

tol = 1e-4; % Tolerance level for the residual
largest_residual = 1; 
NV = 0;
tic; % Compute the solution and time taken
[u0, largest_residual] = FullMG(S, b, u0(1,:));  % Multigrif for inital loop
while largest_residual > tol                    % V loops until residual small
    [u, largest_residual] = FMGV( u0, S, dr, dtheta, r, 1);
    u0 = u;
    NV = NV + 1;
end
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

fig_3D = figure();
surf(x,y,u(1:M+1,:)')
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right

% Method accuracy
Errors = abs(U-u(1:M+1, 1:N+1));
largest_difference = max(max(Errors))
display(largest_residual)
display(NV)
