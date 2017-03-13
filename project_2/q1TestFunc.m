function [avgError, iters, timetaken, u] = q1TestFunc( b, M, N, tol, p)
% Gets the avgError, iterations, timetaken, and solution u for the 
% test function described in the report

if nargin < 4, tol = 1e-8; end  % Set tol
if nargin < 5, p = 2; end       % Set p    

dr = (b - 1) / M;       % Set delta
dtheta = 2 * pi / N;    % Set delta theta

r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array

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
tic;
[u, iters] = Jacobi( u0, S, r, dr, dtheta, tol );
timetaken = toc;
% Compute the error
Error = ( U(1:M+1, 1:N) - u(1:M+1, 1:N) );
avgError = mean(mean(abs(Error)));

end
