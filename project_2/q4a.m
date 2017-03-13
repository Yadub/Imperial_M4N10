% Yadu Bhageria
% 00733164
b = 5;                 % Chosen value of b
M = 2^8;                % Set M
N = 2^8;                % Set N

dr = (b - 1) / M;       % Set delta
dtheta = 2 * pi / N;    % Set delta theta

r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array

u0 = zeros(M+2,N+1);    % Initialize initial u
u0(1,1:N) = sin( theta(1:N) ) .* exp( cos( theta(1:N) ) ); % BC: u(1,theta)
u0(M+2,:) = u0(M,:);    % BC: Insulation (Neumann)
u0(:,N+1) = u0(:,1);    % BC: Periodic

S = zeros(M+1,N+1);     % Initalize the forcing term S

tol = 1e-8; % Tolerance level for the residual
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

fig_contour = figure();
contour(x,y,u(1:M+1,:)');
title('Estimated Solution: Part a'); xlabel('x'); ylabel('y');
colorbar;

fig_3D = figure();
surf(x,y,u(1:M+1,:)')
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
