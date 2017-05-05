% Domain and grid size
b =  2;
M = 2^6;                % Set M
N = 2^8;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
ipic = 100;              % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
istep = 0;           % Number of time-steps 

Q  = zeros(M+1,N);      % Initalize Q
Q = annulusTinit(Q, dr, dtheta, 2);

% For plotting
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);

% Start figure
fig_advect = figure();
% Set residual
res = 1;
% Set tolerance
tol = 1e-5;
% Diffuse forward in time
while res > tol
    % Save current Q for residual comaprison
    Qold = Q;
    % Plot every ipic steps
    if ((mod(istep,ipic)==0))
       plotAnnulusScalarField(Q, xx, yy, fig_advect)
       title(['Q2: Temperture Field for dt = ', num2str(dt), ', tstep = ', num2str(istep) ]);
       drawnow
       pause(0.1);
    end
    % Diffusion step for temperature  
    Q = T_Diffuse(Q, dt, dr, dtheta, r);
    % Compute the residual between iterations
    res = max(max(abs(Q - Qold)));
    % iterate istep
    istep = istep + 1;
end
