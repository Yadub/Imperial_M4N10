% Domain and grid size
k = 0.001;
b = 2;
M = 2^6;                % Set M
N = 2^7;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 100;             % How many steps between pictures
dt = 1e-3;              % Set timestep
istep = 0;              % Number of time-steps 

% Some non uniform vector field that is 0 on the boundaries
[ur, utheta, psi] = annulusUinit( M, N, dr, dtheta, r);

T  = zeros(M+1,N);      % Scalar field
% Start with some non uniform scalar field
T = annulusTinit(T, dr, dtheta, 2.5);

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();
fig_psi = figure();

% Plot stream function
plotAnnulusScalarField(psi, xx, yy, fig_psi);
title('Q2: Stream function')

% Plot velocity field
plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);
title('Q2: Velocity Field function')

% Determine the maximum velocity (for CFL condition)
urMax = max(max(abs(ur(2:M,:))));
uthetaMax = max(max(abs(utheta(2:M,:))));
velmax = max(urMax,uthetaMax);
Courant = velmax * dt / min(dr,dtheta);
dtCFL = .8 * min(dr,dtheta) / max(1,velmax);    
% Print courant number when it becomes smaller than dt
if (dt > dtCFL)
    fprintf('Dangerous Courant number! %i\n', Courant);
    fprintf('dt: %i, ', dt);
    fprintf('dtCFL: %i\n', dtCFL);
end

% Start figure
fig_advect = figure();
% Set residual
res = 1;
% Set tolerance
tol = 1e-5;
% Advect and then Diffuse forward in time
while res > tol
    % Save current Q for residual comaprison
    Told = T;
    % Plot
    if ((mod(istep,ipic)==0)) 
       display(istep)
       plotAnnulusScalarField(T, xx, yy, fig_advect);
    end
    % Advection step for temperature  
    T = annulusAdvect(ur, utheta, T, dt, dr , dtheta);
    % Diffusion step for temperature  
    T = T_Diffuse(T, dt, dr, dtheta, r, k);
    % Compute the residual between iterations
    res = max(max(abs(T - Told)));
    % iterate istep
    istep = istep + 1;
end
% Plot the final result
plotAnnulusScalarField(T, xx, yy, fig_advect);
title(['Q2: Temperature Field converged. Tsteps = ', num2str(istep), ', kappa = ', num2str(k) ]);
