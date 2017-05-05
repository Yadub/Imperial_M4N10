% Parameters (rPrandtl number)
Pr = 1;

% Domain and grid size
b =  2;
M = 2^6;                % Set M
N = 2^8;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 50;              % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
istep = 0;              % Number of time-steps 

% Initalize matrices 
ur = zeros(M+1,N);
utheta = zeros(M+1,N);
psi = zeros(M+1,N);
omega = zeros(M+1,N);
% Set omega
omega = annulusTinit(omega, dr, dtheta, 2.5);

F  = zeros(M+1,N);
% Set the forcing function
F = annulusTinit(F, dr, dtheta, 1);
% Dont need F to have values of 1 on the boundary
F(1,:) = 0;

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();
fig_stream = figure();
fig_omega = figure();
fig_F = figure();

% Plot the starting omega function
plotAnnulusScalarField( omega, xx, yy, fig_omega );
title('Q3: Starting omega');

% Plot the forcing function
plotAnnulusScalarField( F, xx, yy, fig_F );
title('Q3: Forcing function');

% Start figure

% Set residual
res = 1;
% Set tolerance
tol = 1e-10;
% Advect and then Diffuse forward in time
while res > tol
    % Save current omega for residual comaprison
    omegaOld = omega;
    % Advection step for vorticity
    omega = annulusAdvect(ur, utheta, omega, dt, dr , dtheta);   
    % Diffusion step for vorticity
    omega = omega_Diffuse( omega, F, psi, dt, dr , dtheta, r, Pr);
    % Solve Poisson equation for streamfunction, psi
    psi = psi_Eqn(psi, -omega, dr, dtheta, r);
    % Derive the new velocity field from the new stream function
    [ur, utheta] = annulusVelocity(psi, dr, dtheta, r);

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
    
    % Plot
    if ((mod(istep,ipic)==0)) 
       display(istep)
       plotAnnulusScalarField( omega, xx, yy, fig_omega );
       plotAnnulusScalarField(psi, xx, yy, fig_stream);
       plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);
    end
    
    % Compute the residual between iterations
    res = max(max(abs(omega - omegaOld)));
    % Iterate number of steps
    istep = istep + 1;
end

% Plot final vorticity
plotAnnulusScalarField(omega, xx, yy, fig_omega);
title('Q3: Final Vorticity')
% Plot final streamfunction
plotAnnulusScalarField(psi, xx, yy, fig_stream);
title('Q3: Final Stream Function')
% Plot the final velocity Field
plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);
title('Q3: Final Velocity Function')

