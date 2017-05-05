% Parameters (Rayleigh number, Prandtl number, kappa)
Ra = 4000;
Pr = 2;
k = 1;

% Domain and grid size
b =  2;
M = 2^5;                % Set M
N = 2^7;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 50;              % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
istep = 1;          % Number of time-steps 

% Initalize matrices 
ur = zeros(M+1,N);
utheta = zeros(M+1,N);
psi = zeros(M+1,N);
omega = zeros(M+1,N);
T  = zeros(M+1,N);
% Seed a small perturbation at r = 1
T = 0.01 * annulusTinit(T, dr, dtheta, 2.5);
% Set T on the boundary r=1 to be 1
T(1,:) = 1;

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();

% Start figure
fig_advect = figure();
% Set residual
res = 1;
% Set tolerance
tol = 1e-5;
% Advect and then Diffuse forward in time as well compute vorticity
while res > tol
    % Save current omega for residual comaprison
    Told = T;
    
    % Advection step for temperature  
    T = annulusAdvect(ur, utheta, T, dt, dr , dtheta);
    % Diffusion step for temperature  
    T = T_Diffuse(T, dt, dr, dtheta, r, k);
    
    % Advection step for vorticity
    omega = annulusAdvect(ur, utheta, omega, dt, dr , dtheta);
    % Compute the bouancy term use two step Adams-Bashforth for 2nd order
    td = annulusBuoyancy(T, dr, dtheta, r, Pr, Ra);
    % Needs kickstarting on first step
    if(istep==1)
        omegaRhs = td;
    else
        omegaRhs = 1.5 * td - .5 * tdold;
    end
    % Save old bouancy term
    tdold = td;
    % Diffusion step for vorticity
    omega = omega_Diffuse( omega, omegaRhs, psi, dt, dr , dtheta, r, Pr);
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
        fprintf('Dangerous Courant number, %i\n', Courant);
    end
    
    % Plot Temperature Field and Velocity field
    if ((mod(istep,ipic)==0)) 
       display(istep)
       plotAnnulusScalarField(T, xx, yy, fig_advect);
       plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);
    end
    
    % Compute the residual between iterations
    res = max(max(abs(T - Told)));
    % iterate istep
    istep = istep + 1;
end

% Plot final T and velocity
plotAnnulusScalarField(T, xx, yy, fig_advect);
title(['Final Temperature. Ra = ',num2str(Ra), ', Pr = ', num2str(Pr), ', kappa = ', num2str(k)])
plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);
title(['Final Velocity Field. Ra = ',num2str(Ra), ', Pr = ', num2str(Pr), ', kappa = ', num2str(k)])
