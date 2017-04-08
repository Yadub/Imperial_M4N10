% Parameters (Rayleigh number, Prandtl number, kappa)
Ra = 1000;
Pr = .7;

% Domain and grid size
b =  5;
M = 2^5;                % Set M
N = 2^6;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 50;              % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
tstep = 1000;           % Number of time-steps 

% Initalize matrices 
ur = zeros(M+1,N);
utheta = zeros(M+1,N);
psi = zeros(M+1,N);
omega = zeros(M+1,N);
td  = zeros(M+1,N);

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();

% Save old bouancy term
omegaRhs = td;
% Start figure
fig_advect = figure();
% Advect forward in time
for istep=1:tstep
    % Advection step for vorticity
    omega = annulusAdvect(ur, utheta, omega, dt, dr , dtheta);   
    % Diffusion step for vorticity
    omega = annulusDiffuseOmega( omega, omegaRhs, psi, dt, dr , dtheta, r, Pr);
    % Solve Poisson equation for streamfunction, psi
    psi = annulusPsiEqn(psi, -omega, dr, dtheta, r);
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
    if ((mod(istep,ipic)==0)||(istep==1)) 
       display(istep)
       plotAnnulusScalarField(psi, xx, yy, fig_advect);
       plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);
    end
    
end

% Plot final streamfunction
fig_stream = figure(); 
plotAnnulusScalarField(psi, xx, yy, fig_stream);

