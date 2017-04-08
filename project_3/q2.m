% Domain and grid size
k = .1;
b =  5;
M = 2^6;                % Set M
N = 2^7;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 50;              % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
tstep = 10000;           % Number of time-steps 

% Some non uniform vector field that is 0 on the boundaries
[ur, utheta] = annulusUinit( M, N, dr, dtheta, r);

T  = zeros(M+1,N);      % Scalar field
% Start with some non uniform scalar field
T = annulusTinit(T, dr, dtheta, 2);

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();

% Plot velocity field
plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);

% Save original T for error comaprison
T0 = T;
% Start figure
fig_advect = figure();
% Advect forward in time
for istep=1:tstep
    % Plot
    if ((mod(istep,ipic)==0)||(istep==1)) 
       display(istep)
       plotAnnulusScalarField(T, xx, yy, fig_advect);
    end
    % Advection step for temperature  
    T = annulusAdvect(ur, utheta, T, dt, dr , dtheta);
    % Diffusion step for temperature  
    T = annulusDiffuseT(T, dt, dr, dtheta, r, k);
end
