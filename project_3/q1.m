
% Domain and grid size
b =  2;
M = 2^6;                % Set M
N = 2^7;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 1000;            % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
tstep = 20000;          % Number of time-steps 

Q  = zeros(M+1,N);      % Scalar field

% Start with some non uniform scalar field
Q = annulusTinit(Q, dr, dtheta);

% Initalize velocity field
[ur, utheta] = annulusUinit( M, N, dr, dtheta, r);

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();

% Plot velocity field
plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_ufield);

% Save original Q for error comaprison
Q0 = Q;
% Start figure
fig_advect = figure();
% Advect forward in time
for istep=1:tstep
    % Plot every few time steps
    if ((mod(istep,ipic)==0)||(istep==1)) 
       display(istep)
       plotAnnulusScalarField(Q, xx, yy, fig_advect);
    end
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, dt, dr , dtheta);

end
% Advect backwards in time
for istep=1:tstep
    % Plot every few time steps
    if ((mod(istep,ipic)==0)||(istep==1)) 
       display(istep)
       plotAnnulusScalarField(Q, xx, yy, fig_advect);
    end
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, -dt, dr , dtheta);
end
% Check error
Diff = Q0 - Q;
Errors = abs(Diff);
residual = max(max(Diff))
total_error = norm(Errors) / ((M+1) * N)

% TODO COMPUTE ERROR AS FUNCTION OF NUM TIME STEPS!