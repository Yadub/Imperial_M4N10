% Domain and grid size
b =  2;
M = 2^6;                % Set M
N = 2^8;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
ipic = 1000;            % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
tstep = 10000;          % Number of time-steps 

Q  = zeros(M+1,N);      % Scalar field

% Start with some non uniform scalar field
Q = annulusTinit(Q, dr, dtheta);

% Initalize velocity field
[ur, utheta, psi] = annulusUinit( M, N, dr, dtheta, r);

% Determine the maximum velocity (for CFL condition)
urMax = max(max(abs(ur(2:M,:))));
uthetaMax = max(max(abs(utheta(2:M,:))));
velmax = max(urMax,uthetaMax);
Courant = velmax * dt / min(dr,dtheta);
dtCFL = .8 * min(dr,dtheta) / max(1,velmax);    
% Print courant number when it becomes smaller than dt
if (dt > dtCFL)
    fprintf('Dangerous Courant number!\n');
    fprintf('dt: %i, ', dt);
    fprintf('Courant number: %i\n', Courant);
end

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_u = figure();
fig_psi = figure();

% Plot stream function
plotAnnulusScalarField(psi, xx, yy, fig_psi);
title('Q1: Stream function')

% Plot velocity field
plotAnnulusVectorField(ur, utheta, R, Theta, xx, yy, fig_u);
title('Q1: Velocity Field')

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
       title(['Q1: Temperture Field for t = ', num2str(dt * (istep-1))])
       pause(0.1);
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
       title(['Q1: Temperture Field for t = ', num2str(dt * (tstep - (istep-1)) )]);
       pause(0.1)
    end
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, -dt, dr , dtheta);
end
title(['Q1: Temperture Field advected forward to t = ', num2str(dt * (istep-1)), ' and back to t = 0']);
% Check error
Diff = Q0 - Q;
Errors = abs(Diff(2:M,:));
max_error = max(max(Diff))
average_error = norm(Errors) / ((M-1) * N)
