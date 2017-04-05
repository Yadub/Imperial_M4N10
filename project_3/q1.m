
% Domain and grid size
b =  2;
M = 2^5;                % Set M
N = 2^6;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
ipic = 1000;            % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
tstep = 20000;          % Number of time-steps 

ur  = zeros(M+1,N);     % Velocity: r component
utheta  = zeros(M+1,N); % Velocity: theta component
Q  = zeros(M+1,N);      % Scalar field

% Start with some non uniform scalar field
Q = annulusTinit(Q, dr, dtheta);

% Initalize velocity field
[ur,utheta] = annulusUinit(ur, utheta, dr, dtheta);

% For plotting
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
fig_ufield = figure();

% Plot velocity field
plotAnnulusVectorField(ur, utheta, r, theta, fig_ufield);

% Save original Q for error comaprison
Q0 = Q;
% Start figure
fig_advect = figure();
% Advect forward in time
for istep=1:tstep
    % Plot
    if ((mod(istep,ipic)==0)|(istep==1)) 
       figure(fig_advect) 
       Tplot=zeros(M+1,N+1);
       Tplot(:,1:N) = Q(:,1:N);
       Tplot(:,N+1) = Q(:,1);
       contourf(xx',yy',Tplot,10,'r')
       axis image
       xlabel('x','FontSize',18)
       ylabel('y','FontSize',18);
       colorbar
       drawnow
    end
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, dt, dr , dtheta);

end
% Advect backwards in time
for istep=1:tstep
    % Plot
    if ((mod(istep,ipic)==0)|(istep==1)) 
       figure(fig_advect) 
       Tplot=zeros(M+1,N+1);
       Tplot(:,1:N) = Q(:,1:N);
       Tplot(:,N+1) = Q(:,1);
       contourf(xx',yy',Tplot,10,'r')
       axis image
       xlabel('x','FontSize',18)
       ylabel('y','FontSize',18);
       colorbar
       drawnow
    end
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, -dt, dr , dtheta);
end
% Check error
Diff = Q0 - Q;
Errors = abs(Diff);
residual = max(max(Diff))
total_error = norm(Errors) / ((M+1) * N)