
% Domain and grid size
b =  5;
M = 2^5;                % Set M
N = 2^6;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
ipic = 50;              % How many steps between pictures
dt = .01*min(dr,dtheta);% Set timestep
tstep = 1000;           % Number of time-steps 

ur  = zeros(M+1,N);     % Velocity: r component
utheta  = zeros(M+1,N); % Velocity: theta component
Q  = zeros(M+1,N);      % Scalar field
Q(1,:) = 1;             % Start with values of 1 on r=1 and 0 everywhere else

r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array
% Uncomment to remove blob
rc = (b-1)/2;
for m = 2:M
    Q(m,:) = exp( - sqrt( (r(m) - rc)^2 + (theta(1:N) - pi).^2) );
end
% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);

% Save original Q for error comaprison
Q0 = Q;
% Start figure
fig_advect = figure();
% Advect forward in time
for istep=1:tstep
    % Plot
    if ((mod(istep,ipic)==0)|(istep==1)) 
       istep
       figure(fig_advect) 
       Tplot = zeros(M+1,N+1);
       Tplot(:,1:N) = Q(:,1:N);
       Tplot(:,N+1) = Q(:,1);
       contourf(xx',yy',Tplot,10,'r')
       axis image
       xlabel('x','FontSize',18)
       ylabel('y','FontSize',18);
       colorbar
       drawnow
       pause(1)
    end
    % Advection step for temperature  
    Q = annulusDiffuseT(Q, dt, dr, dtheta, r);
end
