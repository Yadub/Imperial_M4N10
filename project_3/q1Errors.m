function [max_error, average_error] = q1Errors(M,N,dt,tstep)

% Domain and grid size
b =  2;
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array

Q  = zeros(M+1,N);      % Scalar field

% Start with some non uniform scalar field
Q = annulusTinit(Q, dr, dtheta);

% Initalize velocity field
[ur, utheta] = annulusUinit( M, N, dr, dtheta, r);

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

% Save original Q for error comaprison
Q0 = Q;
% Advect forward in time
for istep=1:tstep
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, dt, dr , dtheta);
end
% Advect backwards in time
for istep=1:tstep
    % Advection step for temperature  
    Q = annulusAdvect(ur, utheta, Q, -dt, dr , dtheta);
end
% Check error
Diff = Q0 - Q;
Errors = abs(Diff(2:M,:));
max_error = max(max(Diff));
average_error = norm(Errors) / ((M-1) * N);

end
