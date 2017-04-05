function [ Qnew ] = annulusAdvect( ur, utheta, Q, dt, dr, dtheta )
%ANNULUSADVECT: Advection routine in an annulus using Richtmyer's 
%               second-order Lax-Wendroff routine
% INPUTS:   ur,     r components of u
%           utheta, theta components of u
%           Q,      scalar field Q to be advected
%           dt,     time step
%           dr,     grid step size in the r direction
%           dtheta, grid step size in the theta direction
% OUTPUTS:  Qnew,   advected scalar field Q

u = ur;
v = utheta;

[Mp1,N] = size(ur);   % Extract size of grid
M = Mp1 - 1;
% Set derivative fraction
sr = dt/dr; 
stheta = dt/dtheta;

% Setup r integer and intermediate point matrix
ar = 1 + dr * (0:M);
arp = 1 + dr * ((0:M) + .5);
r = zeros(M+1,N);
rrp = zeros(M+1,N);
for i = 1:N
    r(:,i) = ar;
    rrp(:,i) = arp;
end

% Arrays of interest
a2toNp1 = mod(1:N,N) + 1;
aNtoNm1 = [N 1:N-1];
% Initalize arrays used for storage
Qtp = zeros(M+1,N);
Qnew = zeros(M+1,N);

% Find Qpp,upp,vpp at integer + half grid points by local averaging
Qpp = .25 * ( Q(1:M,1:N) + Q(1:M,a2toNp1) ...
                       + Q(2:M+1,1:N) + Q(2:M+1,a2toNp1) );
upp = .25 * ( u(1:M,1:N) + u(2:M+1,1:N) + u(1:M,a2toNp1) + u(2:M+1,a2toNp1) );
vpp = .25 * ( v(1:M,1:N) + v(2:M+1,1:N) + v(1:M,a2toNp1) + v(2:M+1,a2toNp1) );

% set up the velocities at the half-points - periodic in theta
urp = ( u(1:M,:) + u(2:M+1,:) ) / 2;
vtp = ( v(:,a2toNp1) + v(:,1:N) ) / 2; % Periodic in theta

% Define fluxes at integer points
F = Q .* u;
G = Q .* v ./ r;

% Define fluxes at intermediate points
Fpp = Qpp .* upp;
Gpp = Qpp .* vpp ./ rrp(1:M,:);

% First half-timestep for r (periodic in theta)
Qrp = .25 * ( Qpp + Qpp(:,aNtoNm1) + Q(2:M+1,:) + Q(1:M,:) ) ...
     - .5 * ( sr * ( F(2:M+1,:) - F(1:M,:) ) + stheta * ( Gpp - Gpp(:,aNtoNm1) )  );
% First half-timestep for theta (periodic in theta with boundary conditions for r below)
Qtp(2:M,:) = .25 * ( Qpp(2:M,1:N) + Qpp(1:M-1,1:N) + Q(2:M,a2toNp1) + Q(2:M,1:N) ) ...
     - .5 * ( sr * ( Fpp(2:M,1:N) - Fpp(1:M-1,1:N) ) + stheta * ( G(2:M,a2toNp1) - G(2:M,1:N) ) );
% Use condition that v=0=G on these walls
Qtp(1,:) = ( Qpp(1,1:N) + Q(1,a2toNp1) + Q(1,1:N) ) / 3 ...
     - .5 * ( sr * 2 * Fpp(1,1:N) + stheta * ( G(1,a2toNp1) - G(1,1:N) ) );
Qtp(M+1,:) = ( Qpp(M,1:N) + Q(M+1,a2toNp1) + Q(M+1,1:N) ) / 3 ...
     - .5 * ( sr * 2 * ( - Fpp(M,1:N) ) + stheta * ( G(M+1,a2toNp1) - G(M+1,1:N) ) ); 

% Define fluxes at j+1/2 points
Frp = Qrp .* urp;
Gtp = Qtp .* vtp ./ rrp;

% Finally, 2nd half time-step with periodic theta
Qnew(2:M,1:N) = Q(2:M,:) - sr * ( Frp(2:M,:) - Frp(1:M-1,:) ) ...
                           - stheta * ( Gtp(2:M,:) - Gtp(2:M,aNtoNm1) );

% Now y-edges G=0 on boundary: Needed???
% Qnew(1,1:N) = Q(1,:) - sr * ( Frp(1,:) ) - stheta * ( Gtp(1,:) - Gtp(1,aNtoNm1) );
% Qnew(M+1,1:N) = Q(M+1,:) - sr * ( - Frp(M,:) ) - stheta * ( Gtp(M+1,:) - Gtp(M+1,aNtoNm1) );

% Don't change r = 1 and r = b boundary.
Qnew(1,:) = Q(1,:);
Qnew(M+1,:) = Q(M+1,:);

end