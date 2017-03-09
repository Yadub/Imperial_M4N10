function Uout = MultiGridV(Uin, f, dr, dtheta, r)
% Performs 1 multigrid V cycle for the Poisson equation in an annulus,
% D^2 u = f, given a starting guess Uin for delta r, delta theta, vector r
% Editted code taken from Blackboard done by Dr Mestel

iters = 10; % Number of Gauss Seidel Smoothing iterations

[Mp2,Np1] = size(Uin);  % Extract M and N
M = Mp2 - 2;
N = Np1 - 1;

% If we are at the coarsest level take 10 GS iterations - should be enough
if ( (N==2) || (M==2))
  Uout = GaussSeidel( Uin, f, dr, dtheta, r, iters ); 
else
% Otherwise begin the cycle from fine to coarsest

% Start by smoothing input with 10 GS iterations
Usmooth = GaussSeidel( Uin, f, dr, dtheta, r, iters ); 

% Compute the residuals on a coarser grid
res  = residual(Usmooth, f, dr, dtheta, r); 

% Restrict the residual to the grid half the size
[res2, dr2, dtheta2, r2] = restrict(res, dr, dtheta, r); 

% Now call this routine to solve the error equation on the next grid
err = MultiGridV( zeros(size(res2)), res2 , dr2, dtheta2, r2); 

% Now interpolate the course error onto finer grid and add to smoothed
Usmooth = Usmooth + interpolate(err); 

% Finally, smooth out any new high-frequency error (post-smoothing) 
Uout = GaussSeidel( Usmooth, f, dr, dtheta, r, iters );

% This completes a Multigrid V-cycle. It can be called again if need be
end

