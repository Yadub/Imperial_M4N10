function Tout = T_MultiGridV(Tin, f, dt, dr, dtheta, r, k)
% Multigrid V-cycle for Diffusion equation  
% Performs 1 multigrid V cycle for the Poisson equation in an annulus,
% D^2 u = f, given a starting guess Tin for delta r, delta theta, vector r
% ( M has to be of the form 2^k + 1)                                                                       
%       Tin:  guessed solution (M+1 x N)-matrix                        
%       f:    right-hand side  (M+1 x N)-matrix                        
%       dt:   time-step
%       dr:   r-step
%       dtheta: theta-step
%       r:    r values (M+1)-vector
% No need to consider boundary conditions seperately for T and residuals as
% a simple boundary of T = 1 at r = 1 and T = 0 at r = b is present

% Editted code taken from Blackboard done by Dr Mestel

iters = 10; % Number of Gauss Seidel Smoothing iterations

[Mp1,N] = size(Tin);  % Extract M and N
M = Mp1 - 1;

% If we are at the coarsest level take 10 GS iterations - should be enough
if ( (N==2) || (M==2))
    Tout = T_GaussSeidel( Tin, f, dt, dr, dtheta, r, iters, k); 
else
    % Otherwise begin the cycle from fine to coarsest

    % Start by smoothing input with 10 GS iterations
    Tsmooth = T_GaussSeidel( Tin, f, dt, dr, dtheta, r, iters, k); 

    % Compute the residuals on a coarser grid
    res  = T_residual(Tsmooth, f, dt, dr, dtheta, r, k); 

    % Restrict the residual to the grid half the size
    [res2, dr2, dtheta2, r2] = mg_restrict(res, dr, dtheta, r); 

    % Now call this routine to solve the error equation on the next grid
    err = T_MultiGridV( zeros(size(res2)), res2 , dt, dr2, dtheta2, r2, k); 

    % Now interpolate the course error onto finer grid and add to smoothed
    Tsmooth = Tsmooth + mg_interpolate(err); 

    % Finally, smooth out any new high-frequency error (post-smoothing) 
    Tout = T_GaussSeidel( Tsmooth, f, dt, dr, dtheta, r, iters, k); 
end
% This completes a Multigrid V-cycle. It can be called again if need be
end

