function omegaOut = omega_MultiGridV(omegaIn, f, psi, dt, dr, dtheta, r, Pr, Mmax)
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
% Boundary conditions considered for the full problem using omega_BC()

% Editted code taken from Blackboard done by Dr Mestel

iters = 10; % Number of Gauss Seidel Smoothing iterations

[Mp1,N] = size(omegaIn);  % Extract M and N
M = Mp1 - 1;

% Set boundary condition bolean to true if solving the full problem
if Mmax == M
    bc = true;
else
    bc = false;
end

% If we are at the coarsest level take 10 GS iterations - should be enough
if ( (N==2) || (M==2))
    omegaOut = omega_GaussSeidel( omegaIn, f, psi, dt, dr, dtheta, r, iters, Pr, bc); 
else
    % Otherwise begin the cycle from fine to coarsest

    % Start by smoothing input with 10 GS iterations
    omegaSmooth = omega_GaussSeidel( omegaIn, f, psi, dt, dr, dtheta, r, iters, Pr, bc); 
    
    % Compute the residuals on a coarser grid - No need to implement BC
    % condition thats already been done in the previous step within GS
    res  = omega_residual(omegaSmooth, f, dt, dr, dtheta, r, Pr); 

    % Restrict the residual to the grid half the size
    [res2, dr2, dtheta2, r2] = mg_restrict(res, dr, dtheta, r); 
    % Restrict the stream function to half the grid size
    psi2 = mg_restrict(psi);
    
    % Now call this routine to solve the error equation on the next grid
    err = omega_MultiGridV( zeros(size(res2)), res2 , psi2, dt, dr2, dtheta2, r2, Pr, Mmax); 

    % Now interpolate the course error onto finer grid and add to smoothed
    omegaSmooth = omegaSmooth + mg_interpolate(err); 

    % Finally, smooth out any new high-frequency error (post-smoothing) 
    omegaOut = omega_GaussSeidel( omegaSmooth, f, psi, dt, dr, dtheta, r, iters, Pr, bc); 
end
% This completes a Multigrid V-cycle. It can be called again if need be
end

