function Uout = MultiGridV(Uin, f, dr, dtheta, r)
% Performs 1 multigrid V cycle for guessed solution Uin with boundary
% function f

iters = 10;

[Mp2,Np1] = size(Uin);
M = Mp2 - 2;
N = Np1 - 1;

% theta = 0:dtheta:2*pi;  % Initialize theta array

% if we are at the coarsest level take 10 GS iterations - should be enough
if ( (N==2) || (M==2))
  Uout = GaussSeidel( Uin, f, dr, dtheta, r, iters ); 
else
% otherwise begin the cycle from fine to coarsest
%
%   Start by smoothing input with 10 GS iterations
Usmooth = GaussSeidel( Uin, f, dr, dtheta, r, iters ); 

%   compute the residuals on a coarser grid
res  = residual(Usmooth, f, dr, dtheta, r); 

%     and restrict the residual to the grid half the size
[res2, dr2, dtheta2, r2] = restrict(res, dr, dtheta, r); 

%     Now call this routine to solve the error equation on the next grid
err = MultiGridV( zeros(size(res2)), res2 , dr2, dtheta2, r2); 

%     Now interpolate the course error onto finer grid and add to smoothed
Usmooth = Usmooth + interpolate(err); 

%     Finally, smooth out any new high-frequency error (post-smoothing) 
Uout = GaussSeidel( Usmooth, f, dr, dtheta, r, iters );

% This completes a Multigrid V-cycle. If we want, we can call it again
end

