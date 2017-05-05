function Uout = MultigridV(Uin,f,delta_r,delta_theta,r)

    M = size(f,1); 
    N = size(f,2);  

% if we are at the coarsest level take 10 GS iterations - should be enough
    if ((M==3)||(N==2)) 
      Uout = GS(Uin,f,delta_r,delta_theta,r,10); 
    else
% otherwise begin the cycle from fine to coarsest
%
%   Start by smoothing input with 10 GS iterations
      Usmooth = GS(Uin,f,delta_r,delta_theta,r,10); 
%
%   compute the residuals on a coarser grid
      res  = residual(Usmooth,f,delta_r,delta_theta,r); 

%     and restrict the residual to the grid half the size
      reshalf = restrict(res); 

%     Now call this routine to solve the error equation on the next grid
      err = MultigridV(zeros(size(reshalf)),reshalf,2*delta_r,2*delta_theta,r(1:2:M)); 

%     Now interpolate the course error onto finer grid and add to smoothed
      Usmooth = Usmooth + interpolate(err); 

%     Finally, smooth out any new high-frequency error (post-smoothing) 
      Uout = GS(Usmooth,f,delta_r,delta_theta,r,10);

% This completes a Multigrid V-cycle. If we want, we can call it again
    end

