function Uout = MultigridV(Uin,f); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multigrid V-cycle                                                %
%                                                                  %
% Calls itself recursively                                         %
% Homogeneous Dirichlet boundary conditions                        %
% squatr grid dimension N has to be of the form 2^k + 1            %
%                                                                  %
% Ain:  guessed solution (n x m)-matrix                            %
% f:    right-hand side (n x m)-matrix                             %
%                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = size(f,1); 
    m = size(f,2);  

% if we are at the coarsest level take 10 GS iterations - should be enough
    if ((n==3)|(m==3)) 
      Uout = GS(Uin,f,10); 
    else
% otherwise begin the cycle from fine to coarsest
%
%   Start by smoothing input with 10 GS iterations
      Usmooth = GS(Uin,f,10); 
%
%   compute the residuals on a coarser grid
      res  = residual(Usmooth,f); 

%     and restrict the residual to the grid half the size
      reshalf = restrict(res); 

%     Now call this routine to solve the error equation on the next grid
      err = MultigridV(zeros(size(reshalf)),reshalf); 

%     Now interpolate the course error onto finer grid and add to smoothed
      Usmooth = Usmooth + interpolate(err); 

%     Finally, smooth out any new high-frequency error (post-smoothing) 
      Uout = GS(Usmooth,f,10);

% This completes a Multigrid V-cycle. If we want, we can call it again
    end

