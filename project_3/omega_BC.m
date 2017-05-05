function [ omega ] = omega_BC( omega, psi, dr)
%OMEGA_BC: Implement boundary conditions omega = -psi'', perodic in theta

% Get size of the problem
[Mp1, ~] = size(omega);
M = Mp1 - 1;
% Compute boundary conditions using the neumann condition, psi_r = 0, and
% dirichlet condition, psi = 0, on the boundary with a 2nd order centered difference scheme
% omega(1,:)    = 2 * psi(2,:) / dr / dr;  
% omega(M+1,:)  = 2 * psi(M,:) / dr / dr;  

% A more accurate difference scheme. From Lecture 20.
omega(1  ,:)  = ( - 4 * psi(2,:) + .5 * psi(3  ,:) ) / dr / dr ;  
omega(M+1,:)  = ( - 4 * psi(M,:) + .5 * psi(M-1,:) ) / dr / dr ;    

end

