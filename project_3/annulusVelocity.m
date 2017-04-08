function [ ur, utheta ] = annulusVelocity( psi, dr, dtheta, r )
%ANNULUSVELOCITY Computes the velocity in the r and utheta directions using psi

% Extract size of the problem
[Mp1, N] = size(psi);
M = Mp1 - 1;
% Arrays for periodic theta
a0toNm1 = [N 1:N-1];
a2toNp1 = [2:N 1];
% Initalize arrays
ur = zeros( size(psi) );
utheta = zeros( size(psi) );
% Compute the velocities in the r and theta directions
utheta(2:M,:) = - ( psi(3:M+1,:) - psi(1:M-1,:) ) / (2 * dr);
ur(2:M,:) = ( psi(2:M,a2toNp1) - psi(2:M,a0toNm1) ) ./ (2 * r(2:M)' * dtheta);

end