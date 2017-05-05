function [ Tx ] = annulusBuoyancy( T, dr, dtheta, r, Pr, Ra )
%ANNULUSBUOYANCY: Computes buoyancy term, essentially T_theta.
% Periodic in theta

% Extract N
[~,N] = size(T);  
% Arrays of interest
a0toNm1 = [N 1:N-1];
a2toNp1 = [2:N 1];
% Horizontal temperature gradient (i.e. in the theta direction)
Tx = - Ra * Pr * ( T(:,a2toNp1) - T(:,a0toNm1) ) / (2 * dtheta) ./ r(:);

end

