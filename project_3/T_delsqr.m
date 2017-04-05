function [ delsqrT ] = T_delsqr( T, dr, dtheta, r )
%T_DELSQR: Computes Laplacian for T

% Extract size of grid
[Mp1,N] = size(T);
M = Mp1 - 1;
% Initalize array for storage
delsqrT = zeros(M+1, N);
% Compute values of interest
rr = 1 / dr / dr;
rtheta = 1 / dtheta / dtheta;
% Indice arrays of interest     
a0toNm1 = [N 1:N-1];
a2toNp1 = mod(1:N,N) + 1;
% Compute laplacian - 
for i=2:M
    delsqrT(i, :) = rr * ( T(i+1, :) - 2 * T(i, :) + T(i-1, :) ) + .5 / ...
                    r(i) / dr * (T(i+1, :) - T(i - 1, :)) + rtheta / r(i) / ...
                    r(i) * ( T(i, a2toNp1) - 2 * T(i, :) + T(i, a0toNm1) );
end

end

