function [ ur, utheta, psi ] = annulusUinit( M, N, dr, dtheta, r )
%TINIT A velocity field for U

% Compute b and b/2
b = 1 + dr * M;
rc = b / 2;
% Initalize psi
psi = zeros(M+1, N);
% Compute psi
for m = 2:M
    rm = r(m);
    for n = 1:N
        theta = dtheta * (n-1);
        psi(m, n) = cos(theta)^2 * (rm - b)^2 * (rm - 1)^2 * exp(- rm) ;
%         psi(m, n) =  exp(-sqrt( (rm - rc)^2 + (theta-pi)^2 )) * sin(theta) * (rm-b) * (rm-1);
%         psi(m, n) =  - exp( -(rm - rc)^2 - (theta - pi)^2 ) * (rm - b) * (rm - 1);
    end
end
% Normalize psi
psi = psi / max(max(psi));

% Compute the velocities
[ur, utheta] = annulusVelocity(psi, dr, dtheta, r);

end

