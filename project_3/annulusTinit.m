function [ T ] = annulusTinit( T, dr, dtheta, q )
%TINIT A small r-dependent perturbation of T to kickstart instabilities

if nargin < 4, q = 1; end

[Mp1, N] = size(T);
M = Mp1 - 1;

b = 1 + dr * M;
rc = 1 + dr * M/2;

for i = 1:M
    r = 1 + dr * (i-1);
    for j = 1:N
        theta = dtheta * (j-1);
        if q == 1, T(i,j) = .01 * sin(pi*r) * sin(theta); end
        if q == 2, T(i,j) = exp( - sqrt( (r - rc)^2 + (theta - pi).^2) ); end
%         if q == 2, T(i,j) =  sin(theta) * exp( - (r - b) * (r - 1) ); end
    end
end
if q ~= 1, T(1,:) = 1; end

end

