function [coarse, dr2, dtheta2, r2] = mg_restrict(fine, dr, dtheta, r)
% Restriction routine (linear weighting) 
% Editted code taken from Blackboard done by Dr Mestel

% Extract M and N and set their half values
[Mp1,N] = size(fine); 
M = Mp1-1;
M2 = M/2;
N2 = N/2;
% Initalize coarse grid and compute values
coarse = zeros(M2+1,N2);  
coarse(2:M2,1:N2) =  1/4*fine(3:2:M-1,  1:2:N-1) + ...
                      1/8*fine(2:2:M-2, 1:2:N-1) + ...
                      1/8*fine(4:2:M,   1:2:N-1) + ...
                      1/8*fine(3:2:M-1, mod((0:2:N-2)-1,N)+1) + ...
                      1/8*fine(3:2:M-1, 2:2:N)   + ...
                      1/16*fine(2:2:M-2,mod((0:2:N-2)-1,N)+1) + ...
                      1/16*fine(2:2:M-2,2:2:N)   + ...
                      1/16*fine(4:2:M,  mod((0:2:N-2)-1,N)+1) + ...
                      1/16*fine(4:2:M,  2:2:N);
% Get new values for dr, dtheta, and the r vector
if nargout > 1
    dr2 = dr * 2;
    dtheta2 = dtheta * 2;
    r2 = r(1:2:M+1);
end

end
