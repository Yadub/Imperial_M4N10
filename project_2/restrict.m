function [coarse, dr2, dtheta2, r2] = restrict(fine, dr, dtheta, r)
%  restriction routine (full weighting NOT YET)

[Mp2,Np1] = size(fine);
M = Mp2 - 2;
N = Np1 - 1;
M2 = M/2;
N2 = N/2;

coarse = zeros(M2 + 2, N2 + 1); 

% for n = 1:N2
%     coarse(1,n) = fine(1, 2 * n - 1);
% end

coarse(2:M2,2:N2) =  1/4*fine(3:2:M-1, 3:2:N-1) + ...
                     1/8*fine(2:2:M-2, 3:2:N-1) + ...
                     1/8*fine(4:2:M,   3:2:N-1) + ...
                     1/8*fine(3:2:M-1, 2:2:N-2) + ...
                     1/8*fine(3:2:M-1, 4:2:N)   + ...
                     1/16*fine(2:2:M-2,2:2:N-2) + ...
                     1/16*fine(2:2:M-2,4:2:N)   + ...
                     1/16*fine(4:2:M,  2:2:N-2) + ...
                     1/16*fine(4:2:M,  4:2:N);
% coarse(:, N2+1) = coarse(:, 1);
% for m = 2:M2+1
%     for n = 1:N2
%         coarse(m,n) = fine(2*m-1,2*n-1);
%         % Until a better weighting can be found
%     end
%     coarse(m, N2+1) = coarse(m, 1);
% end

% for n = 1:N2
%     coarse(M2+1,n) = fine(M+1, 2*n-1);
% end
% coarse(M2+1,:) = fine(M+1, 2*(1:N2)-1); 
% coarse(M2+1, N2+1) = coarse(M2+1, 1); % Periodic
% coarse(M2+2, :) = coarse(M2, :); % Insulation

dr2 = dr * 2;
dtheta2 = dtheta * 2;
r2 = r(1:2:M+1);
end
