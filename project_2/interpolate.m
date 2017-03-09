function fine = interpolate( coarse )
% Interpolation routine (linear weighting) 
% Editted code taken from Blackboard done by Dr Mestel

[Mp2,Np1] = size(coarse);   % Extract M and N and set their double values
M = Mp2 - 2;
N = Np1 - 1;
M2 = M*2;
N2 = N*2;

fine = zeros(M2 + 2,N2 + 1);% Initalize

fine(1:2:M2+1, 1:2:N2+1) =  coarse(1:M+1,1:N+1);
fine(3:2:M2-1,2:2:N2)    = (coarse(2:M,1:N) + coarse(2:M,2:N+1))/2;
fine(2:2:M2  ,3:2:N2-1)  = (coarse(1:M,2:N) + coarse(2:M+1,2:N))/2;
fine(2:2:M2 ,2:2:N2)     = (coarse(1:M,1:N) + coarse(2:M+1,2:N+1) ...
                          + coarse(2:M+1,1:N) + coarse(1:M,2:N+1)) / 4;
fine(:, N2+1) = fine(:, 1); % BC: Periodic
fine(M+2, :) = fine(M, :);  % BC: Insulating (Neumann)

end
