function fine = mg_interpolate( coarse )
% Interpolation routine (linear weighting) 
% Editted code taken from Blackboard done by Dr Mestel

% Extract M and N and set their double values
[Mp1,N] = size(coarse);   
M = Mp1 - 1;
M2 = M * 2;
N2 = N * 2;
% Arrays of interest
a2toNp1 = mod(1:N,N) + 1;
% Initalize and compute values
fine = zeros(M2 + 1, N2);
fine(1:2:M2+1, 1:2:N2) =  coarse(1:M+1,1:N);
fine(3:2:M2-1,2:2:N2)    = (coarse(2:M,1:N) + coarse(2:M,a2toNp1))/2;
fine(2:2:M2  ,3:2:N2-1)  = (coarse(1:M,2:N) + coarse(2:M+1,2:N))/2;
fine(2:2:M2 ,2:2:N2)     = (coarse(1:M,1:N) + coarse(2:M+1,a2toNp1) ...
                          + coarse(2:M+1,1:N) + coarse(1:M,a2toNp1)) / 4;        

end
