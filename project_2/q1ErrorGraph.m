% Yadu Bhageria
% 00733164
b = 5;                 % Chosen value of b
tol = 1e-8;             % Set tolerance

% Fixed M - Error just decreases with more points overall
% N_vals = [2^6, 2^6, 2^6, 2^6, 2^6];
% M_vals = [2^3, 2^4, 2^5, 2^6, 2^7];
% Fixes N - Error just decreases with more points overall
% M_vals = [2^6, 2^6, 2^6, 2^6, 2^6];
% N_vals = [2^3, 2^4, 2^5, 2^6, 2^7];
% ~Fixed MxN
% N_vals = [2^8, 164, 2^7, 102, 2^6, 2^5];
% M_vals = [2^4, 25,  2^5, 40, 2^6, 2^7];

% Fix Total points and ratios to be computed
T = 4000;
R = 0.05:0.05:1;
R = [R 1.1:0.1:2];
N_vals = sqrt(T./R);
M_vals = T./N_vals;
N_vals = round(N_vals); 
M_vals = round(M_vals);

E_vals = [];
R_vals = [];

p = 2; % Change this value of for the different plots

for i = 1:length(N_vals)
    M = M_vals(i);                % Set M
    N = N_vals(i);                % Set N

    [avgError, iters, timetaken] = q1TestFunc( b, M, N, tol, p);
    E_vals = [E_vals, avgError];
end

ratio = M_vals./N_vals;
figure();
plot(ratio,E_vals,'x-','DisplayName','Largest Error');
xlabel('Ratio (M/N)'); ylabel('Error');
title('Varying M/N for T = 4000');
Best_ratio = min(ratio);
