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
omega = 1.5:0.005:1.95;
omega = [1:0.05:1.5 omega];
No = length(omega);
T = 4000;
R = 0.2:0.1:1.5;
R = [R 1.1:0.1:2];
N_vals = sqrt(T./R);
M_vals = T./N_vals;
N_vals = round(N_vals); 
M_vals = round(M_vals);
Nmn = length(N_vals);

T = zeros(No,Nmn);
I = zeros(No,Nmn);
p = 2; % Change this value of for the different plots

for j = 1:No
    w = omega(j);
    for i = 1:Nmn
        M = M_vals(i);                % Set M
        N = N_vals(i);                % Set N
        [~, iters, timetaken] = q2TestFunc( b, M, N, w, tol, p);
        T(j,i) = timetaken;
        I(j,i) = iters;
    end
    w
end
ratio = M_vals./N_vals;
figure(); hold on;
for i = 1:Nmn 
    plot(omega, I,'DisplayName',['Ratio = ', num2str(ratio(i))]);
end
xlabel('Ratio (M/N)'); ylabel('Iterations');
title('Varying omega for M/N with T = 4000');
legend('show');
