a = 1; % Set a
b = 2; % Set b
D = 0.01; % Diffusivity constant
Q = 0.01;
tmax = 4;
maxdt = 8000;
k = tmax / maxdt; % k = 0.0005

% Set up h values
N_vals = 50:10:300;
num_h = length(N_vals);
h = (b - a) ./ N_vals;

E_h = zeros( num_h, 1);

for i = 1:num_h
    N = N_vals(i); % Set N
    [u, r, t] = solveq2(a,b,Q,D,N,maxdt,k); % Solve

    % Set point
    pr = N/2 + 1;
    r0 = r(pr);
    t0 = t(end);
    u0 = (b-r0)/(b-a) + (r0-b)*(r0-a)*exp(-t0);
    
    E_h(i) = abs(u0 - u(pr,end)); % Save error
end

plot(h,E_h,'x'); % Plot h vs error
xlabel('h'); ylabel('Error');