a = 1; % Set a
b = 2; % Set b
D = 0.01; % Diffusivity constant
Q = 0.01;
tmax = 4;
N = 100;
h = (b - a) / N; % h = 0.01

% Set k values
maxdt = 4000:400:12000;
num_k = length(maxdt);
k = tmax ./ maxdt;

E_k = zeros( num_k, 1); 

for i = 1:num_k
    kval = k(i);
    Nt = maxdt(i);
    [u, r, t] = solveq2(a,b,Q,D,N,Nt,kval); % Solve

    % Set point
    pr = N/2 + 1;
    r0 = r(pr);
    t0 = t(end);
    u0 = (b-r0)/(b-a) + (r0-b)*(r0-a)*exp(-t0);
    
    E_k(i) = abs(u0 - u(pr,end)); % Save error
end

plot(k,E_k,'x'); % Plot k vs error
xlabel('k'); ylabel('Error'); 