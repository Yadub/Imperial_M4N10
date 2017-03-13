b = 5; % Set value of b
M_vals = [25 50 75];
N_vals = M_vals;
n = length(M_vals);
omega = 1.8:0.05:1.95;
No = length(omega);
iters_SOR = zeros(No,n,n);
for k = 1:No 
    w = omega(k)
    for j = 1:n
        M = M_vals(j)
        for i = 1:n
            N = N_vals(i);
            [~, iters_SOR(k, j, i)] = q2TestFunc(b, M, N, w, 1e-8, 2);
        end
    end
end
iters_SOR_R = zeros(No,n,n);
for k = 1:No 
    w = omega(k)
    for j = 1:n
        M = M_vals(j)
        for i = 1:n
            N = N_vals(i);
            [~, iters_SOR_R(k, j, i)] = q2TestFunc(b, M, N, w, 1e-8, 2, true);
        end
    end
end

