b = 5; % Set value of b

M_vals = [16 32 64 128];
N_vals = M_vals;
n = length(M_vals);
Errors = zeros(n);
iters = zeros(n);
timetaken = zeros(n);

for i = 1:n
    M = M_vals(i);
    for j = 1:n
        N = N_vals(j);
        [Errors(i,j), iters(i,j), timetaken(i,j)] = q1TestFunc(b, M, N);
    end
end

display(Errors);
display(iters);
display(timetaken);