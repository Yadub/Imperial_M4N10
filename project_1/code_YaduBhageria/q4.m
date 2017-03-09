figure();
hold on;

D = 0.01;
Qs = [0, 0.25*D, 0.5*D, D, 10*D, 25*D, 64*D, 70*D, 100*D];
num_Q = length(Qs);

for i = 1:num_Q
    Q = Qs(i);
    [u, r] = solveq4(Q, D);
    qod = Q/D;
    plot(r, u,'DisplayName',['Q/D = ',num2str(qod)]);
end

hold off;
xlim([1,16]);
xlabel('r'); ylabel('u_i_n_f(r)'); legend('show');
title('Plot of u_i_n_f(r) for b/a = 17')