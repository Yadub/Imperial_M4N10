M = 2^6;
N = 2^8;
dt = 0.001;
tstep = 1000;
% Compute errors when varying dr
M_vals = [2^3, 2^4, 2^5, 2^6, 2^7];
dr = 1 ./ M_vals;
i = 1;
for M = M_vals
    [max_errorM(i), average_errorM(i)] = q1Errors(M,N,dt,tstep)
    i = i+1;
end
% Compute errors when varying dtheta
N_vals = [2^4, 2^5, 2^6, 2^7, 2^8];
dtheta = 2 * pi ./ N_vals;
i = 1;
for N = N_vals
    [max_errorN(i), average_errorN(i)] = q1Errors(M,N,dt,tstep)
    i = i+1;
end
% Compute errors when varying dt
M = 2^6;
dt_vals = [0.001, 0.00075, 0.0005, 0.00025, 0.0001, 0.00005];
i = 1;
for dt = dt_vals
    [max_errorT(i), average_errorT(i)] = q1Errors(M,N,dt,tstep)
    i = i+1;
end

% Plots results
figure(); hold on;
pr = polyfit(dr,average_errorM,2);
drM = min(dr):min(dr):max(dr);
estimatedEM = polyval(pr,drM);
plot(dr,average_errorM,'o','DisplayName','Average Errors Results');
plot(drM,estimatedEM,'DisplayName','Computed Quadratic Fit');
title('Q1: Average Error vs dr. dtheta and dt fixed'); xlabel('dr'); ylabel('Average Error'); legend('show');

figure(); hold on;
ptheta = polyfit(dtheta,average_errorN,2);
dthetaN = min(dtheta):min(dtheta):max(dtheta);
estimatedEN = polyval(ptheta,dthetaN);
plot(dtheta,average_errorN,'o','DisplayName','Average Errors Results');
plot(dthetaN,estimatedEN,'DisplayName','Computed Quadratic Fit');
title('Q1: Average Error vs dtheta. dr and dt fixed'); xlabel('dtheta'); ylabel('Average Error'); legend('show');

figure(); hold on;
pt = polyfit(dt_vals,average_errorT,2);
dtT = min(dt_vals):min(dt_vals):max(dt_vals);
estimatedET = polyval(pt,dtT);
plot(dt_vals,average_errorT,'o','DisplayName','Average Errors Results');
plot(dtT,estimatedET,'DisplayName','Computed Quadratic Fit');
title('Q1: Average Error vs dt. dr and dtheta fixed'); xlabel('dt'); ylabel('Average Error'); legend('show');
