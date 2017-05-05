% Domain and grid size
b =  2;
M = 2^2;                % Set M
N = 2^4;                % Set N
dr = (b - 1) / M;       % Set delta r
dtheta = 2 * pi / N;    % Set delta theta
r = 1:dr:b;             % Initialize r array
theta = 0:dtheta:2*pi;  % Initialize theta array

% For plotting
% Get x-y grid
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
figure();

plot(xx',yy','.k','MarkerSize',20)
% surf(xx',yy',Qplot)
axis image
xlabel('x','FontSize',18)
ylabel('y','FontSize',18);
title(['Annular Grid. M = ', num2str(M),  ', N = ', num2str(N)])