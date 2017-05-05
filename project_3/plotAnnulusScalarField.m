function [] = plotAnnulusScalarField( Q, xx, yy, fig )
%PLOTANNULUSSCALERFIELD Plots the scalar field Q for position matrices xx and yy
% onto fig

% Set figure if specified
if nargin == 4, figure(fig); end

[Mp1, N] = size(Q);

Qplot = zeros(Mp1,N+1);
Qplot(:,1:N) = Q(:,1:N);
Qplot(:,N+1) = Q(:,1);
contourf(xx',yy',Qplot,10,'r')
% surf(xx',yy',Qplot)
axis image
xlabel('x','FontSize',18)
ylabel('y','FontSize',18);
colorbar
drawnow

end

