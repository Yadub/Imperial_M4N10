function [] = plotAnnulusScalarField( Q, xx, yy, fig )
%PLOTANNULUSSCALERFIELD Plots the scalar field Q for position matrices xx and yy
% onto fig

if nargin < 4, fig = figure(); end

[Mp1, N] = size(Q);

figure(fig) 
Qplot = zeros(Mp1,N+1);
Qplot(:,1:N) = Q(:,1:N);
Qplot(:,N+1) = Q(:,1);
contourf(xx',yy',Qplot,10,'r')
axis image
xlabel('x','FontSize',18)
ylabel('y','FontSize',18);
colorbar
drawnow

end

