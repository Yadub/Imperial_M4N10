function [] = plotAnnulusVectorField( ur, utheta, r, theta, fig )
%PLOTANNULUSVECTORFIELD: 
%       Plots a vector field in the annulus for give r-theta
%           components of the field with grid vector r and theta. 
%           Figure to plot on can be specified in fig

[Mp1, N] = size(ur);
M = Mp1 - 1;

if nargin < 5
    fig = figure();
end
% Set figure
figure(fig);
% Get x-y coordinates
[R, Theta] = meshgrid(r, theta);
xx = R.*cos(Theta);
yy = R.*sin(Theta);
size(ur)
size(Theta)
% Compute x-y Vector Field components
rper = zeros(M+1,N+1);
rper(:,1:N) = ur;
rper(:,N+1) = ur(:,1);
tper = zeros(M+1,N+1);
tper(:,1:N) = utheta;
tper(:,N+1) = utheta(:,1);
xper = rper .* cos(Theta)' + tper .* cos(Theta + pi/2)';
yper = rper .* sin(Theta)' + tper .* sin(Theta + pi/2)';
% Plot components
quiver(xx',yy',xper,yper,1);
% Label
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
drawnow

end

