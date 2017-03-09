function [u,largest_residual] = FullMG(f, b, f_theta, iPlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full geometric multigrid                                        %
%         Solves delsqr(u)=f on a grid given by size of f         %                                                        %
% (homogeneous Dirichlet boundary conditions)                     %
% (dimension N has to be of the form 2^klevel + 1)                %
%                                                                 %
% f:    right-hand side (n x m)-matrix                            %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    iPlot = false;
end

iters = 10;
[Mp1,Np1] = size(f);
M = Mp1 - 1;
N = Np1 - 1;
 
%   determine maximum number of grid levels
klevel = round( log(min(M,N))/log(2));

%    As we will solve on all grids, we need to calculate the
%    restrictions of the right-hand side to all grid levels
%    ff(1) is finest grid ff(klevel) coarsest
ff = cell(klevel,1);
f = [f;f(M,:)]; % BC: Insulating (Neumann) as Ghost points
ff{1} = f;
for k=2:klevel
  ff{k} = restrict(ff{k-1});
end

%   solve the equation at the coarsest level
f0 = ff{klevel};
%    u0 = initu(zeros(size(f0,1),size(f0,2)));
u0 = zeros(size(f0));
[Mp1,Np1] = size(f0);   % Extract M and N
M0 = Mp1-2;
N0 = Np1-1;

dr = (b - 1) / M0;      % Set delta
dtheta = 2 * pi / N0;   % Set delta theta
r = 1:dr:b;             % Initialize r array

u0 = GaussSeidel( u0, f0, dr, dtheta, r, 10 );

%   loop over all higher levels (coarser grids) using V-cycles
for k = klevel-1:-1:1

% Interpolate solution to next-higher grid level
% u1 = initu(prolong(u0));
    u1 = interpolate(u0);
    f1 = ff{k};
    
    dr = dr / 2;
    dtheta = dtheta / 2;
    r = 1:dr:b;             

    if k == 1 % BC Dirichlet u(1, theta) = f(theta)
        u1(1,:) = f_theta(:);
    end    
%     perform one multigrid V-cycle
    u0 = MultiGridV(u1, f1, dr, dtheta, r);
    
    if iPlot
        % Get polar grid
        theta = 0:dtheta:2*pi;  % Initialize theta array
        [R, Theta] = meshgrid(r, theta);
        xx = R.*cos(Theta);
        yy = R.*sin(Theta);
        M0 = M0*2;
        % Plot solution
        figure(k);
        title(['Estimated Solution. Depth level: ', num2str(k)]);
        surf(xx,yy,u0(1:M0+1,:)')
        xlabel('x','FontSize',18)
        ylabel('y','FontSize',18)
        zlabel('u','FontSize',18)
        grid off
        lighting phong
        camlight headlight
        camlight right
        pause
    end
end

largest_residual = max(max(residual(u0, f, dr, dtheta, r)));

%   output solution
u = u0;
end
