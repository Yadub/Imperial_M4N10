     function u = FullMG(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full geometric multigrid                                        %
%         Solves delsqr(u)=f on a grid given by size of f         %                                                        %
% (homogeneous Dirichlet boundary conditions)                     %
% (dimension N has to be of the form 2^klevel + 1)                %
%                                                                 %
% f:    right-hand side (n x m)-matrix                            %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n=size(f,1);
    m=size(f,2);
 
%   determine maximum number of grid levels
    klevel=round(log(min(n,m)-1)/log(2));

%    As we will solve on all grids, we need to calculate the
%    restrictions of the right-hand side to all grid levels
%    ff(1) is finest grid ff(klevel) coarsest
    ff = cell(klevel,1);
    ff{1} = f;
    for k=2:klevel
      ff{k} = restrict(ff{k-1});
    end

%   solve the equation at the coarsest level
    f0 = ff{klevel};
%    u0 = initu(zeros(size(f0,1),size(f0,2)));
    u0=zeros(size(f0,1),size(f0,2));
    u0 = GS(u0,f0,10);

%   loop over all higher levels (coarser grids) using V-cycles
    for k=klevel-1:-1:1

%     interpolate solution to next-higher grid level
%      u1 = initu(prolong(u0));
      u1=interpolate(u0);
      f1 = ff{k};

%     perform one multigrid V-cycle
      u0=MultigridV(u1,f1);

 figure(k)
 [nn,mm] = size(u0);
 xx = linspace(0,1,nn);
 yy = linspace(0,1,mm);
surf(xx,yy,u0)
xlabel('x','FontSize',18)
ylabel('y','FontSize',18)
zlabel('u','FontSize',18)
grid off
lighting phong
camlight headlight
camlight right
% pause

    end

%   output solution
    u=u0;

