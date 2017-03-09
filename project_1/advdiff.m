% Explicit method to solve u_t = V u_x + a u_xx 
% Parameters: -L < x < L, N=number of points
% a = diffusivity, V(x)=velocity
%
L=5;
a=0.01;
Vmax=1.50;
N=600;

maxdt=4000; % Maximum number of timesteps
dt=0.01; % timestep

dx=2*L/N;
time=0;

u=zeros(N+1,maxdt+1);  % Array for u(x,t) initialised to zero
x=zeros(1,N+1);
V=zeros(1,N+1);

for n=1:(N+1)        % Set up grid x(n) and initialise u(n,1)
    x(n)=-L+(n-1)*dx;
    V(n)=Vmax;
    if(abs(x(n))<=0.5)
        u(n,1)=1;
    else
        u(n,1)=0;
    end
end
u(1,1)=0;
u(N+1,1)=0;

time(1)=0;
for j=1:maxdt    % Step in time
    time(j+1)=time(j)+dt;
    for n=2:N
        u(n,j+1)=u(n,j)+dt*(V(n)*(u(n+1,j)-u(n-1,j))/(2*dx) + a*(u(n+1,j)-2*u(n,j)+u(n-1,j))/dx^2);
    end;
    u(1,j+1)=0;
    u(N+1,j+1)=0;
end;

figure(1);    % Contour u(x,t)
contour(time,x,u,30);
%xlim([0,10])
%ylim([-10,0])
colorbar
hold off
figure;    % u(x) for various t in frame moving with Vmax
hold off
for j=1:200:maxdt
    hold on
    plot(x+Vmax*time(j),u(:,j));
    xlim([-3,3])
end;
hold off
