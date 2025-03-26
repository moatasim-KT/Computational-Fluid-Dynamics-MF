% Define the parameters
L = 2*pi;          % Length
T = 3;             % Total time
N = 100;           % Number of spatial grid points
a = 1;             % Wave speed

% Spatial discretization
dx = L/N;

% Stability condition
dt = abs(dx/a);

% Number of iterations
M = int32(T/dt);

% Spatial grid
x = linspace(0,L,N+1);

% Initialize the matrix
u = zeros(N+1,1);

% Set the initial condition
u(:,1) = sin(x);
%%

u_old = u;
u_older = u;

% Apply the finite difference method
figure
for i = 1:500
    for  j=2:N
        u(j) = 2*(1 - a^2)*u_old(j)  + a^2*(u_old(j+1) + u_old(j-1)) - u_older(j);
    end

    u_older = u_old;
    u_old = u;

    plot(x,u)
    xlabel("Position")
    ylabel("U")
    ylim([-1 1])
    xlim([0 2*pi])
    set(gcf,'Position',[400 400 800 380])

    exportgraphics(gcf,'two-way-wave.gif','Append',true)
end