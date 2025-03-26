% Define the parameters
L = 2*pi;          % Length
T = 3;             % Total time
N = 100;           % Number of spatial grid points
a = 1;             % Wave speed

% Spatial discretization
dx = L/N;

% Stability condition
dt = abs(dx/a)*0.5;

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
for i = 1:M
    % Write your finite difference scheme here
    % Don't forget about the boundary conditions
    % And carrying three time levels of the solution





    
end

figure
plot(x,u)
xlabel("Position")
ylabel("U")
ylim([-1 1])
xlim([0 2*pi])
set(gcf,'Position',[400 400 800 380])