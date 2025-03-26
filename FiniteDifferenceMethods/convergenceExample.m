% Define the parameters
L = 1;             % Length of the rod
T = 1;             % Total time
N = 100;           % Number of spatial grid points
alpha = 0.02;      % Thermal diffusivity constant

for k=1:5
    % Refine spacial resolution
    if k > 1
        N = N*2;
    end

    dx = L/N;
    dt = 0.5*dx*dx/alpha;

    M = int32(T/dt);

    x = linspace(0,L,N+1);

    % Initialize the temperature matrix
    u = zeros(N+1,M+1);

    % Set the initial condition
    u(:,1) = sin(pi*x);
    
    time = 0.0;
    toatlError(:) = 0;

    for j = 1:M
        for i = 2:N
            u(i,j+1) = u(i,j) + alpha*dt/dx^2 *(u(i+1,j) - 2*u(i,j) + u(i-1,j));
        end
    
        time = time + dt;
        exact = sin(pi*x)*exp(-alpha*(pi)^2*time);
    
        Error = exact' - u(:,j);
        totalError(j) = sum(Error);
        
    end
    % Error
    E(k) = norm(totalError,inf);
    DX(k) = dx;
end

%%
figure
plot(DX,E,'o-')
xlabel("Delta x")
ylabel("Error")

